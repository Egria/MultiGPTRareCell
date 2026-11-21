"""Queue lifecycle, real scientific stages, HTTP authorization and bounded retries."""
import json
import threading
import time
from urllib import error, request

import pytest

from rarecell import job_queue as jobs
from rarecell.config import config
from rarecell.dashboard import Dashboard, detail
from rarecell.prepare import prepare
from rarecell.run import discover
from rarecell.synthetic import generate


def test_pipeline_plan_and_snapshot_config(tmp_path):
    root = tmp_path / "tasks"
    jobs.initialize(root)
    dataset = generate(tmp_path / "data", 1009)
    tid = jobs.submit(root, tmp_path, {"kind": "pipeline", "mode": "A3", "dataset": str(dataset),
                      "model": "discovery-deployment", "annotation_model": "annotation-deployment",
                      "gt": "data/evaluation_fixture/gt.tsv"})
    task = jobs.get_task(root, tid)
    steps = list(jobs.steps(root, task))
    assert [s[0] for s in steps] == ["prepare", "annotate", "discover", "evaluate"]
    assert "--gt" not in sum(steps[:-1], [])
    assert config(root / tid / "settings.yaml")["provider"]["model"] == "discovery-deployment"
    assert task["spec"]["dataset"] == str(dataset)
    jobs.cancel(root, tid)
    jobs.worker(root)
    assert jobs.get_task(root, tid)["status"] == "stopped"
    assert not (root / tid / "run").exists()
    for form in [{"kind": "resume", "run": str(tmp_path)}, {"kind": "pipeline"},
                 {"kind": "pipeline", "dataset": str(dataset)}, {"kind": "demo", "param": []}]:
        with pytest.raises(ValueError):
            jobs.submit(root, tmp_path, form)


def test_resume_existing_checkpoint_and_evaluate(tmp_path):
    root = tmp_path / "tasks"
    jobs.initialize(root)
    dataset = generate(tmp_path / "data", 1009)
    cfg = config(override={"mode": "A2", "prepare": {"views": ["V_LINEAR"], "graph_backend": "igraph"},
                           "provider": {"kind": "scripted_demo"}})
    snap, run = tmp_path / "snapshot", tmp_path / "run"
    prepare(dataset, snap, cfg)
    assert discover(snap, run, tmp_path / "cache", cfg, policy="scripted_demo", pause_after=1)["detail"] == "paused"
    tid = jobs.submit(root, tmp_path, {"kind": "resume", "run": str(run), "gt": "data/evaluation_fixture/gt.tsv"})
    with pytest.raises(ValueError, match="已有排队"):
        jobs.submit(root, tmp_path, {"kind": "resume", "run": str(run)})
    assert "--retry-unknown" not in next(jobs.steps(root, jobs.get_task(root, tid)))
    jobs.worker(root)
    result = detail(root, tid)
    assert result["status"] == "succeeded" and not result["resumable"]
    assert {"metrics.json", "run_summary.json", "partition.tsv"} <= set(result["files"])
    assert result["metrics"]["n_cells"] == 160
    assert result["events"]


def test_http_restart_cancel_and_download(tmp_path):
    root = tmp_path / "tasks"
    with Dashboard(root, tmp_path, 0) as server:
        base = "http://127.0.0.1:%d" % server.server_port
        def call(path, body=None, key=True):
            headers = {"Content-Type": "application/json"}
            if key:
                headers["Authorization"] = "Bearer " + server.token
            data = None if body is None else json.dumps(body).encode()
            req = request.Request(base + path, data=data, headers=headers)
            with request.urlopen(req, timeout=10) as response:
                return response.read()

        # This record predates startup, and its former worker no longer holds the lease.
        old = jobs.submit(root, tmp_path, {"kind": "demo"})
        with jobs.database(root) as db:
            db.execute("UPDATE tasks SET status='running' WHERE id=?", (old,))
        thr = threading.Thread(target=server.serve_forever, kwargs={"poll_interval": .1}, daemon=True)
        thr.start()
        try:
            assert b'RareCell' in call("/", key=False)
            with pytest.raises(error.HTTPError) as missing:
                call("/api/tasks", {"kind": "demo"}, key=False)
            assert missing.value.code == 401
            tid = json.loads(call("/api/tasks", {"kind":"demo", "name":"HTTP demo"}))["id"]
            deadline = time.monotonic() + 60
            while time.monotonic() < deadline:
                data = json.loads(call("/api/tasks/" + tid))
                if data["status"] not in jobs.ACTIVE:
                    break
                time.sleep(.15)
            assert data["status"] == "succeeded", data["log"]
            assert jobs.get_task(root, old)["status"] == "interrupted"
            assert len(call("/api/tasks/" + tid + "/files/partition.tsv").splitlines()) == 161
            with pytest.raises(error.HTTPError):
                call("/api/tasks/" + tid + "/files/tasks.sqlite")
            # Stop a real child after its stage has been claimed.
            stop_id = json.loads(call("/api/tasks", {"kind":"demo"}))["id"]
            deadline = time.monotonic() + 30
            while not jobs.get_task(root, stop_id)["stage"] and time.monotonic() < deadline:
                time.sleep(.02)
            call("/api/tasks/" + stop_id + "/stop", {})
            deadline = time.monotonic() + 15
            while jobs.get_task(root, stop_id)["status"] in jobs.ACTIVE and time.monotonic() < deadline:
                time.sleep(.1)
            assert jobs.get_task(root, stop_id)["status"] == "stopped"
        finally:
            server.shutdown()
            thr.join()
            if server.worker is not None:
                server.worker.wait(timeout=30)
    with Dashboard(root, tmp_path, 0) as restarted:
        assert jobs.get_task(restarted.root, tid)["status"] == "succeeded"


def test_pipeline_reports_missing_live_credential_without_substitution(tmp_path, monkeypatch):
    root = tmp_path / "queue"
    jobs.initialize(root)
    dataset = generate(tmp_path / "input", 1009)
    monkeypatch.delenv("RARECELL_TEST_DASHBOARD_MISSING_KEY", raising=False)
    (tmp_path / "settings.yaml").write_text(
        "prepare: {views: [V_LINEAR], graph_backend: igraph}\n"
        "provider: {api_key_env: RARECELL_TEST_DASHBOARD_MISSING_KEY}\n")
    tid = jobs.submit(root, tmp_path, {"kind": "pipeline", "mode": "A2", "dataset": str(dataset),
                      "config": "settings.yaml", "model": "test-deployment"})
    jobs.worker(root)
    task = detail(root, tid)
    assert task["status"] == "partial" and task["detail"] == "provider_blocked"
    assert task["resumable"] and task["metrics"]["provider_dispatches"] == 0
    assert "manifest.json" in task["files"] and not (root / tid / "annotations").exists()
