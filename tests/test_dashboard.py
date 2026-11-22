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
            with pytest.raises(error.HTTPError) as unauthorized:
                call("/api/model", {"model": "global"}, key=False)
            assert unauthorized.value.code == 401
            assert json.loads(call("/api/model", {"model": "global"}))["model"] == "global"
            assert json.loads(call("/api/tasks"))["model"] == "global"
            with pytest.raises(error.HTTPError) as invalid:
                call("/api/model", {"model": ["invalid"]})
            assert invalid.value.code == 400
            tid = json.loads(call("/api/tasks", {"kind":"demo", "name":"HTTP demo"}))["id"]
            deadline = time.monotonic() + 60
            while time.monotonic() < deadline:
                data = json.loads(call("/api/tasks/" + tid))
                if data["status"] not in jobs.ACTIVE:
                    break
                time.sleep(.15)
            assert data["status"] == "succeeded", data["log"]
            diag_path = "/api/tasks/" + tid + "/diagnostics"
            with pytest.raises(error.HTTPError) as private:
                call(diag_path + "/turn/1", key=False)
            assert private.value.code == 401
            diag = json.loads(call(diag_path))
            assert diag["trace"]["turns"][0]["actions"][0]["action"] == "list_candidates"
            assert diag["summary"]["submitted"] and diag["evaluation"]["targets"]
            assert "raw" not in diag["trace"]["turns"][0]
            turn = json.loads(call(diag_path + "/turn/1"))
            assert turn["commands"][0]["action"] == "list_candidates" and turn["message"]["raw"]
            assert turn["receipts"][0]["result"]["candidates"]
            with pytest.raises(error.HTTPError) as private:
                call(diag_path + "/transcript", key=False)
            assert private.value.code == 401
            text = call(diag_path + "/transcript").decode()
            assert "T0001" in text and "策略原始输出" in text and "已保存决策" in text and "submit" in text
            linked = json.loads(call("/api/tasks", {"kind":"inspect", "run":data["run_path"]}))["id"]
            assert jobs.get_task(root, linked)["status"] == "linked"
            assert json.loads(call("/api/tasks", {"kind":"inspect", "run":data["run_path"]}))["id"] == linked
            linked_diag = json.loads(call("/api/tasks/" + linked + "/diagnostics"))
            assert linked_diag["trace"] == diag["trace"]
            assert linked_diag["evaluation"] == diag["evaluation"]
            with pytest.raises(error.HTTPError):
                call(diag_path + "/turn/999")
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


def test_diagnostics_keep_recorded_rejections_pages_and_fallback(tmp_path, monkeypatch):
    from rarecell.dashboard_diagnostics import diagnostics
    from rarecell.evaluation import evaluate

    queue, run = tmp_path / "tasks", tmp_path / "run"
    jobs.initialize(queue)
    dataset = generate(tmp_path / "fixture", 1009)
    cfg = config(override={"mode": "A2", "prepare": {"views": ["V_LINEAR"], "graph_backend": "igraph"},
                           "provider": {"kind": "scripted_demo"}, "budget": {"stall_turns": 13, "max_agent_turns": 16}})
    snapshot = tmp_path / "snapshot"
    prepare(dataset, snapshot, cfg)
    raw = '{"action":"expand_candidate","unknown_option":1}'
    monkeypatch.setattr("rarecell.run.offline_message", lambda *args: raw)
    assert discover(snapshot, run, tmp_path / "cache", cfg, policy="scripted_demo")["detail"] == "stalled"
    evaluate(run, tmp_path / "fixture/evaluation_fixture/gt.tsv", tmp_path / "evaluation")
    tid = jobs.submit(queue, tmp_path, {"kind":"inspect", "run":str(run)})
    task = jobs.get_task(queue, tid)
    last = diagnostics(queue, task)
    assert last["evaluation"]["source"] == "fallback"
    assert last["evaluation"]["frozen"]["submitted"] is False
    assert last["trace"]["stall"] == last["stall_limit"] == 13
    assert last["trace"]["page"] == 2 and last["trace"]["turns"][0]["number"] == 13
    first = diagnostics(queue, task, page=1)["trace"]["turns"]
    assert [t["number"] for t in first] == list(range(1, 13))
    rejected = first[0]["actions"][0]
    assert rejected["action"] is None and not rejected["changed"] and "unknown_option" in rejected["error"]
    assert diagnostics(queue, task, turn=1)["message"]["raw"] == raw
    from rarecell.dashboard_transcript import transcript
    tag, body = transcript(queue, task)
    assert body.count("── T".encode()) == 13 and b"T0013" in body and b"unknown_option" in body
    assert transcript(queue, task, tag) == (tag, None)
    with jobs.database(queue) as db:
        assert db.execute("SELECT COUNT(*) FROM tasks WHERE status='queued'").fetchone()[0] == 0
    (run / "state.sqlite").rename(run / "saved.sqlite")
    archived = diagnostics(queue, task)
    assert not archived["trace"].get("turns") and archived["summary"]["detail"] == "stalled"
    assert not (run / "state.sqlite").exists()


def test_transcript_of_pending_and_retried_dispatch(tmp_path):
    from rarecell.dashboard_transcript import transcript
    from rarecell.state import State

    root, run = tmp_path / "queue", tmp_path / "run"
    jobs.initialize(root)
    run.mkdir()
    mid = "run1:T0001"
    (run / "checkpoint.json").write_text(json.dumps({"run_id": "run1", "policy": "live",
                                                    "annotations": str(tmp_path / "labels")}))
    tid = jobs.submit(root, tmp_path, {"kind": "inspect", "run": str(run)})
    task = jobs.get_task(root, tid)
    db = State(run, tmp_path / "cache")
    try:
        db.set_meta("run_id", "run1")
        context_ref = db.artifact({"candidate_id": "Hx", "turn_instruction": "inspect"})
        attempt = mid + "/api/1"
        db.put("ProviderDispatch", attempt, {"id": attempt, "message_id": mid, "context_ref": context_ref})
        tag, body = transcript(root, task)
        assert b"T0001" in body and b"turn_instruction" in body and "尚无已保存的回复".encode() in body
        db.put("ProviderFailure", attempt, {"id": attempt, "status": "temporary_failure", "detail": "HTTP 503"})
        second = mid + "/api/2"
        db.put("ProviderDispatch", second, {"id": second, "message_id": mid, "context_ref": context_ref})
        db.put("Message", mid, {"message_id": mid, "raw": "response text", "complete": True,
                                "usage": {"total_tokens": 53}})
        new_tag, body = transcript(root, task, tag)
        assert tag != new_tag and all(item in body for item in [b"/api/1", b"/api/2", b"HTTP 503", b"response text", b"total_tokens"])
        batch = tmp_path / "labels" / "batches" / ("f" * 64)
        batch.mkdir(parents=True)
        (batch / "request.json").write_text('{"groups":[{"candidate_id":"H-ann"}]}')
        (batch / "response.json").write_text('{"raw_response":"recorded annotation","prompt":"original prompt","rows":[]}')
        (batch / "response.json.failure.json").write_text('{"reason":"former failure"}')
        _, body = transcript(root, task, new_tag)
        assert all(item in body for item in [b"recorded annotation", b"original prompt", b"H-ann", b"former failure"])
    finally:
        db.close()
