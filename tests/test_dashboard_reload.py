"""Exercise a real backend swap with a leased scientific task in flight."""
import fcntl
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from urllib.request import Request, urlopen

from rarecell import job_queue as jobs


def until(read, predicate, timeout=30):
    end = time.monotonic() + timeout
    while True:
        value = read()
        if predicate(value):
            return value
        if time.monotonic() > end:
            raise AssertionError("timed out waiting for dashboard")
        time.sleep(.15)


def test_hot_reload_keeps_listener_login_and_running_task(tmp_path):
    source = Path(jobs.__file__).resolve().parent
    package = tmp_path / "rarecell"
    shutil.copytree(source, package, ignore=shutil.ignore_patterns("__pycache__"))
    root = tmp_path / "queue"
    jobs.initialize(root)
    cfg = tmp_path / "config.yaml"
    cfg.write_text("prepare: {views: [V_LINEAR], graph_backend: igraph}\n")
    dataset = source.parent / "examples/minimal_sparse/input/dataset.yaml"
    tid = jobs.submit(root, tmp_path, {"kind": "prepare", "dataset": str(dataset), "config": str(cfg)})
    errlog = tmp_path / "server.log"
    with (root / tid / ".snapshot.prepare.lock").open("a") as lock, errlog.open("w") as errors:
        fcntl.flock(lock, fcntl.LOCK_EX)
        proc = subprocess.Popen([sys.executable, "-m", "rarecell", "ui", "--workspace", str(root), "--port", "0"],
                                env=dict(os.environ, PYTHONPATH=str(tmp_path)), cwd=tmp_path,
                                stdout=subprocess.PIPE, stderr=errors, text=True, start_new_session=True)
        try:
            assert "RareCell" in proc.stdout.readline()
            url, token = proc.stdout.readline().strip().split("#")
            def info():
                with urlopen(Request(url + "api/tasks", headers={"Authorization": "Bearer " + token}), timeout=5) as reply:
                    return json.load(reply)
            first = until(info, lambda r: r["tasks"][0]["stage"] == "prepare")
            start = first["tasks"][0]["started"]
            original = (package / "dashboard.py").read_bytes()
            (package / "dashboard.py").write_bytes(original + b"\nsyntax error here\n")
            until(errlog.read_text, lambda log: "原版本继续服务" in log)
            assert info()["ui_revision"] == first["ui_revision"]
            (package / "dashboard.py").write_bytes(original + b"\n# valid update\n")
            css = package / "web/app.css"
            css.write_bytes(css.read_bytes() + b"\n/* observed update */\n")
            second = until(info, lambda r: r["ui_revision"] != first["ui_revision"])
            assert second["tasks"][0]["status"] == "running" and second["tasks"][0]["started"] == start
            with urlopen(url + "app.css", timeout=5) as reply:
                assert b"observed update" in reply.read()
            fcntl.flock(lock, fcntl.LOCK_UN)
            end = until(info, lambda r: r["tasks"][0]["status"] not in jobs.ACTIVE)
            assert end["tasks"][0]["status"] == "succeeded" and end["tasks"][0]["started"] == start
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)
            proc.send_signal(signal.SIGINT)
            proc.wait(timeout=15)
