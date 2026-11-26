"""Authenticated loopback UI and lifecycle of detached queue workers."""
import fcntl
import json
import mimetypes
from pathlib import Path
import re
import secrets
import shutil
import socket
import subprocess
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlsplit

import yaml

from . import job_queue as jobs
from .config import ROUTED_BUDGET_STOPS, read_json
from .dashboard_diagnostics import diagnostics
from .dashboard_transcript import transcript
from .dashboard_reload import serve as serve
from .runtime import project_model, save_project_model
from .single_agent_data import COPIED_OUT, RELEASE_OUT

WEB = Path(__file__).with_name("web")
ASSETS = {"/": "index.html", "/app.css": "app.css", "/app.js": "app.js", "/diagnostics.js": "diagnostics.js"}


def imported_views(project):
    choices = []
    for root, verified in ((COPIED_OUT, False), (RELEASE_OUT, True)):
        for path in sorted((project / root).glob("*/*/input/dataset.yaml")):
            folder = path.parent.parent
            if folder.name.startswith(".") or folder.parent.name.startswith(".") or not path.is_file():
                continue
            gt = folder / "evaluation_fixture/gt.tsv"
            label = f"{folder.parent.name} / {folder.name}"
            choices.append({"label": label, "dataset": str(path), "gt": str(gt) if verified and gt.is_file() else ""})
            if verified:
                for filename, title in (("labels_l1.tsv", "作者 L1"), ("labels_l3.tsv", "作者 L3"),
                                        ("targets.tsv", "逐目标 F1")):
                    reference = gt.with_name(filename)
                    if reference.is_file():
                        choices.append({"label": label, "dataset": str(path), "gt": str(reference), "profile": title})
    return choices


def artifacts(root, task):
    dest = root / task["id"]
    run = jobs.run_path(root, task)
    base = dest / "demo" if task["kind"] == "demo" else dest
    locations = {base / "evaluation": ["metrics.json", "targets.tsv"],
                 base / "snapshot": ["manifest.json"],
                 base / ("annotation" if task["kind"] == "demo" else "annotations"): ["annotation_status.json"]}
    if run:
        locations[run] = ["run_summary.json", "usage.json", "partition.tsv", "groups.tsv", "seed_reviews.tsv"]
    return {name: directory / name for directory, names in locations.items()
            for name in names if (directory / name).is_file()}


def tail(path, size=32000):
    try:
        with path.open("rb") as f:
            f.seek(max(0, f.seek(0, 2) - size))
            return f.read(size).decode("utf-8", errors="replace")
    except FileNotFoundError:
        return ""


def detail(root, task_id):
    task = jobs.get_task(root, task_id)
    outputs = artifacts(root, task)
    task["files"] = list(outputs)
    task["output"] = str(root / task_id)
    task["log"] = tail(root / task_id / "task.log")
    run = jobs.run_path(root, task)
    task["events"] = tail(run / "events.jsonl", 6000).splitlines()[-8:] if run else []
    settings = root / task_id / "settings.yaml"
    if not settings.is_file() and run:
        settings = run / "config_resolved.yaml"
    cfg = yaml.safe_load(settings.read_text()) if settings.is_file() else {}
    task["models"] = {section: cfg.get(section, {}).get("model") for section in ("provider", "annotation")}
    task["run_path"] = str(run) if run and (run / "checkpoint.json").is_file() else None
    summary = read_json(outputs["run_summary.json"]) if "run_summary.json" in outputs else {}
    task["resumable"] = (task["kind"] not in {"evaluate", "inspect"} and bool(task["run_path"]) and
                         task["status"] not in jobs.ACTIVE and (summary.get("resumable", True)
                         or (not summary.get("submitted") and summary.get("detail") in ROUTED_BUDGET_STOPS)))
    task["metrics"] = {key: summary[key] for key in ("n_cells", "n_candidates", "n_groups", "unresolved_cells")
                       if key in summary}
    if "usage.json" in outputs:
        usage = read_json(outputs["usage.json"])
        for key in ("provider_dispatches", "reported_total_tokens"):
            if usage.get(key) is not None:
                task["metrics"][key] = usage[key]
    return task


class Dashboard(ThreadingHTTPServer):
    daemon_threads = False  # Finish in-flight requests when a hot update retires this backend.

    def __init__(self, root, project, port=8787, *, listener=None, token=None):
        self.root = Path(root).resolve()
        self.project = Path(project).resolve()
        jobs.initialize(self.root)
        self.token = token or secrets.token_urlsafe(32)
        self.revision = secrets.token_hex(12)
        self.assets = {url: (WEB / name).read_bytes() for url, name in ASSETS.items()}
        self.assets["/"] = self.assets["/"].replace(b"__UI_REVISION__", self.revision.encode())
        self.worker = None
        super().__init__(listener.getsockname() if listener else ("127.0.0.1", port), Handler,
                         bind_and_activate=listener is None)
        if listener is not None:
            self.socket.close()
            self.socket = listener
            self.server_name, self.server_port = listener.getsockname()

    def service_actions(self):
        if self.worker is not None and self.worker.poll() is None:
            return
        with jobs.database(self.root) as db:
            pending = db.execute("SELECT 1 FROM tasks WHERE status IN ('queued','running','stopping') LIMIT 1").fetchone()
        if not pending:
            return
        # A worker from a previous dashboard can still be draining the queue.
        with (self.root / "worker.lock").open("a") as lease:
            try:
                fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return
        with (self.root / "worker.log").open("ab") as log:
            self.worker = subprocess.Popen([sys.executable, "-m", "rarecell.job_queue", "worker", str(self.root)],
                                           cwd=self.project, stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                                           start_new_session=True)


class Handler(BaseHTTPRequestHandler):
    timeout = 30

    # No access logs containing capability URLs, paths or user-supplied headers.
    def log_message(self, fmt, *args):
        pass

    def reply(self, data, status=200, mime="application/json; charset=utf-8"):
        if not isinstance(data, bytes):
            data = json.dumps(data, ensure_ascii=False).encode()
        self.send_response(status)
        self.headers_for(mime, len(data))
        self.end_headers()
        self.wfile.write(data)

    def headers_for(self, mime, size):
        self.send_header("Content-Type", mime)
        self.send_header("Content-Length", str(size))
        self.send_header("Cache-Control", "no-store")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy", "default-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'none'")

    def route(self, method):
        path = urlsplit(self.path).path
        if method == "GET" and path in ASSETS:
            asset = WEB / ASSETS[path]
            return self.reply(self.server.assets[path], mime=mimetypes.guess_type(asset.name)[0] + "; charset=utf-8")
        key = self.headers.get("Authorization", "")
        if not secrets.compare_digest(key.encode(), ("Bearer " + self.server.token).encode()):
            return self.reply({"error": "请使用启动终端显示的完整链接重新打开控制台"}, 401)
        root = self.server.root
        if method == "POST" and path == "/api/model":
            body = self.body()
            if not isinstance(body, dict) or set(body) != {"model"}:
                raise ValueError("请提交模型 ID")
            save_project_model(self.server.project, body["model"])
            return self.reply({"model": body["model"]})
        if method == "GET" and path == "/api/datasets":
            return self.reply(imported_views(self.server.project))
        if path == "/api/tasks":
            if method == "GET":
                return self.reply({"tasks": jobs.list_tasks(root), "host": socket.gethostname(),
                                   "workspace": str(root), "project": str(self.server.project), "ui_revision": self.server.revision,
                                   "model": project_model(self.server.project)})
            task_id = jobs.submit(root, self.server.project, self.body())
            return self.reply({"id": task_id}, 201)
        full = re.fullmatch(r"/api/tasks/([0-9a-f]{32})/diagnostics/transcript", path)
        if method == "GET" and full:
            tag, body = transcript(root, jobs.get_task(root, full[1]), self.headers.get("If-None-Match"))
            self.send_response(304 if body is None else 200)
            if tag:
                self.send_header("ETag", tag)
            if body is not None:
                self.headers_for("text/plain; charset=utf-8", len(body))
            self.end_headers()
            if body is not None:
                self.wfile.write(body)
            return
        diag = re.fullmatch(r"/api/tasks/([0-9a-f]{32})/diagnostics(?:/turn/([1-9][0-9]{0,8}))?", path)
        if method == "GET" and diag:
            task = jobs.get_task(root, diag[1])
            page = parse_qs(urlsplit(self.path).query).get("page", ["0"])[0]
            if not re.fullmatch(r"[0-9]{1,9}", page):
                raise ValueError("页码无效")
            return self.reply(diagnostics(root, task, int(page), int(diag[2]) if diag[2] else None))
        m = re.fullmatch(r"/api/tasks/([0-9a-f]{32})(?:/(stop|files/([a-z_.]+)))?", path)
        if not m:
            return self.reply({"error": "地址不存在"}, 404)
        task_id, action, filename = m.groups()
        if method == "POST" and action == "stop":
            jobs.cancel(root, task_id)
            return self.reply({"ok": True})
        if method == "GET" and not action:
            return self.reply(detail(root, task_id))
        if method == "GET" and filename:
            file = artifacts(root, jobs.get_task(root, task_id)).get(filename)
            if file:
                with file.open("rb") as fh:
                    self.send_response(200)
                    self.headers_for("application/octet-stream", file.stat().st_size)
                    self.send_header("Content-Disposition", 'attachment; filename="%s"' % filename)
                    self.end_headers()
                    shutil.copyfileobj(fh, self.wfile)
                return
        self.reply({"error": "地址不存在"}, 404)

    def body(self):
        size = int(self.headers.get("Content-Length", 0))
        if size <= 0 or size > 16384 or self.headers.get_content_type() != "application/json":
            raise ValueError("请提交 JSON 格式的任务参数（最多 16 KiB）")
        return json.loads(self.rfile.read(size))

    def handle_method(self, method):
        try:
            self.route(method)
        except (BrokenPipeError, ConnectionResetError):
            pass
        except KeyError as exc:
            self.reply({"error": str(exc)}, 404)
        except (ValueError, OSError, TypeError) as exc:
            self.reply({"error": str(exc)}, 400)
        except Exception:
            self.reply({"error": "读取任务失败，请检查任务目录及 YAML 配置"}, 500)

    def do_GET(self):
        self.handle_method("GET")

    def do_POST(self):
        self.handle_method("POST")


def serve_inherited(root, project, listener_fd, control_fd):
    with socket.socket(fileno=control_fd) as control, socket.socket(fileno=listener_fd) as listener:
        # The supervisor alone owns the login across reloads. EOF also stops this
        # process when the supervisor exits unexpectedly.
        with control.makefile("rb", buffering=0) as stream:
            token = stream.readline(128).decode("ascii").strip()
        if not token:
            return
        with Dashboard(root, project, listener=listener, token=token) as server:
            def retire():
                control.recv(1)
                server.shutdown()

            threading.Thread(target=retire, daemon=True).start()
            control.sendall(b"1")
            server.serve_forever(poll_interval=.2)


if __name__ == "__main__":
    serve_inherited(Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
