"""Persistent, single-worker task queue for the local dashboard (Linux/POSIX)."""
import fcntl
import json
import os
from pathlib import Path
import shlex
import signal
import sqlite3
import subprocess
import sys
import time
import traceback
import uuid
from contextlib import closing, contextmanager

import yaml

from .config import ACTIVE_MODES, ANNOTATED_MODES, MODES, ROUTED_BUDGET_STOPS, config, read_json, write_json
from .runtime import model_overrides, project_model

KINDS = {"pipeline", "prepare", "annotate", "discover", "resume", "evaluate", "demo", "inspect"}
ACTIVE = {"queued", "running", "stopping"}


@contextmanager
def database(root):
    with closing(sqlite3.connect(str(root / "tasks.sqlite"), timeout=15)) as db:
        db.row_factory = sqlite3.Row
        with db:
            yield db


def initialize(root):
    root.mkdir(parents=True, mode=0o700, exist_ok=True)
    with database(root) as db:
        db.execute("""CREATE TABLE IF NOT EXISTS tasks (
            id TEXT PRIMARY KEY, name TEXT NOT NULL, kind TEXT NOT NULL, mode TEXT,
            status TEXT NOT NULL, stage TEXT, created REAL NOT NULL, started REAL, ended REAL,
            detail TEXT NOT NULL DEFAULT '', spec TEXT NOT NULL)""")


def get_task(root, task_id):
    with database(root) as db:
        row = db.execute("SELECT * FROM tasks WHERE id=?", (task_id,)).fetchone()
    if row is None:
        raise KeyError("任务不存在")
    task = dict(row)
    task["spec"] = json.loads(task["spec"])
    return task


def list_tasks(root):
    with database(root) as db:
        return [dict(row) for row in db.execute(
            "SELECT id,name,kind,mode,status,stage,created,started,ended,detail FROM tasks ORDER BY created DESC, id")]


def input_path(value, project, directory=False):
    path = (project / Path(value).expanduser()).resolve()
    if not (path.is_dir() if directory else path.is_file()):
        raise ValueError("路径不存在或类型不符：" + str(path))
    return str(path)


def submit(root, project, form):
    if not isinstance(form, dict) or not all(isinstance(v, str) and len(v) <= 4096 for v in form.values()):
        raise ValueError("任务参数格式错误")
    form = {k: v.strip() for k, v in form.items() if v.strip()}
    kind, mode = form.get("kind", ""), form.get("mode", "A3")
    if kind not in KINDS or mode not in MODES:
        raise ValueError("请选择有效的任务类型和实验臂")
    name = form.get("name", "")
    if len(name) > 120:
        raise ValueError("任务名最多 120 个字符")
    spec = {}
    required = {"prepare": ["dataset"], "annotate": ["snapshot"], "discover": ["snapshot"],
                "resume": ["run"], "inspect": ["run"], "evaluate": ["run", "gt"]}.get(kind, [])
    if kind == "pipeline" and bool(form.get("dataset")) == bool(form.get("snapshot")):
        raise ValueError("完整流程需要选择一个 dataset.yaml 或已有快照")
    if any(not form.get(k) for k in required):
        raise ValueError("缺少必填输入：" + "、".join(required))
    path_fields = {
        "pipeline": ("dataset", "snapshot", "annotations", "reuse_annotations", "config", "gt"),
        "prepare": ("dataset", "config"), "annotate": ("snapshot", "config", "reuse_annotations"),
        "discover": ("snapshot", "annotations", "config", "gt"), "evaluate": ("run", "gt"),
        "resume": ("run", "gt"), "inspect": ("run",), "demo": (),
    }
    if kind == "resume" and form.get("retry_unknown") == "yes":
        spec["retry_unknown"] = True
    if kind == "resume":
        for key in ROUTED_BUDGET_STOPS.values():
            if key in form:
                if not form[key].isascii() or not form[key].isdecimal():
                    raise ValueError("续跑上限须为非负整数")
                spec[key] = int(form[key])
    for key in path_fields[kind]:
        if form.get(key):
            spec[key] = input_path(form[key], project, key in {"snapshot", "annotations", "reuse_annotations", "run"})
    if kind in {"resume", "inspect"} and not (Path(spec["run"]) / "checkpoint.json").is_file():
        raise ValueError("运行目录中缺少 checkpoint.json")
    if kind == "evaluate" and not (Path(spec["run"]) / "frozen.json").is_file():
        raise ValueError("运行目录中缺少 frozen.json")
    if kind in {"pipeline", "prepare", "discover", "annotate"}:
        override = model_overrides(project_model(project))
        for field, section in (("model", "provider"), ("annotation_model", "annotation")):
            if form.get(field):
                override[section] = {"model": form[field]}
        if kind in {"discover", "pipeline"}:
            override["mode"] = mode
            setting = form.get("exhaustive_seed_review")
            if setting in {"yes", "no"}:
                override["exhaustive_seed_review"] = setting == "yes"
            elif setting:
                raise ValueError("无效的种子调查选项")
        if kind in {"annotate", "pipeline"} and form.get("panel_mode"):
            override.setdefault("annotation", {})["panel_mode"] = form["panel_mode"]
        if kind in {"annotate", "pipeline"} and form.get("include_singletons"):
            if form["include_singletons"] not in {"yes", "no"}:
                raise ValueError("无效的单细胞注释选项")
            override.setdefault("annotation", {})["include_singletons"] = form["include_singletons"] == "yes"
        cfg = config(spec.get("config"), override)
        if kind in {"pipeline", "discover"} and cfg.get("exhaustive_seed_review") and mode not in ACTIVE_MODES:
            raise ValueError("强制逐种子调查需要 A2、A3 或 A3_adaptive")
        if kind in {"pipeline", "discover"} and not cfg["provider"].get("model"):
            raise ValueError("请填写发现模型，或选择包含 provider.model 的配置文件")
        needs_annotation = kind == "annotate" or (kind == "pipeline" and mode in ANNOTATED_MODES and not spec.get("annotations"))
        if spec.get("reuse_annotations") and (not needs_annotation or not spec.get("snapshot")):
            raise ValueError("复用注释批次需要在已有快照上运行注释")
        if needs_annotation and (not cfg["annotation"].get("model") or not cfg["annotation"].get("enabled")):
            raise ValueError("请启用注释并指定 annotation.model，或选择已有注释目录")
        if kind == "discover" and mode in ANNOTATED_MODES and not spec.get("annotations"):
            raise ValueError("该实验臂需要已有注释目录")
    else:
        cfg = None
    task_id = uuid.uuid4().hex
    dest = root / task_id
    if not name:
        source = spec.get("dataset") or spec.get("snapshot") or spec.get("run")
        name = (Path(source).parent.name if spec.get("dataset") else Path(source).name) if source else "Core demo"
    stored_mode = mode if kind in {"pipeline", "discover"} else None
    with database(root) as db:
        db.execute("BEGIN IMMEDIATE")
        if kind == "inspect":
            linked = db.execute("SELECT id FROM tasks WHERE kind='inspect' AND spec=?", (json.dumps(spec),)).fetchone()
            if linked:
                return linked["id"]
        if kind == "resume":
            active = db.execute("SELECT id,kind,spec FROM tasks WHERE status IN ('queued','running','stopping')")
            for row in active:
                other = dict(row)
                other["spec"] = json.loads(other["spec"])
                if other["kind"] != "evaluate" and str(run_path(root, other)) == spec["run"]:
                    raise ValueError("该运行已有排队或执行中的任务")
        dest.mkdir(mode=0o700)
        if cfg is not None:
            (dest / "settings.yaml").write_text(yaml.safe_dump(cfg, allow_unicode=True), encoding="utf-8")
        db.execute("INSERT INTO tasks(id,name,kind,mode,status,created,spec) VALUES(?,?,?,?,?,?,?)",
                   (task_id, name, kind, stored_mode, "linked" if kind == "inspect" else "queued",
                    time.time(), json.dumps(spec)))
    return task_id


def run_path(root, task):
    if task["kind"] in {"resume", "evaluate", "inspect"}:
        return Path(task["spec"]["run"])
    if task["kind"] in {"pipeline", "discover"}:
        return root / task["id"] / "run"
    if task["kind"] == "demo":
        return root / task["id"] / "demo" / "run"
    return None


def steps(root, task):
    spec, kind, mode = task["spec"], task["kind"], task["mode"]
    dest = root / task["id"]
    snapshot = str(spec.get("snapshot", dest / "snapshot"))
    annotations = str(spec.get("annotations", dest / "annotations"))
    run = str(run_path(root, task))
    cfg = ["--config", str(dest / "settings.yaml")]
    if kind == "demo":
        yield ["demo", "--mode", "core", "--out", str(dest / "demo")]
        return
    if kind == "prepare" or (kind == "pipeline" and "snapshot" not in spec):
        yield ["prepare", "--dataset", spec["dataset"], "--out", snapshot] + cfg
    if kind == "annotate" or (kind == "pipeline" and mode in ANNOTATED_MODES and "annotations" not in spec):
        reuse = ["--reuse-annotations", spec["reuse_annotations"]] if spec.get("reuse_annotations") else []
        yield ["annotate", "--snapshot", snapshot, "--out", annotations] + reuse + cfg
    if kind in {"discover", "pipeline"}:
        args = ["discover", "--snapshot", snapshot, "--mode", mode, "--policy", "live", "--run", run]
        if mode in ANNOTATED_MODES:
            args += ["--annotations", annotations]
        yield args + cfg
    if kind == "resume":
        args = ["resume", "--run", run] + (["--retry-unknown"] if spec.get("retry_unknown") else [])
        for key in ROUTED_BUDGET_STOPS.values():
            if key in spec:
                args.extend(["--" + key.replace("_", "-"), str(spec[key])])
        yield args
    if kind == "evaluate" or (kind in {"pipeline", "resume", "discover"} and spec.get("gt")):
        yield ["evaluate", "--run", run, "--gt", spec["gt"], "--out", str(dest / "evaluation")]


def cancel(root, task_id):
    with database(root) as db:
        cur = db.execute("""UPDATE tasks SET status=CASE status WHEN 'queued' THEN 'stopped' ELSE 'stopping' END,
            ended=CASE status WHEN 'queued' THEN ? ELSE ended END WHERE id=? AND status IN ('queued','running')""",
                         (time.time(), task_id))
        if not cur.rowcount:
            raise ValueError("任务已结束或正在停止")


def finish(root, task_id, status, detail=""):
    with database(root) as db:
        db.execute("""UPDATE tasks SET status=CASE status WHEN 'stopping' THEN 'stopped' ELSE ? END,
                      detail=?,ended=? WHERE id=?""", (status, detail, time.time(), task_id))


def execute_step(result_path, argv):
    # Dispatch through the same parser and scientific implementation as the CLI.
    from .cli import dispatch, parser
    result = dispatch(parser().parse_args(argv))
    report = result.get("summary", result)
    status = report.get("status", "complete")
    detail = report.get("detail") or report.get("issue") or ""
    write_json(result_path, {"status": status, "detail": str(detail)})
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


def stop_process(proc):
    # Stop R descendants in the same, dedicated process group as well.
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        until = time.monotonic() + 3
        while time.monotonic() < until:
            proc.poll()
            os.killpg(proc.pid, 0)
            time.sleep(0.1)
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.wait()


def execute_task(root, task, lease):
    dest = root / task["id"]
    outcome, reason = "succeeded", ""
    with (dest / "task.log").open("a", encoding="utf-8", buffering=1) as log:
        for argv in steps(root, task):
            with database(root) as db:
                changed = db.execute("UPDATE tasks SET stage=? WHERE id=? AND status='running'",
                                     (argv[0], task["id"])).rowcount
            if not changed:
                return "stopped", ""
            print("$ rarecell " + shlex.join(argv), file=log)
            result_path = dest / (argv[0] + "-result.json")
            command = [sys.executable, "-u", "-m", "rarecell.job_queue", "step", str(result_path)] + argv
            with subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                  start_new_session=True, pass_fds=(lease,)) as proc:
                try:
                    while proc.poll() is None:
                        if get_task(root, task["id"])["status"] == "stopping":
                            stop_process(proc)
                            return "stopped", ""
                        time.sleep(0.35)
                except BaseException:
                    stop_process(proc)
                    raise
            if proc.returncode or not result_path.exists():
                return "failed", "阶段退出码 %s，请查看日志" % proc.returncode
            report = read_json(result_path)
            if report["detail"] == "paused":
                return "paused", "paused"
            if report["status"] in {"blocked", "not_run"}:
                return "blocked", report["detail"] or report["status"]
            if report["status"] == "partial":
                outcome, reason = "partial", report["detail"] or argv[0]
    return outcome, reason


def worker(root):
    with (root / "worker.lock").open("a") as lease:
        try:
            fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        # The CLI child inherits the lease. Reaching here implies the old
        # executor and its scientific Python process have both exited.
        with database(root) as db:
            db.execute("""UPDATE tasks SET status='interrupted',ended=?,detail='执行进程已退出'
                          WHERE status IN ('running','stopping')""", (time.time(),))
        while True:
            with database(root) as db:
                db.execute("BEGIN IMMEDIATE")
                row = db.execute("SELECT id FROM tasks WHERE status='queued' ORDER BY created,id LIMIT 1").fetchone()
                if not row:
                    return
                task_id = row["id"]
                db.execute("UPDATE tasks SET status='running',started=? WHERE id=?", (time.time(), task_id))
            try:
                status, detail = execute_task(root, get_task(root, task_id), lease.fileno())
            except Exception as exc:
                with (root / task_id / "task.log").open("a") as log:
                    traceback.print_exc(file=log)
                status, detail = "failed", str(exc)
            finish(root, task_id, status, detail)


if __name__ == "__main__":
    if sys.argv[1] == "worker":
        worker(Path(sys.argv[2]))
    elif sys.argv[1] == "step":
        execute_step(Path(sys.argv[2]), sys.argv[3:])
