"""Render only the persisted, observable model exchanges of a task."""
import hashlib
import json
import re
import sqlite3
from contextlib import closing

from .dashboard_diagnostics import fields, optional_json, task_locations, turn_count, turn_detail

BATCH_FILES = {
    "request.json": "批次请求",
    "dispatch.json": "派发记录",
    "response.json": "上游回复",
    "response.json.failure.json": "失败记录",
}


def file_stamp(path):
    try:
        stat = path.stat()
    except FileNotFoundError:
        return None
    return (stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def block(label, value):
    if isinstance(value, str):
        body = value
    else:
        body = json.dumps(value, ensure_ascii=False, indent=2)
    return label + "\n" + body + "\n"


def annotation_blocks(directory, batches):
    state = optional_json(directory / "annotation_status.json")
    summary = fields(state, "status detail settings annotated annotated_singletons panel_slots partition_issue skipped errors cached_batches "
                     "newly_completed_batches provider_calls_observed_at_wire annotation_cost_usd")
    yield block("OntoAnno 注释", summary or "尚无阶段报告")
    if batches:
        yield "批次按请求保存时间排列，缓存复用数量见阶段报告。"
    for folder in batches:
        yield "── 批次 " + folder.name + " ──"
        for name, title in BATCH_FILES.items():
            path = folder / name
            if path.is_symlink() or path.resolve().parent != folder:
                yield title + "文件位置异常。"
            elif path.is_file():
                try:
                    data = json.loads(path.read_bytes())
                except (FileNotFoundError, ValueError):
                    yield title + "暂不可读。"
                    continue
                if name == "response.json":
                    for key, label in (("prompt", "上游提示"), ("raw_response", "OntoAnno 原始回复")):
                        if key in data:
                            yield block(label, data.pop(key))
                    yield block("注释解析结果", data)
                else:
                    yield block(title, data)


def context_block(run, ref):
    match = re.fullmatch(r"artifacts/([0-9a-f]{64})\.json", ref or "")
    if match:
        path = run / ref
        if path.resolve().parent == run.resolve() / "artifacts" and not path.is_symlink():
            try:
                body = path.read_bytes()
            except FileNotFoundError:
                pass
            else:
                if hashlib.sha256(body).hexdigest() == match[1]:
                    return block("发送给发现模型的任务上下文", json.loads(body))
    return "该次派发的上下文文件缺失或校验失败。"


def discovery_blocks(run, checkpoint):
    yield block("发现运行", fields(checkpoint, "run_id mode policy prompt_version") or "暂无检查点")
    if not run or not (run / "state.sqlite").is_file():
        yield "暂无发现状态库。"
        return
    # Copy the recorded turns in one SQLite snapshot. Release its read lock
    # before reading potentially large context artifacts or serving the client.
    with closing(sqlite3.connect((run / "state.sqlite").resolve().as_uri() + "?mode=ro", uri=True, timeout=1)) as db:
        db.execute("BEGIN")
        turns = [turn_detail(db, i) for i in range(1, turn_count(db) + 1)]
        decisions = {key: json.loads(body) for key, body in db.execute(
            "SELECT id,body FROM record WHERE kind='Decision'")}
        progress, events = {}, []
        for kind, payload in db.execute("SELECT type,payload FROM event WHERE type != 'action' ORDER BY seq"):
            event = json.loads(payload)
            if kind == "turn":
                progress[event["turn"]] = event["progress"]
            else:
                events.append(dict(type=kind, **event))
    yield "动作索引" if turns else "尚无已记录的发现轮次。"
    for number, turn in enumerate(turns, 1):
        for receipt in turn["receipts"]:
            index = receipt["index"]
            command = turn["commands"][index] if index < len(turn["commands"]) else {}
            result = receipt["result"]
            # The controller can replace the first over-budget command with a rejection.
            name = result.get("action") or ("未执行" if result.get("error") else command.get("action"))
            target = command.get("args", {}).get("candidate_id", "")
            yield " ".join(str(x) for x in (f"T{number:04d}", f"[{index}]", name or "未解析", target,
                                            result.get("status", "ok"), result.get("error", "")) if x)
    for number, turn in enumerate(turns, 1):
        yield "── T%04d ──" % number
        yield {True: "本轮有进展", False: "本轮无新进展"}.get(progress.get(number), "尚无轮次结束记录")
        failures = {item["id"]: item for item in turn["failures"]}
        for attempt in turn["dispatches"]:
            yield block("模型请求", attempt)
            yield context_block(run, attempt.get("context_ref"))
            if attempt["id"] in failures:
                yield block("调用失败", failures[attempt["id"]])
        message = turn["message"]
        if message:
            label = "模型原始回复" if checkpoint.get("policy") == "live" else "策略原始输出"
            yield block(label, message.get("raw") or "（无已保存的回复正文）")
            yield block("完成状态与用量", fields(message, "complete usage"))
        else:
            yield "尚无已保存的回复。"
        if turn["seed_followup"]:
            yield block("控制器种子补查", turn["seed_followup"])
        for i, command in enumerate(turn["commands"]):
            yield block("解析请求 [%d]" % i, command)
        for receipt in turn["receipts"]:
            yield block("动作回执 [%d]" % receipt["index"], receipt)
            decision = decisions.get(receipt["result"].get("decision_id"))
            if decision:
                yield block("已保存决策", decision)
    for event in events:
        yield block("运行事件", event)
        if event["type"] == "provisional_submission" and event.get("decision_id") in decisions:
            yield block("停滞时保存的暂定提交", decisions[event["decision_id"]])
    summary = optional_json(run / "run_summary.json")
    if summary:
        yield block("发现结果", fields(summary, "status detail issue submitted submission_status submission_source n_groups n_screened_seeds n_expanded_seeds n_controller_seed_followups seed_review_coverage closeout_trigger n_pending_seeds"))


def transcript(root, task, previous=None):
    run, _, checkpoint, ann_dir = task_locations(root, task)
    batches_dir = ann_dir.resolve() / "batches"
    batches = []
    if not batches_dir.is_symlink():
        for path in batches_dir.glob("*"):
            if re.fullmatch(r"[0-9a-f]{64}", path.name) and path.is_dir() and not path.is_symlink():
                batches.append(path)
    def request_order(folder):
        stamp = file_stamp(folder / "request.json")
        return (stamp[2] if stamp else 0, folder.name)
    batches.sort(key=request_order)
    paths = [ann_dir / "annotation_status.json"]
    paths.extend(p / name for p in batches for name in BATCH_FILES)
    if run:
        paths.extend(run / name for name in ("state.sqlite", "state.sqlite-wal", "checkpoint.json", "run_summary.json"))
    stamps = [(str(p), file_stamp(p)) for p in paths]
    tag = '"' + hashlib.sha256(repr((task["id"], stamps)).encode()).hexdigest() + '"'
    if tag == previous:
        return tag, None
    blocks = ["任务 " + task["id"], *annotation_blocks(ann_dir, batches)]
    try:
        blocks.extend(discovery_blocks(run, checkpoint))
    except sqlite3.DatabaseError as exc:
        # Preserve annotation output even while a new or copied database is unavailable.
        blocks.append("发现状态库暂不可读：" + str(exc))
        tag = None  # The next poll should retry even without a change in file metadata.
    return tag, "\n".join(blocks).encode("utf-8")
