"""Read persisted discovery decisions without opening a controller or running a stage."""
import csv
import json
import sqlite3
from contextlib import closing
from itertools import islice
from pathlib import Path

import yaml

from . import job_queue as jobs
from .config import fingerprint, hash_file, read_json
from .protocol import VERSION

PAGE_SIZE = 12


def optional_json(path):
    try:
        return read_json(path)
    except FileNotFoundError:
        return {}


def fields(record, names):
    return {k: record[k] for k in names.split() if k in record}


def evaluation(path, freeze=None):
    metrics = optional_json(path)
    if not metrics or (freeze is not None and (not freeze.is_file() or metrics.get("freeze_sha256") != hash_file(freeze))):
        return None
    final = metrics["submitted_partition"]
    targets = final["per_target"]
    groups = {r["group_id"] for r in targets if r.get("group_id")} | set(final["extra_groups"])
    claims = sum(g.startswith("claim/") for g in groups)
    residual = sum(g.startswith(("residual/", "unresolved/")) for g in groups)
    source = "unknown"
    if groups and len(groups) == final["n_groups"] == claims + residual:
        source = "claims" if not residual else "fallback" if not claims else "mixed"
    return dict(fields(final, "ari nmi macro_f1"), source=source,
                frozen=fields(metrics["run_summary"], "submitted fallback_view_resolution detail"),
                targets=[fields(r, "target n_target tp fp fn f1") for r in targets])


def record(db, kind, key):
    row = db.execute("SELECT body FROM record WHERE kind=? AND id=?", (kind, key)).fetchone()
    return json.loads(row[0]) if row else {}


def meta(db, key, default=None):
    row = db.execute("SELECT body FROM meta WHERE key=?", (key,)).fetchone()
    return json.loads(row[0]) if row else default


def turn_count(db):
    finished = meta(db, "finished_turn", 0)
    run_id = meta(db, "run_id")
    pending_id = f"{run_id}:T{finished+1:04d}"
    pending = db.execute("SELECT 1 FROM record WHERE (kind='Message' AND id=?) OR "
                         "(kind='ProviderDispatch' AND id=?)", (pending_id, pending_id + "/api/1")).fetchone()
    return finished + int(bool(pending))


def turn_page(db, page):
    total = turn_count(db)
    pages = max(1, (total + PAGE_SIZE - 1) // PAGE_SIZE)
    page = min(page or pages, pages)
    first = (page - 1) * PAGE_SIZE + 1
    last = min(total, page * PAGE_SIZE)
    turns = {i: {"number": i, "progress": None, "actions": []} for i in range(first, last + 1)}
    start_seq = 0
    if first > 1:
        boundary = db.execute("SELECT seq FROM event WHERE type='turn' ORDER BY seq LIMIT 1 OFFSET ?",
                              (first - 2,)).fetchone()
        if boundary:
            start_seq = boundary[0]
    events = db.execute("SELECT type,payload FROM event WHERE seq>? AND type IN ('action','turn') ORDER BY seq",
                        (start_seq,))
    for kind, payload in events:
        event = json.loads(payload)
        number = int(event["message_id"].rsplit(":T", 1)[1])
        if number > last:
            break
        if number not in turns:
            continue
        if kind == "turn":
            turns[number]["progress"] = event["progress"]
        else:
            action = fields(event, "index action status changed")
            if event["status"] in {"blocked", "needs_clarification"}:
                receipt = record(db, "Receipt", event["receipt"])
                action["error"] = receipt.get("result", {}).get("error")
            turns[number]["actions"].append(action)
    return dict(page=page, pages=pages, total=total, stall=meta(db, "stall", 0), turns=list(turns.values()))


def turn_detail(db, number):
    mid = f"{meta(db, 'run_id')}:T{number:04d}"
    msg = record(db, "Message", mid)
    dispatches = [record(db, "ProviderDispatch", mid + f"/api/{i}") for i in (1, 2)]
    dispatches = [item for item in dispatches if item]
    failures = [record(db, "ProviderFailure", item["id"]) for item in dispatches]
    if not msg and not dispatches:
        raise KeyError("该轮次没有记录")
    commands = []
    if msg:
        current = VERSION + ("/seed" if record(db, "SeedBatch", mid) else "")
        for version in (current, "soft-protocol-1.4", "soft-protocol-1.2", "soft-protocol-1.1", "protocol/1"):
            key = fingerprint([version, msg["raw"], msg.get("complete", True)])
            saved = record(db, "ParsedMessage", key)
            if saved:
                commands = saved["commands"]
                break
    receipts = {}
    for key, body in db.execute("SELECT id,body FROM record WHERE kind='Receipt' AND id>=? AND id<?",
                                (mid + "/", mid + "0")):
        receipts[int(key.rsplit("/", 1)[1])] = fields(json.loads(body), "result progress metadata")
    return {"seed_followup": record(db, "SeedFollowup", mid), "message": fields(msg, "message_id raw complete usage"), "commands": commands,
            "receipts": [dict(index=i, **receipts[i]) for i in sorted(receipts)],
            "dispatches": dispatches, "failures": [fields(a, "id status detail") for a in failures if a]}


def task_locations(root, task):
    run = jobs.run_path(root, task)
    dest = root / task["id"]
    base = dest / "demo" if task["kind"] == "demo" else dest
    checkpoint = optional_json(run / "checkpoint.json") if run else {}
    spec = task["spec"]
    ann_dir = Path(spec.get("annotations") or checkpoint.get("annotations") or
                   base / ("annotation" if task["kind"] == "demo" else "annotations"))
    return run, base, checkpoint, ann_dir


def diagnostics(root, task, page=None, turn=None):
    run, base, checkpoint, ann_dir = task_locations(root, task)
    db_path = run / "state.sqlite" if run else None
    trace = {}
    if db_path and db_path.is_file():
        try:
            with closing(sqlite3.connect(db_path.resolve().as_uri() + "?mode=ro", uri=True, timeout=.4)) as db:
                db.execute("BEGIN")
                if turn is not None:
                    return turn_detail(db, turn)
                trace = turn_page(db, page)
                if db.execute("SELECT 1 FROM record WHERE kind='SeedBatch' LIMIT 1").fetchone():
                    counts = dict(db.execute("SELECT kind, count(*) FROM record WHERE kind IN "
                                             "('Seed','SeedMeasurement','SeedAssessment') GROUP BY kind"))
                    trace["seed_reviews"] = {"total": counts.get("Seed", 0), "written": counts.get("SeedAssessment", 0)}
        except sqlite3.DatabaseError as exc:
            trace = {"notice": "暂时无法读取逐轮记录，可刷新重试（%s）。" % exc}
            if turn is not None:
                raise ValueError(trace["notice"]) from exc
    elif turn is not None:
        raise KeyError("逐轮状态库尚不可用")
    snapshot = Path(task["spec"].get("snapshot") or checkpoint.get("snapshot") or base / "snapshot")
    summary = optional_json(run / "run_summary.json") if run else {}
    config_path = run / "config_resolved.yaml" if run else root / task["id"] / "settings.yaml"
    cfg = (yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}) if config_path.is_file() else {}
    manifest = optional_json(run / "manifest.json") if run else {}
    manifest = manifest or optional_json(snapshot / "manifest.json")
    ann = optional_json(ann_dir / "annotation_status.json")
    usage = optional_json(run / "usage.json") if run else {}
    seeds = []
    if summary.get("n_pending_seeds") and (run / "seeds.tsv").is_file():
        with (run / "seeds.tsv").open(encoding="utf-8", newline="") as stream:
            rows = csv.DictReader(stream, delimiter="\t")
            seeds = [fields(r, "candidate_id cell_id") for r in islice((r for r in rows if r["task_status"] == "pending"), 10)]
    inspected = task["kind"] == "inspect"
    scores = evaluation((run.parent if inspected else base) / "evaluation" / "metrics.json",
                        run / "frozen.json" if inspected else None)
    return {"trace": trace, "summary": fields(summary, "status detail issue submitted submission_status submission_source n_cells n_candidates n_groups "
                                           "unresolved_cells n_screened_seeds n_expanded_seeds n_controller_seed_followups seed_review_coverage closeout_trigger n_pending_seeds fallback_view_resolution annotation_status"),
            "pending_seeds": seeds, "annotation_model": (ann.get("settings") or {}).get("model"),
            "stall_limit": (cfg.get("budget") or {}).get("stall_turns"),
            "usage": fields(usage, "turns provider_dispatches protocol_errors action_errors new_queries compute_executions"),
            "prepare": fields(manifest, "status stages"), "annotation": fields(ann, "status detail annotated annotated_singletons partition_issue "
                "newly_completed_batches cached_batches errors skipped annotation_cost_usd provider_calls_observed_at_wire"),
            "evaluation": scores}
