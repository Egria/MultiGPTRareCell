"""Resume a routed discovery in the foreground and evaluate its submitted partition."""
import argparse
from contextlib import closing
import csv
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile


def readonly(path):
    return closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True))


def recent_error(run, previously_shown=None):
    with readonly(run / "state.sqlite") as db:
        receipts = db.execute("SELECT body FROM record WHERE kind='Receipt' ORDER BY id DESC LIMIT 200")
        for (body,) in receipts:
            receipt = json.loads(body)
            error = receipt["result"].get("error")
            if not error:
                continue
            mid = receipt["message_id"]
            if mid == previously_shown:
                return mid
            dispatch = db.execute("SELECT body FROM record WHERE kind='ProviderDispatch' AND id IN (?, ?) "
                                  "ORDER BY id DESC LIMIT 1", (mid + "/api/1", mid + "/api/2")).fetchone()
            phase = receipt["result"].get("action")
            if dispatch and not phase:
                path = run / json.loads(dispatch[0])["context_ref"]
                if path.is_file():
                    phase = json.loads(path.read_text()).get("phase")
            print(f"检查点记录的校验错误 {mid}  {phase or 'unknown'}  {error}", flush=True)
            return mid


def evaluation_gt(task):
    registry = task.parent / "tasks.sqlite"
    if registry.is_file():
        with readonly(registry) as db:
            row = db.execute("SELECT spec FROM tasks WHERE id=?", (task.name,)).fetchone()
        if row:
            gt = json.loads(row[0]).get("gt")
            if gt:
                return Path(gt)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--workers", type=int, help="Parallel independent requests; default uses checkpoint setting")
    parser.add_argument("--gt", type=Path, help="GT for evaluation; default uses the dashboard task's GT")
    args = parser.parse_args()
    run = args.run.resolve()
    spec = json.loads((run / "checkpoint.json").read_text())
    if spec["mode"] != "A3_routed":
        parser.error("需要 A3_routed 的运行目录")
    project = Path(__file__).resolve().parents[1]
    print(f"原地续跑 {spec['run_id']}  {run}\n进度文件 {run / 'events.jsonl'}", flush=True)
    last_error = recent_error(run)
    cmd = [sys.executable, "-u", "-m", "rarecell"]
    options = ["--routed-workers", str(args.workers)] if args.workers is not None else []
    result = subprocess.run([*cmd, "resume", "--traceback", "--run", str(run), *options], cwd=project)
    summary = run / "run_summary.json"
    if result.returncode not in (0, 2) or not summary.is_file():
        return result.returncode or 1
    saved = json.loads(summary.read_text())
    if saved["run_id"] != spec["run_id"] or not saved["submitted"]:
        recent_error(run, last_error)
        return result.returncode or 2
    gt = args.gt or evaluation_gt(run.parent)
    if not gt or not gt.is_file():
        print("已提交。可用 --gt 指定该数据集的参考文件后再次运行，直接评估已提交的结果。", flush=True)
        return result.returncode
    gt = gt.resolve()
    out = Path(tempfile.mkdtemp(prefix="submitted-eval-", dir=run.parent)) / "evaluation"
    print(f"评估参考 {gt}\n评估目录 {out}", flush=True)
    subprocess.run([*cmd, "evaluate", "--run", str(run), "--gt", str(gt), "--out", str(out)],
                   cwd=project, check=True)
    with (out / "targets.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if "mait" in row["target"].casefold():
                print("  ".join(f"{field}={row[field]}" for field in ("target", "group_id", "tp", "fp", "fn", "f1")))
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
