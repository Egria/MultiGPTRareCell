"""Offline metrics for frozen membership partitions. No discovery module imports this service."""
import argparse
import errno
import csv
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import yaml
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


from .config import fingerprint, hash_file as _sha


def _json(path):
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def _read(path, required):
    with path.open(encoding="utf-8", newline="") as f:
        rows = csv.DictReader(f, delimiter="\t")
        if not required.issubset(rows.fieldnames or []):
            raise ValueError("Missing columns in {}: {}".format(path.name, sorted(required)))
        for row in rows:
            if any(row.get(col) in (None, "") for col in required):
                raise ValueError("Missing scientific identifier or status in {}".format(path.name))
            yield row


def _inside(path, root):
    return path == root or root in path.parents


def _write(path, writer):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", prefix=".write-", delete=False,
                                     dir=path.parent, encoding="utf-8", newline="") as f:
        tmp = Path(f.name)
        try:
            writer(f)
            f.flush()
            os.fsync(f.fileno())
        except BaseException:
            tmp.unlink()
            raise
    os.replace(tmp, path)


def _json_write(path, data):
    _write(path, lambda f: json.dump(data, f, indent=2, ensure_ascii=False,
                                      sort_keys=True, allow_nan=False))


def _table(path, rows, fields):
    def write(f):
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    _write(path, write)


def _indexed(rows, key):
    result = {}
    for row in rows:
        name = row[key]
        if name in result:
            raise ValueError("Duplicate {}: {}".format(key, name))
        result[name] = row
    return result


def read_truth(path):
    result = {cell: row["label"] for cell, row in
              _indexed(_read(path, {"cell_id", "label"}), "cell_id").items()}
    if not result:
        raise ValueError("Empty ground truth")
    return result


def _verified(run):
    manifest_path = run / "frozen.json"
    if not _inside(manifest_path.resolve(), run):
        raise ValueError("Freeze manifest resolves outside discovery run")
    frozen = _json(manifest_path)
    files = frozen["files"]
    required = {"partition.tsv", "candidates.jsonl", "candidate_members.tsv",
                "usage.json", "run_summary.json", "artifacts/baselines.tsv"}
    if not required.issubset(files) or not frozen.get("snapshot_hash"):
        raise ValueError("Incomplete freeze manifest; missing scoring files or prepare hash")
    for name, digest in files.items():
        rel = Path(name)
        if rel.is_absolute() or ".." in rel.parts:
            raise ValueError("Unsafe relative name in freeze manifest")
        item = run / rel
        if not _inside(item.resolve(), run):
            raise ValueError("Frozen artifact resolves outside discovery run")
        if not item.is_file() or _sha(item) != digest:
            raise ValueError("Frozen artifact digest mismatch: {}".format(name))
    return frozen, _sha(manifest_path)


def _partition_scores(gt, predictions):
    if gt.keys() != predictions.keys():
        raise ValueError("Ground truth and partition cell_id sets differ")
    cells = sorted(gt)
    targets = sorted(set(gt.values()))
    by_group = sorted(set(predictions.values()))
    true_n, pred_n = Counter(gt.values()), Counter(predictions.values())
    overlap = Counter((gt[cell], predictions[cell]) for cell in cells)
    ti = {target: i for i, target in enumerate(targets)}
    gi = {group: i for i, group in enumerate(by_group)}
    n, m = len(ti), len(gi)
    row, col, weights = [], [], []
    fragment_count, mix = Counter(), Counter()
    for (target, group), tp in sorted(overlap.items()):
        row.append(ti[target])
        col.append(gi[group])
        weights.append(1 + 2 * tp / (true_n[target] + pred_n[group]))
        fragment_count[target] += 1
        mix[group] += 1
    # Each target can instead choose a separate dummy (zero F1).
    row.extend(range(n))
    col.extend(range(m, m + n))
    weights.extend([1.0] * n)
    matrix = coo_matrix((weights, (row, col)), shape=(n, m + n)).tocsr()
    rows, cols = min_weight_full_bipartite_matching(matrix, maximize=True)
    assignment = dict(zip(rows, cols))
    per_target, matches = [], {}
    for i, target in enumerate(targets):
        j = assignment[i]
        group = by_group[j] if j < m else None
        matches[target] = group
        tp = overlap[target, group] if group is not None else 0
        size = pred_n[group] if group is not None else 0
        per_target.append({
            "target": target, "group_id": group, "n_target": true_n[target],
            "n_pred": size, "tp": tp, "fp": size - tp, "fn": true_n[target] - tp,
            "precision": tp / size if size else 0.0, "recall": tp / true_n[target],
            "f1": 2 * tp / (size + true_n[target]), "fragments": fragment_count[target],
            "matched_group_gt_labels": mix[group] if group is not None else 0,
            "unmatched": group is None,
        })
    yt, yp = [gt[c] for c in cells], [predictions[c] for c in cells]
    used = {g for g in matches.values() if g is not None}
    return {
        "ari": adjusted_rand_score(yt, yp),
        "nmi": normalized_mutual_info_score(yt, yp, average_method="arithmetic"),
        "macro_f1": float(np.mean([r["f1"] for r in per_target])),
        "n_groups": m, "n_targets": n, "n_cells": len(cells),
        "unmatched_targets": [t for t in targets if matches[t] is None],
        "extra_groups": sorted(set(by_group) - used),
        "fragmented_targets": sum(k > 1 for k in fragment_count.values()),
        "mixed_groups": sum(k > 1 for k in mix.values()),
        "per_target": per_target,
    }


def _candidates(run, gt):
    with (run / "candidates.jsonl").open(encoding="utf-8") as f:
        meta = _indexed((json.loads(line) for line in f if line.strip()), "id")
    members = {candidate: set() for candidate in meta}
    for row in _read(run / "candidate_members.tsv", {"candidate_id", "cell_id"}):
        cid, cell = row["candidate_id"], row["cell_id"]
        if cid not in members or cell not in gt or cell in members[cid]:
            raise ValueError("Unknown or duplicate candidate membership: {} {}".format(cid, cell))
        members[cid].add(cell)
    if any(not isinstance(c.get("initial"), bool) for c in meta.values()):
        raise ValueError("Every frozen candidate needs an explicit initial flag")
    for cid, info in meta.items():
        if "set_hash" in info and info["set_hash"] != fingerprint(sorted(members[cid])):
            raise ValueError("Candidate membership-version mismatch: {}".format(cid))
    return meta, members


def _best(gt, members, candidates):
    sizes = Counter(gt.values())
    best = {t: {"candidate_id": None, "f1": 0.0, "tp": 0} for t in sizes}
    for cid in sorted(candidates):
        counts = Counter(gt[cell] for cell in members[cid])
        for target, tp in counts.items():
            f1 = 2 * tp / (len(members[cid]) + sizes[target])
            if f1 > best[target]["f1"]:
                best[target] = {"candidate_id": cid, "f1": f1, "tp": tp}
    return best


def evaluate(run: Path, gt: Path, out: Path):
    run, gt, out = map(lambda p: Path(p).resolve(), (run, gt, out))
    if _inside(gt, run) or _inside(out, run) or _inside(run, out):
        raise ValueError("Ground truth and evaluation output must be outside the discovery run")
    frozen, digest = _verified(run)
    gt_hash = _sha(gt)
    truth = read_truth(gt)
    part = _indexed(_read(run / "partition.tsv", {"cell_id", "group_id", "status"}), "cell_id")
    final = _partition_scores(truth, {c: r["group_id"] for c, r in part.items()})
    meta, members = _candidates(run, truth)
    for cell, row in part.items():
        owner = row.get("candidate_id")
        if owner and (owner not in members or cell not in members[owner]):
            raise ValueError("Final assignment contradicts submitted candidate version")
    initial_ids = {cid for cid in meta if meta[cid]["initial"]}
    best_initial, best = _best(truth, members, initial_ids), _best(truth, members, meta)
    matched = {r["target"]: r["group_id"] for r in final["per_target"]}
    recovered = []
    diagnostics = []
    matched_tp = defaultdict(set)
    for cell, label in truth.items():
        if part[cell]["group_id"] == matched[label]:
            matched_tp[label].add(cell)
    for target in sorted(matched):
        initial_id = best_initial[target]["candidate_id"]
        original = {c for c in members.get(initial_id, ()) if truth[c] == target}
        current = matched_tp[target]
        gained = sorted(current - original)
        for cell in gained:
            recovered.append({"target": target, "cell_id": cell,
                              "initial_comparator": initial_id, "final_group": matched[target]})
        diagnostics.append({"target": target, "best_candidate_id": best[target]["candidate_id"],
                            "best_single_candidate_f1": best[target]["f1"],
                            "best_initial_candidate_id": initial_id,
                            "best_initial_f1": best_initial[target]["f1"],
                            "new_true_members_vs_initial": len(gained),
                            "lost_true_members_vs_initial": len(original - current)})
    # Show actual membership gains/losses of each revision against its declared parents.
    revisions = []
    for cid, c in sorted(meta.items()):
        if c["initial"]:
            continue
        parents = c.get("parents", [])
        if any(pid not in members for pid in parents):
            raise ValueError("Frozen revision has missing parent")
        before = set().union(*(members[p] for p in parents))
        for cell in sorted(members[cid] ^ before):
            revisions.append({"candidate_id": cid, "cell_id": cell, "target": truth[cell],
                              "operation": "add" if cell in members[cid] else "remove",
                              "submitted_candidate_for_cell": part[cell].get("candidate_id") == cid})
    baseline_rows = defaultdict(dict)
    for row in _read(run / "artifacts/baselines.tsv", {"view", "resolution", "cell_id", "group_id"}):
        group = baseline_rows[row["view"], row["resolution"]]
        if row["cell_id"] in group:
            raise ValueError("Duplicate cell in baseline")
        group[row["cell_id"]] = row["group_id"]
    baselines = [{"view": v, "resolution": r, **_partition_scores(truth, pred)}
                 for (v, r), pred in sorted(baseline_rows.items())]
    status_counts = Counter(row["status"] for row in part.values())
    summary = {
        "schema_version": "evaluation.v1", "freeze_sha256": digest,
        "snapshot_hash": frozen["snapshot_hash"], "partition_sha256": frozen["files"]["partition.tsv"],
        "gt_sha256": gt_hash, "submitted_partition": final,
        "status_counts": dict(sorted(status_counts.items())),
        "unresolved_fraction": status_counts["unresolved"] / len(part),
        "n_candidates": len(meta), "n_initial_candidates": len(initial_ids),
        "offline_single_candidate_diagnostics": diagnostics, "baselines": baselines,
        "usage": _json(run / "usage.json"), "run_summary": _json(run / "run_summary.json"),
    }
    tables = ["targets.tsv", "candidate_diagnostics.tsv", "new_members.tsv", "revision_changes.tsv", "metrics.json"]
    if any((out / name).resolve() == gt for name in tables):
        raise ValueError("Evaluation output would overwrite the GT table")
    run_snapshot = summary["run_summary"].get("snapshot_hash")
    if run_snapshot is not None and run_snapshot != frozen["snapshot_hash"]:
        raise ValueError("Run and freeze disagree on prepare snapshot")
    if _verified(run)[1] != digest or _sha(gt) != gt_hash:
        raise ValueError("Frozen run or GT changed during scoring")
    out.parent.mkdir(parents=True, exist_ok=True)
    # Stage all five files together; a new result directory appears only on success.
    stage = Path(tempfile.mkdtemp(dir=out.parent, prefix=".evaluation-"))
    try:
        _table(stage / "targets.tsv", final["per_target"], list(final["per_target"][0]))
        _table(stage / "candidate_diagnostics.tsv", diagnostics, list(diagnostics[0]))
        _table(stage / "new_members.tsv", recovered, ["target", "cell_id", "initial_comparator", "final_group"])
        _table(stage / "revision_changes.tsv", revisions,
               ["candidate_id", "cell_id", "target", "operation", "submitted_candidate_for_cell"])
        _json_write(stage / "metrics.json", summary)
        try:
            # Atomic installation also works over an empty directory on POSIX.
            # A competing writer cannot replace a populated result directory.
            os.rename(stage, out)
        except OSError as exc:
            if exc.errno not in (errno.EEXIST, errno.ENOTEMPTY):
                raise
            for name in tables:
                existing = out / name
                if (not existing.is_file() or not _inside(existing.resolve(), out)
                        or _sha(existing) != _sha(stage / name)):
                    raise ValueError("Use a fresh output directory for a different or incomplete evaluation") from exc
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return summary


def aggregate(plan: Path, out: Path):
    plan, out = Path(plan).resolve(), Path(out).resolve()
    spec = yaml.safe_load(plan.read_text(encoding="utf-8"))
    entries = spec["replications"]
    if not isinstance(entries, list) or not entries:
        raise ValueError("Plan needs a nonempty replications list")
    ids, paths, completed_runs = set(), set(), set()
    records, bucket = [], defaultdict(list)
    arm_runs = set()
    for entry in entries:
        rid, arm = str(entry["replication_id"]), str(entry["arm"])
        path = (plan.parent / entry["metrics"]).resolve()
        if rid in ids or path in paths:
            raise ValueError("Replication IDs and metric files must be unique")
        ids.add(rid)
        paths.add(path)
        if not path.is_file():
            records.append({"replication_id": rid, "arm": arm, "status": "missing"})
            continue
        m = _json(path)
        if m.get("schema_version") != "evaluation.v1":
            raise ValueError("Unrecognized metrics in {}".format(path))
        identity = (m["freeze_sha256"], m["gt_sha256"])
        if identity in completed_runs:
            raise ValueError("Duplicate frozen run and GT in repetition plan")
        completed_runs.add(identity)
        run_info = m["run_summary"]
        if run_info.get("replication_id") not in (None, rid):
            raise ValueError("Predeclared replication_id differs from the recorded run")
        run_id = run_info.get("run_id")
        if run_id is not None:
            if (arm, run_id) in arm_runs:
                raise ValueError("The same discovery run was already counted in this arm")
            arm_runs.add((arm, run_id))
        run_status = run_info.get("status")
        values = {name: m["submitted_partition"][name] for name in ("ari", "nmi", "macro_f1")}
        for row in m["submitted_partition"]["per_target"]:
            values["target_f1/" + row["target"]] = row["f1"]
        for name, value in values.items():
            bucket[arm, name].append(value)
        records.append({"replication_id": rid, "arm": arm, "status": "scored",
                        "run_status": run_status, "run_id": run_id,
                        "recorded_replication_id": run_info.get("replication_id"),
                        "freeze_sha256": identity[0], "values": values})
    stats = [{"arm": arm, "metric": metric, "n": len(v), "mean": float(np.mean(v)),
              "median": float(np.median(v)), "min": min(v), "max": max(v),
              "std_sample": float(np.std(v, ddof=1)) if len(v) > 1 else None}
             for (arm, metric), v in sorted(bucket.items())]
    result = {"plan_sha256": _sha(plan), "n_planned": len(records),
              "n_scored": len(completed_runs), "replications": records, "distributions": stats}
    _json_write(out, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    command = sub.add_parser("aggregate")
    command.add_argument("--plan", required=True, type=Path)
    command.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    aggregate(args.plan, args.out)


if __name__ == "__main__":
    main()
