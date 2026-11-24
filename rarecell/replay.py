"""Replay frozen adaptive decisions through the current submit and assignment paths."""
from pathlib import Path
from shutil import copyfile
import tempfile

from . import claim_resolution
from .config import ADAPTIVE_MODE, config, hash_file, read_json, read_lines, write_json
from .discovery import Controller, rows
from .evaluation import _candidates, _partition_scores, _verified, read_truth
from .finalize import assignments, tsv


def replay_submit(run, out, preferred=(), candidate_ids=None, gt=None):
    run, out = Path(run).resolve(), Path(out).resolve()
    if out.exists() or run in out.parents:
        raise ValueError("Use a fresh replay output outside the frozen run")
    frozen, freeze_hash = _verified(run)
    info = read_json(run / "run_summary.json")
    if info["mode"] != ADAPTIVE_MODE or not info["submitted"]:
        raise ValueError("Replay requires a submitted A3_adaptive run")
    needed = {"decisions.jsonl", "annotations.tsv", "manifest.json"}
    if not needed <= frozen["files"].keys():
        raise ValueError("The decision history, annotations and manifest must be frozen")
    decisions = list(read_lines(run / "decisions.jsonl"))
    if not decisions or decisions[-1]["action"] != "submit":
        raise ValueError("Missing final frozen submit decision")
    # Replay only published data; the source DB and live credentials are unnecessary.
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out.parent) as tmp:
        root = Path(tmp)
        snapshot = root / "snapshot"
        snapshot.mkdir()
        manifest = read_json(run / "manifest.json")
        if manifest["snapshot_hash"] != frozen["snapshot_hash"]:
            raise ValueError("Replay manifest differs from the frozen submission")
        write_json(snapshot / "manifest.json", manifest)
        # Only the published cell universe and baseline are needed for final assignment.
        tsv(snapshot / "cells.tsv", ["cell_id"], rows(run / "partition.tsv"))
        tsv(snapshot / "genes.tsv", ["gene_id"], [])
        copyfile(run / "artifacts/baselines.tsv", snapshot / "baselines.tsv")
        ctrl = Controller(snapshot, root/"state", root/"cache", config(override={"mode": ADAPTIVE_MODE}), ADAPTIVE_MODE)
        try:
            meta, members = _candidates(run, ctrl.cell_set)
            for cid, c in meta.items():
                ctrl._index(dict(c, members=sorted(members[cid])))
            for a in rows(run / "annotations.tsv"):
                a["status"] = a.get("status") or "annotated"
                a["selection_history"] = {"genes": [g for g in a["genes"].split(",") if g]}
                ctrl.labels.setdefault(a["candidate_id"], {})[a["evidence_id"]] = a
            # The export retains insertion order, including multiple edits in the final turn.
            last = decisions[-1]
            recorded = last.get("claim_plan_ref")
            if recorded:
                if recorded not in frozen["files"]:
                    raise ValueError("Final decision references an unfrozen assignment plan")
                ctrl.state.set_meta("claim_plan_ref", ctrl.state.artifact(read_json(run / recorded)))
            else:
                for decision in decisions[:-1]:
                    if "selected" in decision:
                        claim_resolution.commit(ctrl, decision["selected"], decision.get("preferred_candidate_ids", ()))
                        ctrl.state.set_meta("selected", decision["selected"])
            ids = last["selected"] if candidate_ids is None else list(candidate_ids)
            preference = list(preferred) if preferred else [cid for cid in last.get("preferred_candidate_ids", ()) if cid in ids]
            receipt = ctrl.decide("submit", {"candidate_ids": ids, "preferred_candidate_ids": preference,
                                  "reason": "Offline replay of frozen candidate memberships"}, None, source="replay")
            partition, conflicts, groups, _ = assignments(ctrl)
            result = root/"result"
            result.mkdir()
            tsv(result/"partition.tsv", ["cell_id", "group_id", "candidate_id", "proposed_label", "status", "assignment_reason"], partition)
            report = {"source_freeze_sha256": freeze_hash, "source_decision": last["id"],
                      "historical_selected": last["selected"], "requested": ids, "preferred": preference,
                      "resolution": receipt["claim_resolution"], "unresolved_conflicts": len(conflicts),
                      "claim_groups": [g for g in groups if g["group_id"].startswith("claim/")],
                      "partition_sha256": hash_file(result/"partition.tsv")}
            if gt:
                # Read reference labels only after assignment and hashing.
                truth = read_truth(Path(gt))
                score = _partition_scores(truth, {r["cell_id"]: r["group_id"] for r in partition})
                tsv(result/"targets.tsv", list(score["per_target"][0]), score["per_target"])
            write_json(result/"replay.json", report)
            result.rename(out)
            return report
        finally:
            ctrl.state.close()
