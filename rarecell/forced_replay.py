"""Repeat forced adjudication on original frozen claims using a content-verified review book."""
import shutil
import tempfile
from pathlib import Path
from types import SimpleNamespace

from scipy import sparse

from .config import fingerprint, hash_file, read_json, read_lines, write_json
from .data import read_dataset
from .discovery import rows
from .finalize import jsonlines, tsv
from .forced_routing import adjudicate
from .prepare import verify_snapshot


def replay(source, snapshot, dataset, book_path, out, cfg):
    source, out = Path(source).resolve(), Path(out).resolve()
    if out.exists():
        raise ValueError("Choose a new output directory")
    summary = read_json(source / "run_summary.json")
    if summary["mode"] != "A3_routed" or not summary["submitted"]:
        raise ValueError("Expected a submitted routed source run")
    book = read_json(book_path)
    if (book["snapshot_hash"] != summary["snapshot_hash"] or book["run_id"] != summary["run_id"]
            or book["route_claims_sha256"] != hash_file(source / "route_claims.jsonl")):
        raise ValueError("Review book does not match this frozen source")
    from .evaluation import _candidates, _verified
    if (source / "frozen.json").exists():
        frozen, source_digest = _verified(source)
        if not {"route_claims.jsonl", "candidate_members.tsv", "candidates.jsonl"} <= frozen["files"].keys():
            raise ValueError("Routed source inputs are missing from its freeze")
    else:
        source_digest = None
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out.parent) as temp:
        stage = Path(temp) / "submission"
        stage.mkdir()
        if snapshot:
            snap = Path(snapshot).resolve()
            manifest = verify_snapshot(snap)
            if manifest["snapshot_hash"] != summary["snapshot_hash"]:
                raise ValueError("Wrong prepare snapshot")
            baseline = snap / "baselines.tsv"
            expression_source = "verified_snapshot"
        else:
            # For exported runs where only the exact pre-prepare matrix is available.
            import yaml
            path = Path(dataset).resolve()
            path = path / "dataset.yaml" if path.is_dir() else path
            spec = yaml.safe_load(path.read_text())
            inputs = {key: hash_file(path.parent / spec[key]) for key in ("counts", "genes", "cells")}
            inputs["dataset"] = hash_file(path)
            if book.get("input_sha256") != inputs:
                raise ValueError("Review book is not bound to these dataset files")
            ds = read_dataset(path, book["target_sum"])
            snap = Path(temp) / "expression"
            snap.mkdir()
            sparse.save_npz(snap / "data.npz", ds.x, compressed=False)
            tsv(snap / "cells.tsv", ds.cells[0].keys(), ds.cells)
            tsv(snap / "genes.tsv", ds.genes[0].keys(), ds.genes)
            manifest = dict(ds.info, snapshot_hash=summary["snapshot_hash"])
            del ds
            baseline = source / "artifacts/baselines.tsv"
            expression_source = "verified_input_export"
        cells = [c["cell_id"] for c in rows(snap / "cells.tsv")]
        meta, members = _candidates(source, set(cells))
        claims = [dict(row, members=sorted(members[row["candidate_id"]]))
                for row in read_lines(source / "route_claims.jsonl") if row.get("candidate_id") and row["kind"] == "lineage"]
        ctrl = SimpleNamespace(snapshot=snap, snapshot_hash=summary["snapshot_hash"], manifest=manifest,
                             genes=rows(snap / "genes.tsv"), cells=cells, budget=cfg["budget"])
        returned = {}
        def answer(ctx, validate):
            phase = ctx["phase"]
            if phase == "conflict_label_audit":
                reply = {"reviews": {k: book["reviews"][k] for k in ctx["candidates"]}}
            elif phase == "conflict_panel_nomination":
                reply = {"labels": {k: book["labels"][k] for k in ctx["labels"]}}
            else:
                returned.update(ctx["pairs"])
                reply = {"rules": {k: book["rules"][k] for k in ctx["pairs"]}}
            validate(reply)
            return reply
        result = adjudicate(ctrl, claims, answer)
        if (fingerprint(returned) != book["measurements_sha256"]
                or fingerprint(result["frozen_claims"]) != book["members_sha256"]):
            raise ValueError("Measurements or original memberships changed since model review")
        assignments = {cell: group["id"] for group in result["groups"] for cell in group["members"]}
        if fingerprint(assignments) != book["assignments_sha256"]:
            raise ValueError("Final owners differ from the archived adjudication")
        artifacts = stage / "artifacts"
        artifacts.mkdir()
        write_json(artifacts / "forced_choice.json", result)
        original = {r["cell_id"]: r for r in rows(source / "partition.tsv")}
        if set(original) != set(cells):
            raise ValueError("Source partition and expression cells differ")
        owner = {}
        for group in result["groups"]:
            assigned = group["members"]
            if not assigned:
                continue
            cid = "R" + fingerprint([result["method"], group["id"], assigned])[:24]
            c = dict(id=cid, proposed_label=group["label"], parents=group["parents"], initial=False,
                   set_hash=fingerprint(assigned), sources=[{"kind": "forced_single_choice"}])
            meta[cid] = c
            members[cid] = set(assigned)
            for cell in assigned:
                if cell in owner:
                    raise ValueError("Nonexclusive replay partition")
                owner[cell] = cid
        if len(owner) != result["n_claimed"]:
            raise ValueError("Incomplete forced projection")
        partition = []
        for cell in cells:
            if cell in owner:
                cid = owner[cell]
                partition.append(dict(cell_id=cell, group_id="claim/" + cid, candidate_id=cid,
                                      proposed_label=meta[cid]["proposed_label"], status="provisional",
                                      assignment_reason="forced_single_choice"))
            else:
                row = original[cell]
                if row["assignment_reason"] != "unclaimed":
                    raise ValueError("Originally unclaimed universe inconsistent with frozen nominations")
                partition.append(row)
        tsv(stage / "partition.tsv", partition[0].keys(), partition)
        jsonlines(stage / "candidates.jsonl", (meta[cid] for cid in sorted(meta)))
        tsv(stage / "candidate_members.tsv", ["candidate_id", "cell_id"],
            (dict(candidate_id=cid, cell_id=c) for cid in sorted(members) for c in sorted(members[cid])))
        if baseline.is_file():
            shutil.copyfile(baseline, artifacts / "baselines.tsv")
        else:
            tsv(artifacts / "baselines.tsv", ["view", "resolution", "cell_id", "group_id"], [])
        write_json(stage / "usage.json", dict(replay_provider_calls=0, review_source=book.get("review_source", "archived")))
        report = dict(mode="A3_routed", status="complete", detail="forced_reconciliation_replay",
                    source_run=summary["run_id"], source_freeze_sha256=source_digest,
                    expression_source=expression_source, full_baselines_available=baseline.is_file(),
                    prepare_status=summary["prepare_status"], annotation_status=summary["annotation_status"],
                    submitted=True, snapshot_hash=summary["snapshot_hash"], n_entries=len(partition),
                    n_claimed=len(owner), n_conflicted=0, n_unclaimed=len(cells) - len(owner),
                    review_book_sha256=hash_file(book_path))
        write_json(stage / "run_summary.json", report)
        files = {p.relative_to(stage).as_posix(): hash_file(p) for p in sorted(stage.rglob("*")) if p.is_file()}
        write_json(stage / "frozen.json", dict(version=1, snapshot_hash=summary["snapshot_hash"], files=files))
        stage.rename(out)
    return report
