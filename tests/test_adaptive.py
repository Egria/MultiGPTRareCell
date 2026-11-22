"""Regression for the recorded singleton request and its downstream assignment."""
import csv
import json
from pathlib import Path
from contextlib import closing

import pytest

from rarecell import job_queue
from rarecell.cli import parser
from rarecell.config import ADAPTIVE_MODE, config
from rarecell.discovery import Controller
from rarecell.finalize import assignments
from rarecell.prepare import prepare
from rarecell.protocol import parse_message
from rarecell.run import explore, process_message
from test_science import expression_fixture


def act(ctrl, name, **params):
    raw = json.dumps(dict(action=name, **params))
    parsed, = parse_message(raw)
    return ctrl.execute(parsed, "test")


def test_permissions_scope_reuse_and_exact_rules(tmp_path):
    snap = expression_fixture(tmp_path / "fixture", 4)
    with (snap / "cells.tsv").open() as stream:
        cells = [r["cell_id"] for r in csv.DictReader(stream, delimiter="\t")]
    with (snap.parent / "discovery_input/cells.tsv").open() as stream:
        seed_cell = next(csv.DictReader(stream, delimiter="\t"))["cell_id"]
    objects = [{"id": cid, "members": group, "sources": [], "markers": []}
               for cid, group in [("seed", [seed_cell]), ("other", cells[1:2]), ("broad", cells)]]
    (snap / "initial.jsonl").write_text("".join(json.dumps(row) + "\n" for row in objects))
    (snap / "seeds.jsonl").write_text(json.dumps({"id": "S", "cell_id": seed_cell, "candidate_id": "seed"}) + "\n")
    manifest = json.loads((snap / "manifest.json").read_text())
    manifest.update(snapshot_hash="test-snapshot", data_version="test-data")
    (snap / "manifest.json").write_text(json.dumps(manifest))
    cfg = config(override={"search":{"max_cells": 2}, "budget":{"max_expansion_depth": 2}})
    for mode in ("A3", ADAPTIVE_MODE):
        ctrl = Controller(snap, tmp_path / mode, tmp_path / "cache", cfg, mode)
        with closing(ctrl.state):
            ctrl.initialize()
            first = act(ctrl, "expand_candidate", candidate_id="seed", background="all")
            assert first["background_rule"] == "multiview_local"
            repeat = act(ctrl, "expand_candidate", candidate_id="seed")
            assert first["query_key"] == repeat["query_key"] and ctrl.state.usage()["compute_executions"] == 1
            with pytest.raises(ValueError):
                act(ctrl, "expand_candidate", candidate_id="seed", background=["all"])
            complement = act(ctrl, "expand_candidate", candidate_id="seed", background_mode="complement")
            assert complement["query_key"] != first["query_key"]
            ev = first["evidence_id"]
            args = {"candidate_id": "seed", "evidence_ids":[ev], "metric":"score", "threshold":0, "reason":"Compare all scored cells"}
            if mode != ADAPTIVE_MODE:
                assert "apply_boundary" not in ctrl.context(1)["actions"]
                with pytest.raises(ValueError):
                    act(ctrl, "apply_boundary", **args)
                with pytest.raises(ValueError):
                    act(ctrl, "coexpress", scope="all", genes=["g020","g021"], min_detected=2)
                continue
            assert ctrl.context(1)["rarity_search"]["resolvable_global_bands_pct"] == [5]
            too_many = act(ctrl, "apply_boundary", **args)
            assert too_many["status"] == "over_limit" and too_many["added_count"] > 2 and "candidate_id" not in too_many
            with pytest.raises(ValueError, match="original candidate"):
                act(ctrl, "apply_boundary", **dict(args, candidate_id="other"))
            ids = [ev]
            page = act(ctrl, "get_evidence", evidence_ids=ids, metric="positive_expression", offset=2, limit=3)
            ranked = page["sections"]["ranking"]
            assert ranked["total"] == len(cells) and len(ranked["rows"]) == 3
            panel = ["g020", "g021", "g006"]
            both = act(ctrl, "coexpress", scope="all", genes=panel, min_detected=3)
            same = act(ctrl, "coexpress", scope="all", genes=panel, logic="AND")
            assert same["evidence_id"] == both["evidence_id"]
            k = act(ctrl, "coexpress", scope="all", genes=panel, min_detected=2)
            assert k["summary"]["logic"] == "k_of_m" and k["summary"]["n_matching"] >= both["summary"]["n_matching"]
            with pytest.raises(ValueError, match="distinct resolved"):
                act(ctrl, "coexpress", scope="all", genes=panel, min_detected=4)
            with pytest.raises(ValueError, match="coexpress"):
                act(ctrl, "propose_revision", groups=[{"matched_evidence_id": ev}], parents=["seed"],
                    evidence_ids=ids, reason="invalid source")


def test_gui_cli_and_slurm_entry_point_share_mode(tmp_path):
    root = tmp_path / "queue"
    job_queue.initialize(root)
    data = Path(__file__).resolve().parents[1] / "examples/minimal_sparse/input/dataset.yaml"
    task_id = job_queue.submit(root, tmp_path, {"kind":"pipeline", "mode": ADAPTIVE_MODE,
                              "dataset":str(data), "model":"test-deployment", "annotation_model":"test-deployment"})
    args = list(job_queue.steps(root, job_queue.get_task(root, task_id)))
    assert [a[0] for a in args] == ["prepare", "annotate", "discover"]
    parsed = parser().parse_args(args[-1])
    assert parsed.mode == config(parsed.config)["mode"] == ADAPTIVE_MODE and parsed.annotations
    a = parser().parse_args(["ablate", "--snapshot", str(tmp_path), "--out", str(tmp_path/"out"),
                             "--arms", "A3", ADAPTIVE_MODE])
    assert a.arms == ["A3", ADAPTIVE_MODE]


@pytest.fixture(scope="module")
def three_views(tmp_path_factory):
    # Generate the real three-view snapshot. A Python-only installation can run the other tests.
    import shutil
    import subprocess
    rscript = shutil.which("Rscript")
    if not rscript or subprocess.run([rscript, "--vanilla", "-e",
            'quit(status=as.integer(!all(sapply(c("Seurat","SeuratObject","leidenbase"), requireNamespace, quietly=TRUE))))'],
            capture_output=True).returncode:
        pytest.skip("requires external Seurat and leidenbase for the three-view fixture")
    root = Path(__file__).resolve().parents[1] / "examples/minimal_sparse"
    snap = tmp_path_factory.mktemp("adaptive_views") / "snapshot"
    prepare(root / "input/dataset.yaml", snap, config(override={"mode": ADAPTIVE_MODE}))
    return snap


def test_minimal_sparse_adaptive_seed_to_final_partition(tmp_path, three_views):
    root = Path(__file__).resolve().parents[1] / "examples/minimal_sparse"
    snap = three_views
    cfg = config(override={"mode": ADAPTIVE_MODE})
    ctrl = Controller(snap, tmp_path / "run", tmp_path / "cache", cfg, ADAPTIVE_MODE)
    seed = "H235dc8d4ae01746f2bfffbb0"  # exact anchor from the reported T0003
    with closing(ctrl.state):
        ctrl.initialize()
        ctrl.state.set_meta("run_id", "run")
        raw = json.dumps({"action":"expand_candidate", "candidate_id":seed, "background":"all",
                          "evidence_ids":ctrl.candidate(seed)["evidence_ids"], "reason":"Search for outlier neighbors."})
        msg = {"message_id":"run:T0001", "raw":raw, "complete":True}
        process_message(ctrl, msg, 1)
        expansion = ctrl.state.get("Receipt","run:T0001/0")["result"]
        assert expansion["added_count"] == 2 and len(ctrl.candidate(expansion["candidate_id"])["members"]) == 3
        ev = expansion["evidence_id"]
        ranked = act(ctrl, "get_evidence", evidence_ids=[ev])["sections"]["ranking"]["rows"]
        assert all(0 <= row[metric] <= 1 for row in ranked
                   for metric in ("score", "expression", "positive_expression") if row[metric] is not None)
        gap = expansion["summary"]["rank_gaps"][0]
        cutoff = (gap["upper"] + gap["lower"]) / 2
        args = dict(candidate_id=seed, metric="score", threshold=cutoff, evidence_ids=[ev], reason="Use largest rank gap")
        boundary = act(ctrl, "apply_boundary", **args)
        small = boundary["candidate_id"]
        assert act(ctrl, "apply_boundary", **args)["candidate_id"] == small
        assert ctrl.state.usage()["new_queries"] == 1
        assert ctrl.candidate(small)["boundary_rule"]["threshold"] == cutoff
        broad = [c["id"] for c in ctrl.ordered() if c["initial"] and len(c["members"]) == 80]
        rejected = act(ctrl, "preview_submission", candidate_ids=broad, candidate_id=small, reference_id=broad[0])
        assert rejected["candidates_reviewed"][0]["assigned_to_other"] == 8
        preview = act(ctrl, "preview_submission", candidate_ids=[small,*broad], candidate_id=small, reference_id=broad[0])
        item = preview["candidates_reviewed"][0]
        assert preview["conflicted_cells"] == 8 and item["retained"] == 0
        assert item["proposed_local"]["n"] == 4 and item["proposed_local"]["denominator"] == 80
        co = act(ctrl, "coexpress", genes=expansion["summary"]["positive_genes"], scope="all", logic="AND")
        groups = [{"matched_evidence_id":co["evidence_id"], "status":"provisional",
                   "biological_rarity":"uncertain", "rarity_reason":"Synthetic context; no population prior"}]
        groups.extend({"candidate_id":cid, "exclude_candidate_ids":[small]} for cid in broad)
        rev = act(ctrl, "propose_revision", groups=groups, parents=[small,*broad], evidence_ids=[ev], reason="Separate the shared program")
        assert len(rev["candidate_ids"]) == 3
        done = act(ctrl, "submit", reason="Accept trial memberships")
        decision = ctrl.state.get("Decision",done["decision_id"])
        assert decision["assignment_preview"]["exclusive_cells"] == 160
        partition, conflicts, _, _ = assignments(ctrl)
        assert not conflicts and all(row["assignment_reason"] == "selected" for row in partition)
        # Only now open truth and check membership, including both coarse backgrounds.
        with (root / "evaluation_fixture/gt.tsv").open() as stream:
            truth = {r["cell_id"] for r in csv.DictReader(stream, delimiter="\t") if r["label"] == "synthetic_rare"}
        assert set(ctrl.candidate(small)["members"]) == truth
        rare_group, = [ctrl.candidate(cid) for cid in done["candidate_ids"] if set(ctrl.candidate(cid)["members"]) == truth]
        assert rare_group["biological_rarity"]["assessment"] == "uncertain" and rare_group["status"] == "provisional"

@pytest.mark.parametrize("ending,expected_n,source", [
    ("wider", 8, "agent"), ("core", 3, "agent"), ("decline", 0, "agent"),
    ("repeat", 3, "controller_stall"), ("no_turn", 3, "controller_stall"),
    ("no_tokens", 3, "controller_stall"),
])
def test_adaptive_review_and_stall_submission(tmp_path, three_views, ending, expected_n, source):
    cfg = config(override={"mode": ADAPTIVE_MODE, "budget": {"max_agent_turns": 6 if ending == "no_turn" else 12}})
    run, cache = tmp_path / "run", tmp_path / "cache"
    ctrl = Controller(three_views, run, cache, cfg, ADAPTIVE_MODE)
    anchor = "H235dc8d4ae01746f2bfffbb0"
    def save(turn, command):
        mid = f"audit:T{turn:04}"
        ctrl.state.put("Message", mid, {"message_id": mid, "raw": json.dumps(command), "complete": True,
                                        "usage": {"provider": "transcript_replay", "total_tokens": 1}})
    with closing(ctrl.state):
        ctrl.initialize()
        ctrl.state.set_meta("run_id", "audit")
        save(1, {"action": "expand_candidate", "candidate_id": anchor})
        assert explore(ctrl, "live", pause_after=1)[1] == "paused"
        eid = ctrl.state.get("Receipt", "audit:T0001/0")["result"]["evidence_id"]
        read = {"action": "get_evidence", "evidence_ids": [eid], "offset": 0, "limit": 20}
        cut = {"action": "apply_boundary", "candidate_id": anchor, "evidence_ids": [eid],
               "metric": "score", "threshold": 0.82, "reason": "Keep the measured core."}
        for turn, command in enumerate([read, cut, read, cut, read], 2):
            save(turn, command)
        status, phase = explore(ctrl, "live", pause_after=6)
        review = ctrl.context(7)["boundary_reviews"][0]
        core, wider = review["requested"], review["wider"]
        assert (core["n"], wider["n"], core["selected"], wider["selected"]) == (3, 8, True, False)
        assert review["cohorts"]["additional"]["scored"] == 5
        assert review["cohorts"]["additional"]["mean"]["positive_expression"] > 0.9
        assert [core["partition_unclaimed"], wider["partition_unclaimed"]] == [157, 152]
        if ending != "no_turn":
            assert phase == "paused" and set(ctrl.context(7)["actions"]) == {"submit", "propose_revision", "select_boundary"}
        assert ctrl.context(7)["recent_evidence_pages"][0]["offset"] == 0
        with pytest.raises(ValueError, match="one boundary"):
            act(ctrl, "submit", candidate_ids=[core["candidate_id"], wider["candidate_id"]])
    # Reopen SQLite so the draft and the comparison must survive a process boundary.
    ctrl = Controller(three_views, run, cache, cfg, ADAPTIVE_MODE)
    with closing(ctrl.state):
        if ending == "no_tokens":
            ctrl.budget["max_reported_tokens"] = 6
        elif ending != "no_turn":
            last = {"action": "submit", "reason": "Accept core."}
            if ending == "wider":
                last = [{"action": "select_boundary", "candidate_id": wider["candidate_id"], "reason": "Include margin."},
                        {"action": "submit", "reason": "Accept whole group."}]
            elif ending == "decline":
                last["candidate_ids"] = []
            elif ending == "repeat":
                last = read
            save(7, last)
        if ending != "no_turn":
            status, phase = explore(ctrl, "live")
        assert (status, phase) == (("complete", "submitted") if source == "agent" else ("partial", "stalled_provisional"))
        assert ctrl.state.get_meta("submission_source") == source
        decision = ctrl.state.all("Decision")[-1]
        if source == "controller_stall":
            assert decision["message_id"] is None and decision["status"] == "provisional"
        assert ctrl.state.get_meta("finished_turn") == (6 if ending in {"no_turn", "no_tokens"} else 7)
        assert not ctrl.state.all("ProviderDispatch") and ctrl.state.usage()["new_queries"] == 1
        partition, conflicts, _, _ = assignments(ctrl)
        owned = {r["cell_id"] for r in partition if r["assignment_reason"] == "selected"}
        gt = Path(__file__).resolve().parents[1] / "examples/minimal_sparse/evaluation_fixture/gt.tsv"
        with gt.open() as stream:
            truth = {r["cell_id"] for r in csv.DictReader(stream, delimiter="\t") if r["label"] == "synthetic_rare"}
        assert len(owned) == expected_n and owned <= truth and not conflicts and len(partition) == 160


def test_review_options_preserve_unrelated_selections_and_revise_held_wider(tmp_path, three_views):
    cfg = config(override={"mode": ADAPTIVE_MODE})
    ctrl = Controller(three_views, tmp_path/"run", tmp_path/"cache", cfg, ADAPTIVE_MODE)
    with closing(ctrl.state):
        ctrl.initialize()
        seed = "H235dc8d4ae01746f2bfffbb0"
        exp = act(ctrl, "expand_candidate", candidate_id=seed)
        eid = exp["evidence_id"]
        args = dict(candidate_id=seed, evidence_ids=[eid], metric="score", threshold=.82, reason="Core")
        review = act(ctrl, "apply_boundary", **args)["boundary_review"]
        small, wide = (review[k]["candidate_id"] for k in ("requested", "wider"))
        assert set(ctrl.candidate(small)["members"]) < set(ctrl.candidate(wide)["members"])
        outside = sorted(ctrl.cell_set - set(ctrl.candidate(wide)["members"]))
        parent = next(c["id"] for c in ctrl.ordered() if c["initial"] and outside[0] in c["members"])
        separate = act(ctrl, "propose_revision", parents=[parent], groups=[{"members":outside[:1]}],
                       reason="Unrelated provisional group", evidence_ids=[eid])
        unrelated, = set(separate["candidate_ids"]) - {small}
        act(ctrl, "select_boundary", candidate_id=wide, reason="Accept additional cells")
        assert set(ctrl.state.get_meta("selected")) == {wide, unrelated}
        act(ctrl, "apply_boundary", **args)
        assert set(ctrl.state.get_meta("selected")) == {wide, unrelated}
        act(ctrl, "select_boundary", candidate_id=small, reason="Retain core while reviewing")
        revised = act(ctrl, "propose_revision", parents=[wide], groups=[{"candidate_id":wide}],
                      reason="Accept the wider boundary after review", evidence_ids=[eid])
        assert unrelated in revised["candidate_ids"] and small not in revised["candidate_ids"]
        assert len(revised["candidate_ids"]) == 2
        assert ctrl.state.usage()["new_queries"] == 1
