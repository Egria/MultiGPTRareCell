"""Submit competing measured drafts through the same controller and exporter as a live run."""
from contextlib import closing
import json

import pytest

from rarecell import claim_resolution
from rarecell.config import ADAPTIVE_MODE, config, contract, hash_file, read_json, write_json
from rarecell.discovery import Controller, rows
from rarecell.finalize import assignments, export
from rarecell.replay import replay_submit
from test_adaptive import act
from test_science import expression_fixture


@pytest.fixture
def ctrl(tmp_path):
    snap = expression_fixture(tmp_path / "fixture", 4)
    write_json(snap / "manifest.json", dict(read_json(snap / "manifest.json"), snapshot_hash="snapshot", data_version="v"))
    (snap / "seeds.jsonl").touch()
    anchor = {"id": "anchor", "members": [rows(snap / "cells.tsv")[0]["cell_id"]], "sources": [], "markers": []}
    (snap / "initial.jsonl").write_text(json.dumps(anchor) + "\n")
    (snap / "baselines.tsv").write_text("view\tresolution\tcell_id\tgroup_id\n")
    run = tmp_path / "run"
    cfg = config(override={"mode": ADAPTIVE_MODE, "exhaustive_seed_review": True})
    c = Controller(snap, run, tmp_path / "cache", cfg, ADAPTIVE_MODE)
    with closing(c.state):
        c.initialize()
        write_json(run / "manifest.json", read_json(snap / "manifest.json"))
        write_json(run / "config_resolved.yaml", cfg)
        write_json(run / "project_contract.yaml", contract())
        yield c


def draft(c, label, members, why):
    return c.make_revision(members, ["anchor"], label, "provisional", [], "revision", why)["id"]


def test_earlier_viable_group_and_explicit_priority_survive_export_and_replay(ctrl, tmp_path):
    cells = list(ctrl.cells)
    wide = draft(ctrl, "epsilon", cells[:6], "wide")
    twin = draft(ctrl, "epsilon", cells[:6], "identical")
    near = draft(ctrl, "epsilon", cells[1:7], "similar")
    core = draft(ctrl, "epsilon", cells[1:4], "core")
    other = draft(ctrl, "other", cells[8:12], "separate")
    first = act(ctrl, "propose_revision", parents=[wide], evidence_ids=ctrl.candidate(wide)["evidence_ids"],
                groups=[{"candidate_id": wide, "proposed_label": "epsilon"}], reason="Register initial boundary")
    start, = first["candidate_ids"]
    all_ids = [start, near, twin, core, other]
    preview = act(ctrl, "preview_submission", candidate_ids=all_ids)
    checkout = act(ctrl, "submit", candidate_ids=all_ids, reason="Preserve viable draft")
    assert preview["claim_resolution"] == checkout["claim_resolution"]
    assert checkout["claim_resolution"]["retained"] == {start: 6, twin: 0, near: 0, core: 0, other: 4}
    assert {r["rule"] for r in checkout["claim_resolution"]["held"].values()} == {"alternative_boundary", "identical"}
    export(ctrl, "complete", "submitted", "replay", "success")
    run = ctrl.run
    original = hash_file(run / "partition.tsv")
    again = replay_submit(run, tmp_path / "again")
    assert again["partition_sha256"] == original
    # Legacy runs store the decision sequence and member sets without an assignment plan.
    history = [json.loads(line) for line in (run/"decisions.jsonl").read_text().splitlines()]
    for d in history:
        d.pop("claim_plan_ref", None)
        d.pop("claim_resolution", None)
    (run/"decisions.jsonl").write_text("".join(json.dumps(d)+"\n" for d in history))
    manifest = read_json(run/"frozen.json")
    manifest["files"]["decisions.jsonl"] = hash_file(run/"decisions.jsonl")
    write_json(run/"frozen.json", manifest)
    assert replay_submit(run, tmp_path/"legacy")["partition_sha256"] == original
    changed = replay_submit(run, tmp_path / "changed", preferred=[near, twin])
    assert changed["resolution"]["retained"][near] == 6
    assert {r["cell_id"] for r in rows(tmp_path/"changed/partition.tsv") if r["candidate_id"] == near} == set(cells[1:7])
    assert hash_file(run / "partition.tsv") == original


def test_competing_partial_boundaries_and_ambiguous_names(ctrl):
    cells = list(ctrl.cells)
    a = draft(ctrl, "a", cells[:5], "first")
    b = draft(ctrl, "b", cells[3:8], "second")
    rough = act(ctrl, "preview_submission", candidate_ids=[a, b])
    assert rough["conflicted_cells"] == 2
    p = act(ctrl, "preview_submission", candidate_ids=[a, b], preferred_candidate_ids=[b, a])
    assert p["conflicted_cells"] == 0 and p["claim_resolution"]["retained"] == {a: 3, b: 5}
    near = draft(ctrl, "different", cells[1:6], "other rationale")
    p = act(ctrl, "preview_submission", candidate_ids=[a, near])
    assert p["conflicted_cells"] == 0 and p["claim_resolution"]["provisional_overlap_cells"] == 4
    duplicate = draft(ctrl, "different", cells[:5], "same cells")
    act(ctrl, "submit", candidate_ids=[a, duplicate])
    partition, conflicts, _, _ = assignments(ctrl)
    assert not conflicts and {r["proposed_label"] for r in partition if r["candidate_id"]} == {"unknown"}
    act(ctrl, "submit", candidate_ids=[a, duplicate], preferred_candidate_ids=[a])
    assert {r["proposed_label"] for r in assignments(ctrl)[0] if r["candidate_id"]} == {"a"}
    # A held intermediate cannot join genuinely distinct cohorts through transitive overlap.
    close = draft(ctrl, "a", cells[2:7], "bridge")
    end = draft(ctrl, "a", cells[4:9], "end")
    plan = claim_resolution.plan(ctrl, [a, close, end], [a])
    assert close in plan["held"] and end in plan["order"]


def test_separate_pool_programs_bind_exact_global_members(ctrl):
    seed_gene, second = "g020", "g021"
    evidence = act(ctrl, "compare_groups", members=list(ctrl.cells[:3]), genes=[seed_gene])["evidence_id"]
    c = ctrl.make_revision(ctrl.cells[:3], ["anchor"], "glial", "provisional", [evidence], "revision", "pool")
    c["ontology_pool"] = {"marker_counts": [(seed_gene, 2)]}
    ctrl.state.put("PoolBatch", "test", {"candidate_ids": [c["id"]]})
    receipt = act(ctrl, "investigate_pool", candidate_id=c["id"], seed_genes=[seed_gene], support_genes=[second],
                  counter_genes=["g006"], min_detected=1, reason="Measure joint programs")
    assert receipt["status"] == "ok" and receipt["n_background"] == len(ctrl.cells)-3
    assert receipt["programs"]["support"]["genes"] == [second]
    assert receipt["programs"]["counter"]["curve"][0]["outside_hits"] == len(ctrl.cells)-3
    assert sum(receipt["programs"]["seed"]["curve"][0][k] for k in ("outside_hits", "seed_hits")) == 4
    measured = act(ctrl, "coexpress", scope="all", genes=[seed_gene])
    revised = act(ctrl, "propose_revision", parents=[c["id"]], evidence_ids=[receipt["evidence_id"]],
                  groups=[dict(matched_evidence_id=receipt["evidence_id"], program="seed", min_detected=1,
                               proposed_label="glial")], reason="Select exact measured group")
    cid, = revised["candidate_ids"]
    from rarecell.adaptive import coexpression_members
    assert set(ctrl.candidate(cid)["members"]) == set(coexpression_members(ctrl, measured["evidence_id"]))
