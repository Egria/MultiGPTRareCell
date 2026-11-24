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


def test_requested_pool_result_reaches_next_dispatch_before_catalogue(ctrl):
    from rarecell.config import json_bytes
    from rarecell.run import obtain, process_message
    from rarecell import providers
    from unittest.mock import patch

    state = ctrl.state
    row = {"candidate_id": "anchor", "label": "rare", "n_seeds": 3, "n_background": 77,
           "status": "ok", "evidence_id": ctrl.candidate("anchor")["evidence_ids"][0],
           "contrasts": {"g020": {"effect": 4., "high_seeds": 3, "high_outside": 1, "means_per_10k": [0.1, 5.]}},
           "programs": {"seed": {"genes": ["g020"], "curve": [{"k": 1, "seed_hits": 3, "outside_hits": 1}]}}}
    state.put("PoolFollowup", "anchor", row)
    state.set_meta("run_id", "test")
    # Force the full page above the context cap. A compressed numeric result must arrive.
    original = ctrl.read_tool
    def padded(action, args):
        return dict(original(action, args), full_payload="a" * (ctrl.budget["context_bytes"] + 1000))
    with patch.object(ctrl, "read_tool", padded):
        process_message(ctrl, {"message_id": "test:T0001", "raw": json.dumps(
            {"action": "get_evidence", "candidate_id": "anchor"})}, 1)
    rid, = state.get_meta("requested_receipts")
    for i in range(20):
        state.put("PoolFollowup", str(i), dict(row, candidate_id=str(i)))
    ctx = ctrl.context(2)
    delivered, = ctx["requested_evidence"]
    assert delivered["receipt_id"] == rid
    pool = delivered["result"]["pool_followup"]
    assert pool["evidence_id"] == row["evidence_id"]
    assert pool["contrasts"]["g020"]["means_per_10k"] == [0.1, 5.]
    assert pool["contrasts"]["g020"]["high_outside"] == 1
    assert pool["programs"]["seed"]["seed_hits_by_k"] == [3]
    assert len(json_bytes(ctx)) <= ctrl.budget["context_bytes"] and state.get_meta("requested_receipts") == [rid]
    with patch.object(providers, "diagnose", return_value={"status":"ok"}), patch.object(
            providers, "complete", return_value=('{"action":"submit","candidate_ids":[]}', {"api_calls":1})):
        message = obtain(ctrl, "live", 2, False)
    assert message["delivered_receipts"] == [rid]
    process_message(ctrl, message, 2)
    assert state.get_meta("requested_receipts") == []


def test_independent_disagreement_uses_exact_expansion_and_one_reserved_query(ctrl, tmp_path):
    from rarecell import pool_followup
    from rarecell.config import read_lines
    from rarecell.run import process_message
    state = ctrl.state
    cid = draft(ctrl, "unknown", ctrl.cells[:4], "seed expanded")
    for bg,symbol,cl in [("neighbor","g020","CL:0000001"),("global","g006","CL:0000002")]:
        eid = ctrl.evidence("anchor", "interpretation", ctrl.candidate("anchor")["evidence_ids"], "input_expression", ctrl.candidate("anchor")["set_hash"],
                            {"label": cl}, [])
        row = {"candidate_id":"anchor", "label":cl, "ontology_id":cl, "background_type":bg,
               "status":"annotated", "evidence_id":eid, "selection_history":{"genes":[symbol]}}
        state.put("Annotation",bg,row)
        ctrl.labels.setdefault("anchor",{})[eid] = row
    eid = ctrl.candidate("anchor")["evidence_ids"][0]
    state.put("SeedMeasurement","anchor",{"candidate_id":"anchor", "evidence_id":eid, "status":"expanded", "summary":{"positive_genes":["g020"]}})
    state.put("SeedAssessment","anchor",{"candidate_id":"anchor","seed_id":"S","message_id":"seed","measurement_status":"expanded",
                                        "proposal_id":cid,"reason":"Check the competing lineage","disposition":"follow_up"})
    ctrl.budget["max_new_queries"] = 0
    state.set_meta("run_id", "test")
    mid = "test:T0001"
    assert pool_followup.prepare_batch(ctrl, mid)
    item, = ctrl.context(1)["pools"]
    assert item["anchor_id"] == "anchor" and set(item["members"]) == set(ctrl.cells[:4])
    assert item["annotations"]["neighbor"]["label"] == "CL:0000001"
    cmd = json.dumps({"action":"investigate_pool","candidate_id":cid,"reason":"Resolve ambiguity",
                      "seed_genes":["g020"],"support_genes":["g021"],"counter_genes":["g006"]})
    process_message(ctrl, {"message_id": mid, "raw":cmd}, 1)
    assert state.usage()["new_queries"] == ctrl.query_limit == 1
    run = ctrl.run
    reopened = Controller(ctrl.snapshot,run,tmp_path/"cache",ctrl.cfg,ADAPTIVE_MODE)
    with closing(reopened.state):
        assert not pool_followup.prepare_batch(reopened,"test:T0002")
        cov = pool_followup.coverage(reopened)
        assert cov["conflict_measured"] == 1 and cov["eligible_pools"] == 0
        followup = reopened.read_tool("get_evidence",{"candidate_id":cid})["pool_followup"]
        assert followup["n_seeds"] == 4 and followup["n_background"] == len(ctrl.cells)-4
        assert followup["written_questions"][0]["reason"] == "Check the competing lineage"
        assert len(reopened.state.all("ConflictQuerySlot")) == 1
        assert reopened.context(2)["pool_followups"][0]["candidate_id"] == cid
        with pytest.raises(ValueError, match="max_new_queries"):
            act(reopened,"coexpress",scope="all",genes=["g006"])
        export(reopened,"partial","paused","live","success")
        assert next(read_lines(run/"pool_followups.jsonl"))["follow_up_seed_ids"] == ["anchor"]


def test_deferred_scientific_receipt_survives_small_context_and_reopen(ctrl, tmp_path):
    from rarecell import handoff
    from rarecell.run import process_message
    msg = {"message_id":"test:T0001","raw":json.dumps({"action":"compare_groups",
           "candidate_id":"anchor","genes":["g020", "g021"]})}
    process_message(ctrl,msg,1)
    rid, = ctrl.state.get_meta("requested_receipts")
    frame = {}
    assert not handoff.requested(frame,ctrl.state,75)
    assert frame["pending_evidence_count"] == 1
    reopened = Controller(ctrl.snapshot,ctrl.run,tmp_path/"cache",ctrl.cfg,ADAPTIVE_MODE)
    with closing(reopened.state):
        next_context = reopened.context(2)
        item, = next_context["requested_evidence"]
        assert item["receipt_id"] == rid and next_context["pending_evidence_count"] == 0
        assert item["result"]["evidence_id"] == reopened.state.get("Receipt",rid)["result"]["evidence_id"]
        assert set(item["result"]["nb_high_counts"]) == {"g020","g021"}


def test_named_seed_guides_search_and_reviews_exact_additions(ctrl, tmp_path):
    from rarecell import exhaustive, expansion_review, handoff, pool_followup
    from rarecell.run import process_message

    s = ctrl.state
    s.set_meta("run_id", "named")
    ctrl.budget["max_new_queries"] = 0
    seed = read_json(ctrl.snapshot.parent / "evaluation_fixture" / "truth.json")["seed"]
    anchor = draft(ctrl, "unknown", [seed], "singleton to audit")
    ctrl.seed_ids[anchor] = "seed"
    obs = ctrl.candidate(anchor)["evidence_ids"][0]
    ctrl.labels[anchor] = {obs: {"candidate_id": anchor, "label": "named lineage",
        "evidence_id": obs, "status": "annotated", "background_type": "global",
        "selection_history": {"genes": ["g020"]}}}
    s.put("SeedMeasurement", anchor, {"candidate_id": anchor, "evidence_id": obs, "status": "no_expansion"})
    s.put("SeedBatch", "assessment", {"candidate_ids": [anchor]})
    note = dict(candidate_id=anchor, disposition="retain", reason="The measured genes support this name")
    with pytest.raises(ValueError, match="label_reliability"):
        exhaustive.record(ctrl, "assessment", note)
    with pytest.raises(ValueError, match="proposed names"):
        exhaustive.record(ctrl, "assessment", dict(note, guide_label="invented", label_reliability="supported"))
    exhaustive.record(ctrl, "assessment", dict(note, guide_label="Named Lineage", label_reliability="supported"))
    assert pool_followup.prepare_batch(ctrl, "named:T0001")
    case, = ctrl.context(1)["pools"]
    assert case["kind"] == "guided_named_seed" and case["seed_name_assessment"]["label_reliability"] == "supported"
    def send(i, cmd):
        mid = f"named:T{i:04d}"
        process_message(ctrl, {"message_id": mid, "raw": json.dumps(cmd)}, i)
        return s.get("Receipt", mid + "/0")["result"]

    receipt = send(1, dict(action="investigate_pool", candidate_id=anchor, seed_genes=["g020"],
        support_genes=["g020", "g021"], counter_genes=["g006"], reason="Compare additional cells across the dataset"))
    assert receipt["status"] == "ok"
    option = next(b for b in receipt["boundaries"] if {"program": "seed", "min_detected": 1} in b["rules"])
    assert option["n_added"] == 3 and option["added_hits_by_k"]["counter"] == [3]
    assert handoff.index(receipt)["programs"]["seed"]["outside_hits_by_k"] == [3]
    assert not s.get_meta("selected", [])
    assert expansion_review.prepare_batch(ctrl, "named:T0002")
    packet = ctrl.context(2)
    assert packet["cases"][0]["boundaries"] == receipt["boundaries"]
    reopened = Controller(ctrl.snapshot, ctrl.run, tmp_path/"cache", ctrl.cfg, ADAPTIVE_MODE)
    with closing(reopened.state):
        assert reopened.context(2) == packet
    bad = send(2, dict(action="review_expansion", candidate_id=anchor, choice="foreign", reason="Unbound choice"))
    assert bad["status"] == "needs_clarification" and not s.get_optional("ExpansionReview", anchor)
    assert expansion_review.prepare_batch(ctrl, "named:T0003")
    saved = send(3, dict(action="review_expansion", candidate_id=anchor, choice=option["candidate_id"],
        include_core=True, reason="The added cells share the positive markers despite counter signal"))
    assert saved["choice"] == option["candidate_id"] and saved["include_core"]
    assert ctrl.candidate(saved["working_candidate_id"])["members"] == ctrl.candidate(option["candidate_id"])["members"]
    assert not expansion_review.prepare_batch(ctrl, "named:T0004")
    assert ctrl.read_tool("get_evidence", {"candidate_id":anchor})["pool_followup"]["expansion_review"]["choice"] == saved["choice"]
