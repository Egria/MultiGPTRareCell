"""End-to-end checks of freeze, checkpoint, mode gates, and total assignment."""
import json

import pytest

from rarecell import providers
from rarecell.config import config, fingerprint, read_json
from rarecell.discovery import Controller
from rarecell.evaluation import evaluate
from rarecell.finalize import assignments
from rarecell.prepare import prepare
from rarecell.run import discover, process_message, resume
from rarecell.synthetic import generate


@pytest.fixture
def app(tmp_path):
    dataset = generate(tmp_path / "fixture", 1009)
    cfg = config(override={"mode": "A2", "prepare": {"views": ["V_LINEAR", "V_TENSOR2"], "graph_backend": "igraph"},
                           "provider": {"kind": "scripted_demo"}})
    snap = tmp_path / "snapshot"
    prepare(dataset, snap, cfg)
    return snap, cfg, tmp_path


def start(app, name, **kwargs):
    snap, cfg, root = app
    return discover(snap, root / name, root / "cache", cfg, policy="scripted_demo", **kwargs)


def test_resume_after_real_queries_freeze_and_isolated_evaluation(app):
    _, _, root = app
    run = root / "run"
    assert start(app, "run", pause_after=2)["detail"] == "paused"
    before = read_json(run / "usage.json")
    assert before["compute_executions"] == 2
    initial_freeze = (run/"frozen.json").read_bytes()
    result = resume(run)
    after = read_json(run / "usage.json")
    assert result["status"] == "complete" and result["submitted"]
    assert before["compute_executions"] == after["compute_executions"]
    old = next(p for p in (run/"history").iterdir() if p.is_dir())
    assert (old/"frozen.json").read_bytes() == initial_freeze
    assert evaluate(old, root/"fixture/evaluation_fixture/gt.tsv", root/"eval-paused")["run_summary"]["detail"] == "paused"
    assert after["resume_count"] == 1 and after["deterministic_format_repairs"] == 1
    stable = (run / "frozen.json").read_bytes()
    assert resume(run) == result
    assert stable == (run / "frozen.json").read_bytes()
    metrics = evaluate(run, root/"fixture/evaluation_fixture/gt.tsv", root/"eval")
    assert metrics["run_summary"]["n_cells"] == 160
    assert "planted" not in (run/"manifest.json").read_text()
    with pytest.raises(ValueError):
        evaluate(run, root/"fixture/evaluation_fixture/gt.tsv", run/"eval")
    with (run/"partition.tsv").open("a") as out:
        out.write("altered\n")
    with pytest.raises(ValueError, match="[Cc]hecksum|digest"):
        resume(run)
    with pytest.raises(ValueError, match="digest"):
        evaluate(run, root/"fixture/evaluation_fixture/gt.tsv", root/"eval2")


def test_warm_cache_does_not_buy_more_queries_or_inherit_revisions(app):
    _, _, root = app
    a = start(app, "first")
    b = start(app, "second")
    u1, u2 = [read_json(root / n / "usage.json") for n in ("first", "second")]
    assert u1["new_queries"] == u2["new_queries"] == u2["global_cache_hits"] == 2
    assert u2["compute_executions"] == 0
    assert a["snapshot_hash"] == b["snapshot_hash"] and a["run_id"] != b["run_id"]
    assert (root/"first/partition.tsv").read_bytes() == (root/"second/partition.tsv").read_bytes()
    for mode in ["A0", "A1", "A3"]:
        res = start(app, mode, mode=mode, annotations=root/"no-annotation")
        counts = read_json(root/mode/"usage.json")
        assert counts["compute_executions"] == 0, mode
        if mode == "A0":
            assert res["status"] == "complete" and counts["new_queries"] == 0
        else:
            assert res["detail"] == "annotation_blocked" and not res["submitted"] and counts["turns"] == 0
        assert res["n_cells"] == 160


def test_total_conflicts_and_seed_can_be_reassigned(app):
    snap, cfg, root = app
    ctrl = Controller(snap, root/"logic", root/"cache", cfg, "A2")
    try:
        ctrl.initialize()
        parent = max(ctrl.candidates.values(), key=lambda c: len(c["members"]))
        first, second, third = parent["members"][:3]
        eid = parent["evidence_ids"][0]
        def edit(mid, groups, parents):
            return ctrl.execute({"action":"propose_revision", "args":{"parents":parents, "groups":groups,
                                "reason":"numeric membership inspection", "evidence_ids":[eid]}}, mid)
        with pytest.raises(ValueError, match="overlap"):
            edit("invalid", [{"members":[first, second]}, {"members":[second, third]}], [parent["id"]])
        assert not ctrl.state.all("Decision") and all(c["initial"] for c in ctrl.candidates.values())
        seed = edit("seed", [{"members":[first]}], [parent["id"]])["candidate_ids"][0]
        revised = edit("release", [{"members":[second, third]}], [seed])["candidate_ids"][0]
        ctrl.execute({"action":"submit", "args":{"candidate_ids":[revised]}}, "finish")
        part, _, groups, _ = assignments(ctrl)
        by_cell = {r["cell_id"]:r for r in part}
        assert by_cell[first]["status"] == "unresolved"
        assert by_cell[second]["candidate_id"] == revised
        assert len(by_cell) == len(part) == len(ctrl.cells)
        assert {g["cell_set_hash"] for g in groups} == {
            fingerprint(sorted(r["cell_id"] for r in part if r["group_id"] == gid)) for gid in {r["group_id"] for r in part}}
        ctrl.state.set_meta("selected", [parent["id"], revised])
        part, conflicts, _, _ = assignments(ctrl)
        assert {r["cell_id"] for r in conflicts} == {second, third}
        assert all(r["status"] == "unresolved" for r in part if r["cell_id"] in {second, third})
    finally:
        ctrl.state.close()


def test_mixed_batch_invalid_item_and_duplicate_submit(app):
    snap, cfg, root = app
    ctrl = Controller(snap, root/"batch", root/"cache", cfg, "A0")
    try:
        ctrl.initialize()
        cid = ctrl.ordered()[0]["id"]
        payload = [{"action":"file_read","path":"../../evaluation"}, {"action":"list_candidates","limit":1},
                   {"action":"submit","candidate_ids":[cid]}, {"action":"submit","candidate_ids":[]}]
        process_message(ctrl, {"message_id":"M1","raw":json.dumps(payload)}, 1)
        assert ctrl.state.get_meta("selected") == [cid]
        assert len(ctrl.state.all("Decision")) == 1
        assert len(ctrl.state.all("Receipt")) == 4
        assert ctrl.state.get_meta("protocol_errors") == 1
    finally:
        ctrl.state.close()


def test_budget_or_stall_terminates(app, monkeypatch):
    snap, cfg, root = app
    cfg["budget"].update(stall_turns=2, max_agent_turns=9)
    # Explicit scripted wire input stays inside the test; actual registered tools execute.
    monkeypatch.setattr("rarecell.run.offline_message", lambda *args: '{"action":"list_candidates","limit":1}')
    res = start(app, "stall")
    assert res["detail"] == "stalled" and read_json(root/"stall/usage.json")["turns"] == 3
    cfg["budget"]["max_agent_turns"] = 1
    assert start(app, "budget")["detail"] == "max_agent_turns"
    cfg["budget"]["context_bytes"] = 100
    res = start(app, "too-small")
    assert res["detail"] == "context_blocked" and res["n_cells"] == 160


def test_active_evidence_pages_and_stall_warning(app):
    snap, cfg, root = app
    ctrl = Controller(snap, root / "pages", root / "cache", cfg, "A3")
    try:
        ctrl.initialize()
        candidate = ctrl.ordered()[0]
        def step(turn, **cmd):
            mid = f"step-{turn}"
            process_message(ctrl, {"message_id": mid, "raw": json.dumps(cmd)}, turn)
            return ctrl.state.get("Receipt", mid + "/0")

        comparison = step(1, action="compare_groups", candidate_id=candidate["id"])["result"]
        read = dict(action="get_evidence", evidence_ids=[comparison["evidence_id"]], limit=1)
        assert step(2, **read)["progress"]
        assert step(3, offset=1, **read)["progress"]
        assert not step(4, offset=1, **read)["progress"]
        assert not step(5, offset=0, **read)["progress"]
        before_stop = ctrl.context(6)
        assert before_stop["remaining"]["consecutive_no_progress_turns_before_stop"] == 1
        assert "One more turn without progress" in before_stop["turn_instruction"]
        seed, = ctrl.state.all("Seed")
        assert before_stop["pending_seeds"] == {"count": 1, "candidate_ids": [seed["candidate_id"]]}
        assert step(6, offset=2, **read)["progress"]
        assert ctrl.context(7)["remaining"]["consecutive_no_progress_turns_before_stop"] == cfg["budget"]["stall_turns"]
        assert ctrl.state.usage()["new_queries"] == 1
        assert not ctrl.state.get_meta("submitted", False)
    finally:
        ctrl.state.close()


def test_cited_expansion_validates_then_reuses_the_same_query(app):
    snap, cfg, root = app
    ctrl = Controller(snap, root / "citations", root / "cache", cfg, "A3")
    try:
        ctrl.initialize()
        seed, other = ctrl.ordered()[:2]
        expansion = {"action": "expand_candidate", "candidate_id": seed["id"]}
        def expand(turn, **extra):
            mid = f"expand-{turn}"
            process_message(ctrl, {"message_id": mid, "raw": json.dumps(dict(expansion, **extra))}, turn)
            return ctrl.state.get("Receipt", mid + "/0")
        invalid = expand(1, evidence_ids=["missing"])
        assert invalid["result"]["status"] == "blocked" and "Unknown evidence ID" in invalid["result"]["error"]
        assert ctrl.state.usage()["new_queries"] == 0
        plain = expand(2)
        eid = other["evidence_ids"][0]
        cited = expand(3, evidence_ids=[eid])
        assert not cited["progress"]
        assert cited["result"]["evidence_id"] == plain["result"]["evidence_id"]
        use = ctrl.state.get("QueryUse", cited["result"]["query_use_id"])
        assert use["cited_evidence_ids"] == [eid]
        assert ctrl.state.usage()["new_queries"] == ctrl.state.usage()["compute_executions"] == 1
    finally:
        ctrl.state.close()


def test_ambiguous_api_is_not_automatically_replayed(app, monkeypatch):
    snap, cfg, root = app
    cfg["provider"].update(kind="openai_compatible", model="test", api_key_env=None)
    attempts = []
    def transport(context, options):
        attempts.append(1)
        raise providers.ProviderError("test transport lost after send", "unknown_execution")
    monkeypatch.setattr(providers, "complete", transport)
    run = root/"unknown"
    assert discover(snap, run, root/"cache", cfg)["detail"] == "unknown_execution"
    assert resume(run)["detail"] == "unknown_execution" and len(attempts) == 1
    assert resume(run, retry_unknown=True)["detail"] == "unknown_execution" and len(attempts) == 2
    assert resume(run, retry_unknown=True)["detail"] == "unknown_execution" and len(attempts) == 2


def test_interrupted_export_recovers_from_state(app, monkeypatch):
    _, _, root = app
    import rarecell.run as loop
    original = loop.export
    def interrupted(ctrl, *args):
        (ctrl.run/"partition.tsv").write_text("incomplete")
        raise OSError("injected interrupted export")
    monkeypatch.setattr(loop, "export", interrupted)
    with pytest.raises(OSError):
        start(app, "repair")
    monkeypatch.setattr(loop, "export", original)
    res = resume(root/"repair")
    assert res["submitted"] and res["n_cells"] == 160
    assert read_json(root/"repair/usage.json")["compute_executions"] == 2


def test_contract_keys_and_secret_config():
    from rarecell.config import contract
    assert len(contract()["requirements"]) == 17
    for values in [{"provider":{"api_key":"credential-fixture"}}, {"annotation":{"base_url":"https://user:pass@host"}},
                   {"provider":{"options":{"Authorization":"not-a-real-token"}}}]:
        with pytest.raises(ValueError):
            config(override=values)


def test_truncated_reply_spends_a_turn_and_can_recover(app, monkeypatch):
    snap, cfg, root = app
    cfg["provider"].update(kind="openai_compatible", model="test", api_key_env=None)
    previous = []
    def call(ctx, _):
        previous.append(ctx)
        if len(previous) == 1:
            raise providers.ProviderError("test token cap", "incomplete")
        assert ctx["recent_receipts"][0]["status"] == "needs_clarification"
        return '{"action":"submit", "candidate_ids":[]}', {"api_calls":1}
    monkeypatch.setattr(providers, "complete", call)
    res = discover(snap, root/"short", root/"cache", cfg)
    assert res["submitted"]
    usage = read_json(root/"short/usage.json")
    assert usage["provider_dispatches"] == 2 and usage["new_queries"] == 0 and usage["protocol_errors"] == 1


def test_context_deduplicates_repeated_read_receipts(app):
    snap, cfg, root = app
    ctrl = Controller(snap, root/"dedup", root/"cache", cfg, "A2")
    try:
        ctrl.initialize()
        receipt = ctrl.read_tool("list_candidates", {"limit":2})
        ctxs = []
        for n in (1, 3, 10):
            ctrl.state.set_meta("receipts", [receipt] * n)
            ctxs.append(ctrl.context(2))
        assert ctxs[0] == ctxs[1] == ctxs[2]
        shown = {c["id"] for c in receipt["candidates"]}
        assert not shown.intersection(c["id"] for c in ctxs[0]["current_candidates"])
    finally:
        ctrl.state.close()


@pytest.mark.parametrize("mode", ["A0", "A1", "A2", "A3", "A3_adaptive"])
def test_live_deadline_rejects_recorded_browse_then_accepts_submission(app, mode):
    from rarecell.run import explore
    snap, cfg, root = app
    cfg["budget"].update(max_agent_turns=12, max_actions_per_turn=1)
    ctrl = Controller(snap, root/"deadline", root/"cache", cfg, mode)
    try:
        ctrl.initialize()
        s = ctrl.state
        s.set_meta("run_id", "deadline")
        s.set_meta("finished_turn", 10)
        cid = ctrl.ordered()[0]["id"]
        # Exact last reply from the PGC trace, followed by a mixed correction.
        browse = {"action": "list_candidates", "offset": 0, "limit": 50}
        for turn, payload in [(11, browse), (12, [browse, {"action": "submit", "candidate_ids": [cid]}])]:
            mid = f"deadline:T{turn:04d}"
            s.put("Message", mid, {"message_id": mid, "raw": json.dumps(payload), "usage": {}})
        assert explore(ctrl, "live", pause_after=11) == ("partial", "paused")
        assert "list_candidates" not in ctrl.context(12)["actions"]
        result = s.get("Receipt", "deadline:T0011/0")["result"]
        assert result["status"] == "blocked" and "candidates" not in result
        assert not s.all("SeedFollowupPlan") and not s.all("CanonicalQuery")
        assert explore(ctrl, "live") == ("complete", "submitted")
        assert s.get_meta("selected") == [cid] and s.get_meta("submission_source") == "agent"
        assert s.get_meta("finished_turn") == 12
    finally:
        ctrl.state.close()


@pytest.mark.parametrize("has_draft", [True, False])
def test_ordinary_live_closeout_after_two_nondecisions(app, has_draft):
    from rarecell.run import explore
    snap, cfg, root = app
    cfg["budget"].update(max_agent_turns=3)
    ctrl = Controller(snap, root/"draft", root/"cache", cfg, "A3")
    try:
        ctrl.initialize()
        s = ctrl.state
        s.set_meta("run_id", "draft")
        c = ctrl.ordered()[0]
        first = ({"action": "propose_revision", "parents": [c["id"]], "groups": [{"candidate_ids": [c["id"]]}],
                  "reason": "retain measured membership", "evidence_ids": c["evidence_ids"]} if has_draft
                 else {"action": "get_evidence", "candidate_id": c["id"]})
        for turn, action in enumerate([first, {"action":"expand_candidate", "candidate_id":c["id"]},
                                      {"action":"list_candidates", "offset":0, "limit":50}], 1):
            mid = f"draft:T{turn:04d}"
            s.put("Message", mid, {"message_id": mid, "raw": json.dumps(action), "usage": {}})
        assert explore(ctrl, "live") == ("partial", "max_agent_turns" if has_draft else "decision_missing")
        assert s.get_meta("closeout_trigger") == "max_agent_turns"
        assert not s.all("CanonicalQuery") and not s.all("SeedFollowup")
        assert bool(s.get_meta("submitted")) == has_draft
        if has_draft:
            assert s.get_meta("submission_source") == "controller_budget"
            assert s.get_meta("submission_status") == "provisional"
    finally:
        ctrl.state.close()


def test_seed_screen_and_journalled_followup_on_browse(app):
    from rarecell.run import follow_browsing_with_seed, explore
    snap, cfg, root = app
    run = root/"seed-review"
    ctrl = Controller(snap, run, root/"cache", cfg, "A3")
    try:
        ctrl.initialize()
        s = ctrl.state
        s.set_meta("run_id", "seeds")
        for i in range(1, 3):
            mid = f"seeds:T{i:04d}"
            s.put("Message", mid, {"message_id": mid, "raw": '{"action":"list_candidates"}', "usage": {}})
        assert explore(ctrl, "live", pause_after=1) == ("partial", "paused")
        seed, = ctrl.screened_seeds
        record = s.get("SeedFollowup", "seeds:T0001")
        assert record["source"] == "controller_seed_review" and record["anchor_id"] == seed["candidate_id"]
        assert record["result"]["query_key"] and s.usage()["new_queries"] == 1
        assert not ctrl.context(2)["seed_review"]["priority"]
    finally:
        ctrl.state.close()
    ctrl = Controller(snap, run, root/"cache", cfg, "A3")
    try:
        s = ctrl.state
        assert follow_browsing_with_seed(ctrl, "seeds:T0001") == record
        assert explore(ctrl, "live", pause_after=2) == ("partial", "paused")
        assert len(s.all("SeedFollowupPlan")) == len(s.all("Investigation")) == 1
        assert s.usage()["compute_executions"] == 1
    finally:
        ctrl.state.close()


@pytest.mark.parametrize("setting", ["max_new_queries", "max_seed_followups", "max_expansion_depth", "max_actions_per_turn"])
def test_seed_followup_respects_all_caps(app, setting):
    snap, cfg, root = app
    cfg["budget"][setting] = 1 if setting == "max_actions_per_turn" else 0
    ctrl = Controller(snap, root/"caps", root/"cache", cfg, "A2")
    try:
        ctrl.initialize()
        ctrl.state.set_meta("policy", "live")
        process_message(ctrl, {"message_id": "cap:T0001", "raw": '{"action":"list_candidates"}'}, 1)
        assert not ctrl.state.all("SeedFollowupPlan") and not ctrl.state.all("CanonicalQuery")
    finally:
        ctrl.state.close()


def test_live_ordinary_stall_enters_decision_with_no_draft(app):
    from rarecell.run import explore
    snap, cfg, root = app
    cfg["budget"].update(stall_turns=1, max_seed_followups=0)
    ctrl = Controller(snap, root/"early-stall", root/"cache", cfg, "A3")
    try:
        ctrl.initialize()
        s = ctrl.state
        s.set_meta("run_id", "stall")
        for i, payload in enumerate([{"action":"list_candidates"}, {"action":"list_candidates"},
                                      {"action":"submit", "candidate_ids":[]}], 1):
            mid = f"stall:T{i:04d}"
            s.put("Message", mid, {"message_id":mid, "raw":json.dumps(payload), "usage":{}})
        assert explore(ctrl, "live", pause_after=2) == ("partial", "paused")
        assert set(ctrl.context(3)["actions"]) == {"submit", "propose_revision"}
        assert explore(ctrl, "live") == ("complete", "submitted")
        assert s.get_meta("decision_start_turn") == 3 and s.get_meta("selected") == []
    finally:
        ctrl.state.close()


def test_frozen_seed_and_decision_diagnostics(app, monkeypatch):
    import csv
    from rarecell.dashboard_transcript import discovery_blocks
    snap, cfg, root = app
    cfg["provider"].update(kind="openai_compatible", model="test", api_key_env=None)
    cfg["budget"]["max_agent_turns"] = 3
    calls = []
    def reply(ctx, _):
        calls.append(ctx)
        if len(calls) > 1:
            assert "list_candidates" not in ctx["actions"] and ctx["phase"] == "decision"
        return '{"action":"list_candidates","offset":0,"limit":50}', {"api_calls":1}
    monkeypatch.setattr(providers, "complete", reply)
    run = root/"frozen-review"
    summary = discover(snap, run, root/"cache", cfg)
    assert len(calls) == 3 and calls[0]["seed_review"]["screened"] == 1
    assert summary["detail"] == "decision_missing" and summary["closeout_trigger"] == "max_agent_turns"
    assert not summary["submitted"] and summary["n_expanded_seeds"] == summary["n_controller_seed_followups"] == 1
    assert summary["n_screened_seeds"] == 1 and summary["unresolved_cells"] == summary["n_cells"]
    with (run/"seeds.tsv").open() as stream:
        row, = csv.DictReader(stream, delimiter="\t")
        assert row["screened"] == "True" and row["task_status"] == "investigated" and row["controller_followup"]
    display = "\n".join(discovery_blocks(run, read_json(run/"checkpoint.json")))
    assert "控制器种子补查" in display and "controller_seed_review" in display and "decision_missing" in display
