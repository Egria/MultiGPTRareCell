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
