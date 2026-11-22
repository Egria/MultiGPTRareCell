"""Exercise the mandatory phase with real matrices and recorded provider-interface replies."""
import csv
import json

import pytest

from rarecell import exhaustive, job_queue, providers
from rarecell.cli import parser
from rarecell.config import config, read_json
from rarecell.discovery import Controller
from rarecell.prepare import prepare
from rarecell.protocol import parse_message
from rarecell.run import discover, explore, resume
from rarecell.synthetic import generate


@pytest.fixture(scope="module")
def frozen(tmp_path_factory):
    root = tmp_path_factory.mktemp("seed-input")
    dataset = generate(root/"input", 1009)
    cfg = config(override={"prepare": {"views": ["V_LINEAR"], "graph_backend": "igraph", "outlier_mad": 1.5}})
    prepare(dataset, root/"snapshot", cfg)
    return root/"snapshot"


def strict(**budget):
    return config(override={"mode": "A2", "exhaustive_seed_review": True,
                            "provider": {"model": "test-model", "api_key_env": None},
                            "budget": {"max_agent_turns": 3, "max_new_queries": 1, "seed_review_batch_size": 4, **budget}})


def written(row):
    return {"seed_id": row["seed_id"], "recommendation": "follow up", "explanation":
            f"Measured {row['measurement']['status']} for {row['candidate_id']}. Compare its peer before selection."}


def test_all_seeds_measured_and_written_after_partial_batch_and_resume(frozen, tmp_path, monkeypatch):
    calls = []
    def reply(ctx, _cfg):
        calls.append(ctx)
        if ctx.get("phase") == "exhaustive_seed_review":
            notes = [written(seed) for seed in ctx["seeds"]]
            if len(calls) == 1:
                notes = notes[:1] + [{"action":"submit", "candidate_ids":[]}]
            response = {"reviews": notes}
        else:
            assert ctx["remaining"]["turns_including_this"] == 3
            assert ctx["remaining"]["new_queries"] == 1 and ctx["written_seed_reviews"]
            response = {"action":"submit", "candidate_ids":[]}
        return json.dumps(response), {"api_calls": 1, "total_tokens": 10}
    monkeypatch.setattr(providers, "complete", reply)
    run = tmp_path/"run"
    paused = discover(frozen, run, tmp_path/"cache", strict(), pause_after=1)
    assert paused["detail"] == "paused" and not paused["submitted"]
    initial = paused["seed_review_coverage"]
    assert initial["total"] > 4 and initial["measured"] == 4 and initial["written"] == 1
    first_batch = [s["candidate_id"] for s in calls[0]["seeds"]]
    result = resume(run)
    counts = result["seed_review_coverage"]
    assert counts["total"] == counts["written"] == counts["completed"] == counts["measured"]
    assert result["submitted"] and counts["remaining"] == 0
    assert first_batch[0] not in [s["candidate_id"] for s in calls[1]["seeds"]]
    assert len(calls) == counts["review_turns"] + 1
    usage = read_json(run/"usage.json")
    assert usage["new_queries"] == usage["compute_executions"] == counts["total"]
    assert usage["action_errors"] == 1 and usage["provider_dispatches"] == len(calls)
    with (run/"seed_reviews.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len({r["candidate_id"] for r in rows}) == len(rows) == counts["total"]
    assert all(r["query_key"] and r["evidence_id"] and r["reason"] and r["message_id"] for r in rows)
    assert "seed_reviews.tsv" in read_json(run/"frozen.json")["files"]
    assert resume(run) == result


@pytest.mark.parametrize("ending", ["stall", "token"])
def test_missing_written_coverage_never_submits(frozen, tmp_path, monkeypatch, ending):
    calls = []
    def reply(ctx, _):
        calls.append(ctx)
        seed = ctx["seeds"][0]
        return json.dumps([{"action":"submit", "candidate_ids":[]}, dict(written(seed), evidence_ids=["wrong"])]), {"api_calls": 1, "total_tokens": 10}
    monkeypatch.setattr(providers, "complete", reply)
    result = discover(frozen, tmp_path/"run", tmp_path/"cache",
                      strict(stall_turns=2, max_reported_tokens=1 if ending=="token" else 0))
    assert result["detail"] == ("max_reported_tokens" if ending=="token" else "seed_review_incomplete")
    cov = result["seed_review_coverage"]
    assert not result["submitted"] and cov["written"] == 0 and cov["remaining"] == cov["total"]
    if ending == "token":
        assert cov["measured"] == 4 and len(calls) == 1
    else:
        assert cov["measured"] == cov["total"] and len(calls) == (cov["total"] + 3)//4 + 2
        order = [seed["candidate_id"] for ctx in calls for seed in ctx["seeds"]]
        assert len(set(order[:cov["total"]])) == cov["total"]


def test_seed_shorthand_only_in_assigned_phase():
    row = {"seed_id": "S1", "recommendation": "keep", "explanation": "Coherent positive signal"}
    raw = json.dumps({"assessments": [row]})
    assert parse_message(raw, seed_batch=True)[0]["args"] == {
        "candidate_id": "S1", "disposition": "retain", "reason": row["explanation"]}
    assert "error" in parse_message(raw)[0]
    assert "error" in parse_message(json.dumps(dict(row, explanation=" ")), seed_batch=True)[0]


def test_cli_queue_and_mode_boundary(frozen, tmp_path):
    args = ["discover", "--run", "r", "--snapshot", "s"]
    assert parser().parse_args(args + ["--exhaustive-seed-review"]).exhaustive_seed_review is True
    assert parser().parse_args(args + ["--no-exhaustive-seed-review"]).exhaustive_seed_review is False
    assert not config()["exhaustive_seed_review"]
    queue = tmp_path/"queue"
    job_queue.initialize(queue)
    form = {"kind":"discover", "snapshot":str(frozen), "mode":"A2", "model":"test", "exhaustive_seed_review":"yes"}
    task = job_queue.submit(queue, tmp_path, form)
    assert config(queue/task/"settings.yaml")["exhaustive_seed_review"]
    with pytest.raises(ValueError):
        job_queue.submit(queue, tmp_path, dict(form, mode="A0"))
    with pytest.raises(ValueError, match="requires live"):
        discover(frozen, tmp_path/"invalid", tmp_path/"cache", strict(), policy="fixed", mode="A3")


def test_crash_after_first_note_in_paid_batch(frozen, tmp_path, monkeypatch):
    requested = []
    def reply(ctx, _):
        requested.append(ctx["seeds"])
        return json.dumps([written(row) for row in ctx["seeds"]]), {"api_calls": 1}
    monkeypatch.setattr(providers, "complete", reply)
    def open_controller():
        return Controller(frozen, tmp_path/"run", tmp_path/"cache", strict(), "A2")
    ctrl = open_controller()
    try:
        s = ctrl.state
        s.set_meta("run_id", "interrupted")
        ctrl.initialize()
        original = s.event
        def stop(kind, event):
            if kind == "action" and event.get("index") == 1:
                raise RuntimeError("simulated process interruption")
            return original(kind, event)
        monkeypatch.setattr(s, "event", stop)
        with pytest.raises(RuntimeError, match="interruption"):
            explore(ctrl, "live", pause_after=1)
        assert len(s.all("SeedAssessment")) == 1 and not s.get_meta("finished_turn")
    finally:
        ctrl.state.close()
    ctrl = open_controller()
    try:
        assert explore(ctrl, "live", pause_after=1) == ("partial", "paused")
        assert len(requested) == 1 and len(ctrl.state.all("SeedAssessment")) == 4
        assert ctrl.state.usage()["compute_executions"] == 4 and ctrl.review_turns == 1
    finally:
        ctrl.state.close()


def test_measurement_failure_audited_and_other_seeds_still_reviewed(frozen, tmp_path, monkeypatch):
    from rarecell.science import Science
    original = Science.execute
    failed = []
    def fail_one(self, name, args):
        if not failed:
            failed.append(tuple(args["members"]))
        if tuple(args["members"]) == failed[0]:
            raise ValueError("test measurement failure")
        return original(self, name, args)
    monkeypatch.setattr(Science, "execute", fail_one)
    def reply(ctx, _):
        return json.dumps([written(row) for row in ctx["seeds"]]), {"api_calls":1}
    monkeypatch.setattr(providers, "complete", reply)
    result = discover(frozen, tmp_path/"run", tmp_path/"cache", strict())
    assert result["detail"] == "seed_measurement_blocked" and not result["submitted"]
    cov = result["seed_review_coverage"]
    assert cov["total"] == cov["written"] and cov["measurement_blocked"] == cov["remaining"] == 1
    with (tmp_path/"run/seed_reviews.tsv").open() as f:
        row, = [r for r in csv.DictReader(f, delimiter="\t") if r["measurement_status"] == "blocked"]
        assert "test measurement failure" in row["measurement_error"] and row["reason"]


def test_inapplicable_seed_batch_and_explicit_empty_set(frozen, tmp_path):
    ctrl = Controller(frozen, tmp_path/"run", tmp_path/"cache", strict(), "A2")
    try:
        ctrl.initialize()
        ctrl.state.set_meta("run_id", "empty")
        with pytest.raises(ValueError, match="every seed"):
            ctrl.decide("submit", {"candidate_ids":[]}, "early")
        with pytest.raises(ValueError, match="assigned"):
            exhaustive.record(ctrl, "wrong", {"candidate_id":"unassigned"})
        # An initialized controller can legitimately have an empty frozen seed inventory.
        ctrl.state.db.execute("DELETE FROM record WHERE kind='Seed'")
    finally:
        ctrl.state.close()
    ctrl = Controller(frozen, tmp_path/"run", tmp_path/"cache", strict(), "A2")
    try:
        mid = "empty:T0001"
        ctrl.state.put("Message", mid, {"message_id":mid, "raw":'{"action":"submit","candidate_ids":[]}', "usage":{}})
        assert explore(ctrl, "live") == ("complete", "submitted")
        assert ctrl.review_turns == 0 and not ctrl.state.all("SeedBatch")
    finally:
        ctrl.state.close()


@pytest.mark.parametrize("mode", ["A3", "A3_adaptive"])
def test_annotated_controller_requires_seed_review_before_closeout(frozen, tmp_path, monkeypatch, mode):
    ctrl = Controller(frozen, tmp_path/"run", tmp_path/"cache", strict(max_agent_turns=1), mode)
    calls = []
    def reply(ctx, _):
        calls.append(ctx)
        if ctx["phase"] == "exhaustive_seed_review":
            return json.dumps([written(r) for r in ctx["seeds"]]), {"api_calls":1}
        assert ctx["phase"] == "decision" and "submit" in ctx["actions"]
        return '{"action":"submit","candidate_ids":[]}', {"api_calls":1}
    monkeypatch.setattr(providers, "complete", reply)
    try:
        ctrl.initialize()
        ctrl.state.set_meta("run_id", "annotated")
        assert explore(ctrl, "live") == ("complete", "submitted")
        assert exhaustive.coverage(ctrl)["remaining"] == 0 and len(calls) == ctrl.review_turns + 1
    finally:
        ctrl.state.close()


def test_parent_annotation_join_binds_both_group_and_resolution(tmp_path):
    from scripts.audit_seed_annotations import memberships
    (tmp_path/"seeds.jsonl").write_text('{"cell_id":"x","candidate_id":"seed-x"}\n')
    (tmp_path/"baselines.tsv").write_text("view\tresolution\tcell_id\tgroup_id\n"
        "V_HVG\t0.4\tx\tH1\nV_HVG\t1\tx\tH1\nV_LINEAR\t0.4\tx\tH2\nV_HVG\t1\ty\tH2\n")
    ann = [{"resolution": .4, "candidate_id": "H1", "label":"unknown", "selection_history":{"genes":["m","n"]}},
           {"resolution": 1, "candidate_id": "H1", "label":"b", "selection_history":{"genes":["p"]}}]
    rows = list(memberships(tmp_path, ann, {"x", "y"}))
    assert [(r["cell_id"], r["parent_label"], r["seed_candidate_id"], r["marker_input"]) for r in rows] == [
        ("x","unknown","seed-x","m,n"), ("x","b","seed-x","p"), ("y","","","")]
