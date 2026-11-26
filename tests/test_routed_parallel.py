from collections import Counter
from contextlib import closing
import json
from threading import Barrier

import pytest

from rarecell import providers
from rarecell.config import read_json, read_lines
from rarecell.routed import RoutedDiscovery
from rarecell.run import discover, resume
from rarecell.state import State
from test_routed import setup_input


def reviews(ctx):
    return dict(action="route_seeds", reviews=[
        dict(candidate_id=s["candidate_id"], guide_label="L" if s["candidate_id"] in ("a", "c") else "M",
             reliability="supported", reason="Audited seed")
        for s in ctx["seeds"]
    ])


def test_independent_triage_and_names_have_isolated_context_and_ordered_publication(tmp_path, monkeypatch):
    cfg, snapshot, annotations, _ = setup_input(tmp_path)
    cfg["budget"]["seed_review_batch_size"] = 1
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **k: [])
    barriers = {"seed_triage": Barrier(4), "blind_name_search": Barrier(2)}

    def reply(ctx, cfg):
        phase = ctx["phase"]
        if cfg.get("routed_workers") == 4:
            barriers[phase].wait(timeout=5)
        if phase == "seed_triage":
            assert len(ctx["seeds"]) == 1
            action = reviews(ctx)
        else:
            assert not {"seed_review", "scopes", "annotations"}.intersection(ctx)
            gene = {"L": "G", "M": "H"}[ctx["target"]]
            if ctx["history"]:
                previous, = ctx["history"]
                assert set(previous["response"]["genes"]) == {gene}
                action = dict(action="freeze", support=80, reason="measured",
                              rule=[gene, 3] if gene == "G" else None)
            else:
                action = dict(action="survey", genes=[gene])
        return json.dumps(action), {"api_calls": 1, "total_tokens": 10}

    monkeypatch.setattr(providers, "complete", reply)
    for workers in (1, 4):
        run = tmp_path / f"run{workers}"
        result = discover(snapshot, run, tmp_path / "cache", dict(cfg, routed_workers=workers), annotations=annotations)
        assert result["submitted"]
        claims = read_lines(run / "route_claims.jsonl")
        assert [(r["label"], r["n"]) for r in claims] == [("L", 2), ("M", 0)]
        with closing(State(run, tmp_path / "cache")) as st:
            first = st.get_optional("RouteTurnPlan", "1")
            if workers == 1:
                assert first is None
            else:
                assert first["last_turn"] == 4
                assert st.get("RouteTurnPlan", "5")["last_turn"] == 6
                assert st.get("RouteTurnPlan", "7")["last_turn"] == 8
                ordered = [r["message_id"] for r in st.all("Receipt") if r["result"].get("action") == "freeze"]
                assert [x.split(":")[-1] for x in ordered] == ["T0007", "T0008"]


def test_paid_wave_replayed_with_one_worker_after_budget_stop(tmp_path, monkeypatch):
    cfg, snapshot, annotations, _ = setup_input(tmp_path)
    cfg["routed_workers"] = 4
    cfg["budget"].update(max_new_queries=2, max_agent_turns=9)
    cfg["potential_names"] = ["L", "M"]
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **k: [])
    called = Counter()

    def reply(ctx, cfg):
        if ctx["phase"] == "seed_triage":
            action = reviews(ctx)
        elif ctx.get("target") == "L":
            called["L"] += 1
            action = [dict(action="survey", genes=["G"]), dict(action="survey", genes=["H"]),
                      dict(action="views", rules=["G"], genes=["H"])] if not ctx["history"] else dict(
                          action="freeze", support=70, reason="done", rule=None)
        else:
            called["M"] += 1
            action = dict(action="freeze", support=70, reason="done", rule=None)
        return json.dumps(action), {"api_calls": 1, "total_tokens": 8}

    monkeypatch.setattr(providers, "complete", reply)
    run = tmp_path / "run"
    first = discover(snapshot, run, tmp_path / "cache", cfg, annotations=annotations)
    assert first["detail"] == "routed_query_budget" and not first["submitted"]
    assert read_json(run / "usage.json")["new_queries"] == 2
    with closing(State(run, tmp_path / "cache")) as st:
        assert len(st.all("Message")) == 3 and st.get_meta("finished_turn") == 1
    assert called == {"L": 1, "M": 1}
    after = resume(run, max_new_queries=10, routed_workers=1)
    assert after["submitted"] and called == {"L": 2, "M": 1}
    assert read_json(run / "usage.json")["new_queries"] == 3


def test_search_errors_do_not_leak_between_agents(tmp_path, monkeypatch):
    cfg, snapshot, annotations, _ = setup_input(tmp_path)
    cfg["routed_workers"] = 2
    cfg["budget"]["stall_turns"] = 2
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **k: [])

    def reply(ctx, cfg):
        if ctx["phase"] == "seed_triage":
            action = reviews(ctx)
        elif ctx["target"] == "L":
            action = dict(action="wrong")
        else:
            assert "previous_error" not in ctx
            action = dict(action="survey", genes=["H"])
        return json.dumps(action), {"api_calls": 1}

    monkeypatch.setattr(providers, "complete", reply)
    run = tmp_path / "run"
    result = discover(snapshot, run, tmp_path / "cache", cfg, annotations=annotations)
    assert result["detail"] == "routed_action_blocked" and "name/l" in result["issue"]
    with closing(State(run, tmp_path / "cache")) as st:
        assert st.get_meta("action_errors") == 2
        plan = st.get("RouteTurnPlan", "4")
        assert "previous_error" in st.read_artifact(plan["context_ref"])


@pytest.mark.parametrize("pause,token_cap", [(1, 0), (None, 1)])
def test_global_budget_gates_dispatch(tmp_path, monkeypatch, pause, token_cap):
    cfg, snapshot, annotations, _ = setup_input(tmp_path)
    cfg.update(routed_workers=4)
    cfg["budget"].update(seed_review_batch_size=1, max_reported_tokens=token_cap)
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **k: [])
    monkeypatch.setattr(providers, "complete", lambda ctx, cfg: (json.dumps(reviews(ctx)), {"total_tokens": 5}))
    run = tmp_path / "run"
    result = discover(snapshot, run, tmp_path / "cache", cfg, annotations=annotations, pause_after=pause)
    assert result["detail"] == ("paused" if pause else "max_reported_tokens")
    assert read_json(run / "usage.json")["provider_dispatches"] == 1


def test_interrupted_serial_turn_replays_before_parallel_work(tmp_path, monkeypatch):
    cfg, snapshot, annotations, _ = setup_input(tmp_path)
    cfg["budget"]["seed_review_batch_size"] = 1
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **k: [])
    called = Counter()

    def reply(ctx, cfg):
        if ctx["phase"] == "seed_triage":
            called.update(s["candidate_id"] for s in ctx["seeds"])
            action = reviews(ctx)
        else:
            action = dict(action="freeze", support=70, reason="done", rule=None)
        return json.dumps(action), {"api_calls": 1}

    monkeypatch.setattr(providers, "complete", reply)
    run = tmp_path / "run"
    with monkeypatch.context() as crash:
        def stop(*args):
            raise SystemExit()
        crash.setattr(RoutedDiscovery, "apply", stop)
        with pytest.raises(SystemExit):
            discover(snapshot, run, tmp_path / "cache", cfg, annotations=annotations)
    assert called == {"a": 1}
    result = resume(run, routed_workers=4)
    assert result["submitted"]
    assert called == {"a": 1, "c": 1, "d": 1, "f": 1}
    with closing(State(run, tmp_path / "cache")) as st:
        assert st.get_optional("RouteTurnPlan", "1") is None
        assert st.get("RouteTurnPlan", "2")["last_turn"] == 4


def test_parallel_seed_searches_nominate_one_blind_name_with_all_four_parents(tmp_path, monkeypatch):
    cfg, snapshot, annotations, _ = setup_input(tmp_path)
    cfg["routed_workers"] = 4
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **k: [])
    barrier = Barrier(4)
    seeds, names = [], []

    def reply(ctx, cfg):
        if ctx["phase"] == "seed_triage":
            req = dict(action="route_seeds", reviews=[dict(candidate_id=row["candidate_id"], guide_label="",
                       reliability="unknown", reason="Uncertain identity") for row in ctx["seeds"]])
        elif ctx["phase"] == "seed_expansion":
            seeds.append(ctx["seed_review"]["candidate_id"])
            assert ctx["scopes"]["seed"] == 1 and "retrieval" in ctx
            barrier.wait(timeout=5)
            req = dict(action="freeze", kind="state", rule={"scope": "seed"}, proposed_label="Novel",
                       support=88, reason="Suggest independent search by name")
        else:
            names.append(ctx["target"])
            assert not {"seed_review", "seed_markers", "scopes", "retrieval"}.intersection(ctx)
            req = dict(action="freeze", rule="G", support=78, reason="Measured G expression")
        return json.dumps(req), {"api_calls": 1}

    monkeypatch.setattr(providers, "complete", reply)
    run = tmp_path / "run"
    result = discover(snapshot, run, tmp_path / "cache", cfg, annotations=annotations)
    assert result["submitted"] and names == ["Novel"] and sorted(seeds) == ["a", "c", "d", "f"]
    with closing(State(run, tmp_path / "cache")) as state:
        freezes = state.all("RouteFreeze")
        seeds = {r["candidate_id"] for r in freezes if r["job"].startswith("seed/")}
        named, = [r for r in freezes if r["job"] == "name/novel"]
        assert set(state.get("Candidate", named["candidate_id"])["parents"]) == seeds
