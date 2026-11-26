from collections import Counter
from contextlib import closing
import json
import threading

import pytest

from rarecell import forced_routing, providers
from rarecell.config import ROUTED_MODE
from rarecell.discovery import Controller
from rarecell.routed import RoutedDiscovery
from rarecell.run import obtain
from test_routed import setup_input


def test_forced_review_parallel_siblings_survive_invalid_reply_and_interruption(tmp_path, monkeypatch):
    cfg, snapshot, _, _ = setup_input(tmp_path)
    cfg.update(routed_workers=2, routed_conflict_policy="single")
    ctrl = Controller(snapshot, tmp_path / "run", tmp_path / "cache", cfg, ROUTED_MODE)
    route = RoutedDiscovery(ctrl)
    old_batches = forced_routing._batches
    monkeypatch.setattr(forced_routing, "_batches", lambda rows, prompt, limit, size: old_batches(rows, prompt, limit, 1))
    barrier = threading.Barrier(2)
    calls = Counter()

    def reply(ctx, config):
        phase = ctx["phase"]
        assert phase in {"conflict_label_audit", "conflict_panel_nomination"}
        if phase == "conflict_label_audit":
            key, = ctx["candidates"]
            calls[key] += 1
            if "previous_error" not in ctx:
                barrier.wait(timeout=5)
            data = {"reviews": {} if key == "a" and calls[key] == 1 else {key: "supported"}}
        else:
            # Nominations must see all decisions, including the repaired first audit.
            assert calls["a"] == 2 and calls["c"] == 1
            assert ctx["vocabulary"] == ["cell A", "cell C"]
            name, = ctx["labels"]
            calls[name] += 1
            data = {"labels": {name: dict(canonical=name, panels=[["G"]])}}
        return json.dumps(data), {"api_calls": 1, "total_tokens": 3}

    monkeypatch.setattr(providers, "complete", reply)
    with closing(ctrl.state) as s:
        s.set_meta("run_id", "test")
        ctrl.initialize()
        for cid in ("a", "c"):
            s.put("RouteFreeze", cid, dict(kind="lineage", job="name/" + cid, candidate_id=cid,
                                           label="cell " + cid.upper(), reason="example", support=85))
        original = route.complete_turn
        def interrupt(mid, receipt):
            if mid.endswith("T0002"):
                raise KeyboardInterrupt()
            original(mid, receipt)
        with monkeypatch.context() as crash:
            crash.setattr(route, "complete_turn", interrupt)
            with pytest.raises(KeyboardInterrupt):
                forced_routing.finish(route, obtain, None, False)
        assert s.get_meta("finished_turn") == 1
        assert len(s.all("Message")) == 2 and not s.all("ForcedReview")
    cfg["routed_workers"] = 1
    ctrl = Controller(snapshot, tmp_path / "run", tmp_path / "cache", cfg, ROUTED_MODE)
    with closing(ctrl.state) as s:
        assert forced_routing.finish(RoutedDiscovery(ctrl), obtain, None, False) == ("complete", "submitted")
        assert s.get_meta("finished_turn") == 5 and not s.get_meta("forced_review_wave")
        assert len(s.all("ForcedReview")) == 4 and len(s.all("ProviderDispatch")) == 5
        assert len(s.all("ForcedReviewError")) == 1 and calls["c"] == 1


def test_token_cap_serializes_audits_and_allows_paid_reply_after_crash(tmp_path, monkeypatch):
    cfg, snapshot, _, _ = setup_input(tmp_path)
    cfg.update(routed_workers=3, routed_conflict_policy="single")
    cfg["budget"]["max_reported_tokens"] = 1
    ctrl = Controller(snapshot, tmp_path / "run", tmp_path / "cache", cfg, ROUTED_MODE)
    route = RoutedDiscovery(ctrl)
    batching = forced_routing._batches
    monkeypatch.setattr(forced_routing, "_batches", lambda rows, prompt, limit, size: batching(rows, prompt, limit, 1))
    def reply(ctx, _):
        assert list(ctx["candidates"]) == ["a"]
        return json.dumps({"reviews": {"a": "supported"}}), {"api_calls": 1, "total_tokens": 20}
    monkeypatch.setattr(providers, "complete", reply)
    with closing(ctrl.state) as s:
        s.set_meta("run_id", "capped")
        ctrl.initialize()
        for cid in ("a", "c"):
            s.put("RouteFreeze", cid, dict(kind="lineage", job="name/" + cid, candidate_id=cid,
                                           label="cell " + cid, reason="example", support=75))
        def crash(*args):
            raise KeyboardInterrupt
        with monkeypatch.context() as interrupted:
            interrupted.setattr(route, "complete_turn", crash)
            with pytest.raises(KeyboardInterrupt):
                forced_routing.finish(route, obtain, None, False)
        assert len(s.all("ProviderDispatch")) == 1 and not s.all("ForcedReview")
        assert forced_routing.finish(route, obtain, None, False) == ("partial", "max_reported_tokens")
        assert len(s.all("ProviderDispatch")) == len(s.all("ForcedReview")) == s.get_meta("finished_turn") == 1


@pytest.mark.parametrize("workers", [1, 2])
def test_resume_terminal_review_from_previous_backend(tmp_path, monkeypatch, workers):
    from rarecell import run as runner
    from rarecell.config import RESOURCES, fingerprint, hash_file, read_json
    from rarecell.forced_replay import replay

    cfg, snap, ann, _ = setup_input(tmp_path)
    cfg.update(routed_workers=workers, routed_conflict_policy="single", potential_names=list("ABCD"))
    cfg["budget"]["max_agent_turns"] = 80
    monkeypatch.setattr(runner, "verify_annotation", lambda *a, **kw: [])
    batches = forced_routing._batches
    monkeypatch.setattr(forced_routing, "_batches", lambda rows, prompt, limit, size: batches(rows, prompt, limit, min(size, 2)))
    active = []
    def reply(ctx, config):
        phase = ctx["phase"]
        active.append((phase, ctx.get("labels")))
        if phase == "seed_triage":
            result = dict(action="route_seeds", reviews=[
                dict(candidate_id=row["candidate_id"], reliability="uncertain", reason="uncertain", guide_label="")
                for row in ctx["seeds"]])
        elif phase == "blind_name_search":
            result = dict(action="freeze", rule={"A": "G", "B": "H", "C": "U", "D": ["G", 3]}[ctx["target"]],
                          reason="Measured", support=86, label_reliability="supported")
        elif phase == "seed_expansion":
            result = dict(action="freeze", rule=None, reason="abstain", support=0)
        elif phase == "conflict_panel_nomination":
            labels = ctx["labels"]
            result = {"labels": {k: dict(canonical=k, panels=[["G"]]) for k in labels if k != "B"}}
            if labels == ["C", "D"] and "previous_error" not in ctx:
                result = {"labels": {}}
        else:
            assert phase == "conflict_measured_review"
            result = {"rules": {k: dict(default="geometry") for k in ctx["pairs"]}}
        return json.dumps(result), {"api_calls": 1, "total_tokens": 10}

    monkeypatch.setattr(providers, "complete", reply)
    def legacy_finish(route, obtain, pause_after, retry_unknown):
        # Write the old on-disk layout, including the consumed wave without an accepted review.
        s = route.s
        claims = [dict(c, members=route.ctrl.candidate(c["candidate_id"])["members"])
                  for c in s.all("RouteFreeze") if c.get("kind") == "lineage" and c.get("candidate_id")]
        def old_wave(requests):
            pending = requests[:workers]
            wave = []
            for attempt in (1, 2):
                for ctx, validate in pending:
                    assert ctx["phase"] == "conflict_panel_nomination"
                    key = fingerprint(ctx)
                    prev = s.get_optional("ForcedReviewError", key)
                    sent = dict(ctx, previous_error=prev["error"]) if prev else ctx
                    turn = s.get_meta("finished_turn") + 1
                    msg = obtain(route.ctrl, "live", turn, False, context=sent)
                    mid = msg["message_id"]
                    if attempt == 2:
                        wave.append(dict(key=key, turn=turn, previous_error=prev["error"]))
                    parsed = json.loads(msg["raw"])
                    try:
                        validate(parsed)
                    except ValueError as exc:
                        if not prev:
                            s.put("ForcedReviewError", key, {"error": str(exc)})
                        route.complete_turn(mid, {"error": str(exc)})
                    else:
                        s.put("ForcedReview", key, dict(context_ref=s.artifact(ctx), reply=parsed, message_id=mid))
                        route.complete_turn(mid, dict(action=ctx["phase"]))
            if workers > 1:
                s.set_meta("forced_review_wave", wave)
            raise forced_routing.StopReview("routed_action_blocked")
        with pytest.raises(forced_routing.StopReview):
            forced_routing.adjudicate(route.ctrl, claims, None, ask_many=old_wave)
        s.set_meta("current_issue", "routed_action_blocked")
        return "partial", "routed_action_blocked"

    run = tmp_path / "run"
    with monkeypatch.context() as legacy:
        legacy.setattr(forced_routing, "finish", legacy_finish)
        legacy.setattr(runner, "code_fingerprint", lambda: read_json(RESOURCES / "routed_checkpoint_upgrade.json")["also_from"][-1])
        old = runner.discover(snap, run, tmp_path / "cache", cfg, annotations=ann)
    assert not old["submitted"] and not old["resumable"] and old["issue"] == "routed_action_blocked"
    digest = hash_file(run / "frozen.json")
    before_claims = (run / "route_claims.jsonl").read_bytes()
    active.clear()
    result = runner.resume(run)
    assert result["submitted"] and result["submission_source"] == "forced_single_choice"
    assert old["run_id"] == result["run_id"] and (run / "route_claims.jsonl").read_bytes() == before_claims
    assert (run / "history" / digest / "frozen.json").is_file()
    assert set(p for p, _ in active) <= {"conflict_panel_nomination", "conflict_measured_review"}
    panels = [labels for phase, labels in active if phase == "conflict_panel_nomination"]
    assert ["A", "B"] not in panels and panels.count(["B"]) == 2
    if workers > 1:
        assert ["C", "D"] not in panels
    forced = result["forced_choice"]
    detail = read_json(run / forced["artifact_ref"])
    assert forced["rejected_review_batches"] == 2 and forced["review_fallback_items"] == 1
    assert detail["review_recoveries"][0]["salvaged"] == ["A"]
    assert detail["review_recoveries"][1]["fallback"] == ["B"]
    assert set(detail["panels"]) == {"A", "C", "D"}
    assert any(g["label"] == "B" and not g["reliable"] for g in detail["groups"])
    replay(run, snap, None, None, tmp_path / "replayed", cfg)
    assert read_json(tmp_path / "replayed/artifacts/forced_choice.json")["groups"] == detail["groups"]


def test_rejected_relationship_repair_preserves_uncontested_parents_and_aliases():
    from rarecell.review_recovery import conservative_relations
    from rarecell.single_choice import validate_panels
    specs = {name: dict(canonical=name, panels=[["G"]]) for name in ("T", "MAIT", "alias", "X", "Y")}
    specs["alias"]["canonical"] = "MAIT"
    specs["MAIT"]["broader_labels"] = ["T"]
    specs["X"]["broader_labels"] = ["Y", "T"]
    specs["Y"]["broader_labels"] = ["X"]
    mapping = conservative_relations(specs)
    fixed = {n: dict(spec, **mapping[n]) for n, spec in specs.items()}
    validate_panels(specs, fixed)
    assert fixed["alias"]["canonical"] == "MAIT" and fixed["MAIT"]["broader_labels"] == ["T"]
    assert fixed["X"]["broader_labels"] == ["T"] and not fixed["Y"]["broader_labels"]
