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
