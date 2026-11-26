"""Concurrent worker I/O and single-writer checkpoints, without external calls."""
from contextlib import closing
from email.message import Message
import io
import json
import threading
from types import SimpleNamespace
from urllib.error import HTTPError

import pytest

from rarecell import providers
from rarecell.run import obtain, obtain_parallel
from rarecell.state import State


def test_parallel_persists_siblings_and_reuses_paid_turns(tmp_path, monkeypatch):
    entered = threading.Barrier(3)
    calls = []
    main_thread = threading.get_ident()

    def complete(context, cfg):
        name = context["target"]
        assert threading.get_ident() != main_thread
        calls.append(name)
        if len(calls) < 4:
            entered.wait(timeout=2)
        if name == "A" and calls.count("A") == 1:
            raise providers.ProviderError("Interrupted after dispatch", "unknown_execution")
        if name == "B":
            exc = providers.ProviderError("Truncated answer", "incomplete")
            exc.usage = {"total_tokens": 123}
            raise exc
        return json.dumps({"action": "freeze", "target": name}), {"api_calls": 1}

    monkeypatch.setattr(providers, "complete", complete)
    requests = [(i + 1, {"target": name}) for i, name in enumerate("ABC")]
    cfg = dict(routed_workers=3, provider=dict(model="test", base_url="http://localhost", api_key_env=None))
    run, cache = tmp_path / "run", tmp_path / "cache"
    with closing(State(run, cache)) as state:
        state.set_meta("run_id", "r")
        ctrl = SimpleNamespace(state=state, cfg=cfg)
        with pytest.raises(providers.ProviderError, match="Interrupted after dispatch"):
            obtain_parallel(ctrl, requests, False)
        assert {item["message_id"] for item in state.all("Message")} == {"r:T0002", "r:T0003"}
        b = state.get("Message", "r:T0002")
        assert not b["complete"] and b["usage"]["total_tokens"] == 123
        assert len(state.all("ProviderDispatch")) == 3
        assert not state.get_optional("Message", "r:T0001")
    with closing(State(run, cache)) as state:
        ctrl = SimpleNamespace(state=state, cfg=cfg)
        with pytest.raises(providers.ProviderError, match="Prior dispatch"):
            obtain_parallel(ctrl, requests, False)
        assert len(calls) == 3
        obtain_parallel(ctrl, requests, True)
        obtain_parallel(ctrl, requests, False)
        assert calls.count("A") == 2 and len(calls) == 4
        for turn, context in requests:
            assert obtain(ctrl, "live", turn, False, context=context) == state.get("Message", "r:T%04d" % turn)


@pytest.mark.parametrize("code,header,automatic,status", [
    ("insufficient_quota", "1", False, "blocked"),
    ("rate_limit_exceeded", None, True, "temporary_failure"),
    ("other", None, False, "temporary_failure"),
    ("other", "12", True, "temporary_failure"),
])
def test_429_classification(code, header, automatic, status):
    headers = Message()
    if header:
        headers["Retry-After"] = header
    content = json.dumps({"error": {"code": code}}).encode()
    with HTTPError("http://localhost", 429, "Too many requests", headers, io.BytesIO(content)) as response:
        exc = providers._http_failure(response)
    assert exc.status == status and hasattr(exc, "retry_after") == automatic
    if automatic:
        assert exc.retry_after == float(header or 0)


def test_transient_retry_is_bounded_and_shares_cooldown(monkeypatch):
    clock = [100.]
    sleeps = []
    def sleep(dt):
        sleeps.append(dt)
        clock[0] += dt
    monkeypatch.setattr(providers.time, "sleep", sleep)
    monkeypatch.setattr(providers.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(providers.random, "uniform", lambda a, b: .125)
    calls = []
    def complete(ctx, cfg):
        calls.append(clock[0])
        if len(calls) <= 3:
            exc = providers.ProviderError("HTTP 429", "temporary_failure")
            exc.retry_after = 5
            raise exc
        return "ok", {"api_calls": 1, "total_tokens": 10}
    monkeypatch.setattr(providers, "complete", complete)
    window = providers.RateLimitWindow()
    with pytest.raises(providers.ProviderError) as failed:
        window.complete({}, {})
    assert calls == [100., 105.125, 110.25] and failed.value.usage["api_calls"] == 3
    assert window.complete({}, {})[1]["total_tokens"] == 10
    assert calls[-1] == 115.25 and sleeps[-1] == 5.
