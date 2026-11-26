"""Plan independent model calls. The controller owns all checkpoint and expression writes."""

from itertools import islice

from .state import BudgetExhausted, IntegrityError, QueryBlocked


def independent_tasks(route, width):
    pending = route.pending_seeds()
    if pending:
        for _ in range(width):
            if not pending:
                break
            context = route.triage_context(pending)
            yield None, context
            pending = pending[len(context["seeds"]):]
        return
    frozen = {r["job"] for r in route.s.all("RouteFreeze")}
    # Reserve the last few query slots for applying paid answers in order.
    room = route.ctrl.query_limit - route.s.usage()["new_queries"]
    if room < width:
        width = max(1, room)
    available = (j for j in route.jobs() if j["id"] not in frozen)
    for job in islice(available, width):
        yield job, None


def planned_turn(route, turn, pause_after):
    state, ctrl = route.s, route.ctrl
    planned = state.get_optional("RouteTurnPlan", str(turn))
    if planned:
        return planned
    workers = ctrl.cfg.get("routed_workers", 1)
    end = ctrl.budget["max_agent_turns"]
    if pause_after is not None:
        end = min(end, pause_after)
    width = min(workers, end - turn + 1)
    if width < 2 or ctrl.budget.get("max_reported_tokens"):
        return None
    mid = f"{state.get_meta('run_id')}:T{turn:04d}"
    if state.get_optional("Message", mid) or state.get_optional("ProviderDispatch", mid + "/api/1"):
        # Resume any turn allocated by the serial controller under its original routing.
        return None

    tasks = list(independent_tasks(route, width))
    if len(tasks) < 2:
        return None
    ready = []
    for job, ctx in tasks:
        if (ready and job and job["kind"] == "seed"
                and ctrl.query_limit - state.usage()["new_queries"] <= len(ready)):
            break
        try:
            scopes = route.scope_members(job) if job else {}
            ctx = route.context_for(job, ctx, scopes)
        except (BudgetExhausted, QueryBlocked):
            # The first task is handled by the normal serial stop policy.
            break
        if ctx is None:
            break
        ready.append((job, ctx))
    if len(ready) < 2:
        return None

    last_turn = turn + len(ready) - 1
    with state.transaction():
        for n, (job, ctx) in enumerate(ready, turn):
            plan = dict(job=job, context_ref=state.artifact(ctx), last_turn=last_turn)
            state.put("RouteTurnPlan", str(n), plan)
        state.event("routed_wave_planned", {"first_turn": turn, "last_turn": last_turn})
    return state.get("RouteTurnPlan", str(turn))


def wave_message(route, turn, planned, retry_unknown):
    from .run import obtain_parallel

    s = route.s
    run_id = s.get_meta("run_id")
    mid = f"{run_id}:T{turn:04d}"
    msg = s.get_optional("Message", mid)
    if msg is None:
        # The persisted plan controls membership even if the execution width changes on resume.
        width = route.ctrl.cfg.get("routed_workers", 1)
        stop = min(planned["last_turn"] + 1, turn + width)
        requests = [(n, s.read_artifact(s.get("RouteTurnPlan", str(n))["context_ref"]))
                    for n in range(turn, stop)
                    if s.get_optional("Message", f"{run_id}:T{n:04d}") is None]
        obtain_parallel(route.ctrl, requests, retry_unknown)
        msg = s.get("Message", mid)
    if "artifacts/" + msg["context_hash"] + ".json" != planned["context_ref"]:
        raise IntegrityError("Saved provider message does not belong to its routed turn plan")
    return msg
