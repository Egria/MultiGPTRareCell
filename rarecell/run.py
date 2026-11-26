"""Checkpointed, single-writer discovery policy loop."""
import fcntl
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import nullcontext
import json
import os
import shutil
import tempfile
import time
import uuid
from pathlib import Path

import yaml

from . import routed, __version__, exhaustive, pool_followup, expansion_review, providers
from .annotation import verify_annotation
from .config import ACTIVE_MODES, ANNOTATED_MODES, MODES, ROUTED_BUDGET_STOPS, ROUTED_MODE, RESOURCES, atomic, contract, fingerprint, hash_file, read_json, validate_routed_workers, validate_settings, write_json
from .discovery import Controller, SCIENTIFIC, prompt_version
from .finalize import export
from .prepare import verify_snapshot
from .protocol import VERSION as PROTOCOL_VERSION, ProtocolError, parse_message


def code_fingerprint():
    module = Path(__file__).parent
    files = ("discovery.py", "run.py", "protocol.py", "science.py", "state.py", "providers.py",
             "config.py", "seed_review.py", "exhaustive.py", "adaptive.py", "base_evidence.py", "finalize.py", "annotation.py",
             "contrast_panels.py", "singleton_markers.py", "pool_followup.py", "nb_mixture.py", "claim_resolution.py", "handoff.py", "expansion_review.py", "routed.py", "rule_search.py",
             "single_choice.py", "forced_routing.py", "lineage_review.py", "unclaimed_review.py", "routed_parallel.py")
    return fingerprint({name: hash_file(module / name) for name in files})


def disjoint(a, b):
    return a != b and a not in b.parents and b not in a.parents


def new_spec(snapshot, run, cache, annotations, cfg, mode, policy, replication_id):
    if cfg.get("exhaustive_seed_review") and (policy != "live" or mode not in ACTIVE_MODES):
        raise ValueError("Exhaustive seed review requires live A2, A3, or A3_adaptive")
    if mode == ROUTED_MODE and policy != "live":
        raise ValueError("A3_routed requires live discovery")
    if policy not in {"live", "fixed", "scripted_demo"}:
        raise ValueError("Unknown policy")
    if policy == "fixed" and mode != "A3":
        raise ValueError("The matched fixed-policy control uses A3 with shared annotations")
    if policy == "scripted_demo" and cfg["provider"]["kind"] != "scripted_demo":
        raise ValueError("Scripted messages require an explicit demo configuration")
    if mode not in MODES:
        raise ValueError("Unknown discovery mode")
    snapshot, run, cache = map(lambda p: Path(p).resolve(), (snapshot, run, cache))
    if not disjoint(snapshot, run) or not disjoint(cache, run) or not disjoint(snapshot, cache):
        raise ValueError("Use separate snapshot, run and cache directories")
    if mode in ANNOTATED_MODES:
        annotations = Path(annotations or run.parent / "annotation").resolve()
        if any(not disjoint(annotations, x) for x in (snapshot, run, cache)):
            raise ValueError("Use a separate annotation directory")
    else:
        annotations = None
    return {"snapshot": str(snapshot), "cache": str(cache), "annotations": str(annotations) if annotations else None,
            "mode": mode, "policy": policy, "run_id": uuid.uuid4().hex, "replication_id": replication_id,
            "software": __version__, "code_sha256": code_fingerprint(), "prompt_version": prompt_version(mode, cfg.get("exhaustive_seed_review", False)), "cfg_hash": fingerprint(cfg), "contract_hash": fingerprint(contract())}


def fixed_priority(ctrl, c):
    kinds = {str(src.get("kind", "")).lower() for src in c["sources"]}
    priority = 0 if any("outlier" in k for k in kinds) else 1 if len(c["members"]) == 1 else 2
    return (priority, *ctrl.priority(c))


def offline_message(ctrl, strategy, turn):
    """The scripts know only the controller's discovered objects and run budgets."""
    if strategy == "scripted_demo" and turn == 1:
        return "Candidate inventory\n```json\n{\"operation\":\"list_candidates\",\"limit\":\"4\",}\n```"
    investigated = {q["candidate_id"] for q in ctrl.state.all("Investigation")}
    initial = sorted((c for c in ctrl.candidates.values() if c["initial"]), key=lambda c: fixed_priority(ctrl, c))
    queue = [c for c in initial if c["id"] not in investigated]
    remaining = ctrl.budget["max_new_queries"] - ctrl.state.usage()["new_queries"]
    room = min(ctrl.budget["max_actions_per_turn"], remaining, 2 if strategy == "scripted_demo" else remaining)
    can_expand = ctrl.active and room > 0 and turn < ctrl.budget["max_agent_turns"]
    if strategy == "scripted_demo" and turn > 2:
        can_expand = False
    if can_expand and queue:
        return json.dumps([{"action": "expand_candidate", "candidate_id": c["id"]} for c in queue[:room]])
    expanded = sorted([c for c in ctrl.candidates.values() if any(x.get("kind") == "expansion" for x in c["sources"])],
                      key=lambda c: fixed_priority(ctrl, c))
    return json.dumps({"action": "submit", "candidate_ids": [c["id"] for c in expanded],
                       "reason": "Submit threshold-qualified expansions. Competing assignments remain unresolved."})


def dispatch_id(state, mid, retry_unknown):
    first, second = mid + "/api/1", mid + "/api/2"
    if state.get_optional("ProviderDispatch", first) is None:
        return first
    last = second if state.get_optional("ProviderDispatch", second) is not None else first
    if last == second or not retry_unknown:
        failure = state.get_optional("ProviderFailure", last)
        raise providers.ProviderError("Prior dispatch has no usable response; bounded explicit retry available once",
                                      failure["status"] if failure else "unknown_execution")
    return second


def record_dispatch(state, mid, aid, context):
    ref = state.artifact(context)
    state.put("ProviderDispatch", aid, {"id": aid, "message_id": mid, "context_ref": ref})


def record_reply(state, mid, context, raw, usage, complete=True):
    item = {"message_id": mid, "complete": complete, "raw": raw, "usage": usage, "context_hash": fingerprint(context)}
    if context.get("requested_evidence") and complete:
        item["delivered_receipts"] = [r["receipt_id"] for r in context["requested_evidence"]]
    state.put("Message", mid, item)
    return item


def record_failure(state, mid, aid, context, exc):
    failure = {"id": aid, "status": exc.status, "detail": str(exc)}
    if hasattr(exc, "usage"):
        failure["usage"] = exc.usage
    with state.transaction():
        state.put("ProviderFailure", aid, failure)
        if exc.status == "incomplete":
            # Preserve paid unfinished text as a consumed turn, without parsing a partial action.
            usage = getattr(exc, "usage", {"provider": "openai_compatible", "api_calls": 1})
            return record_reply(state, mid, context, "", usage, complete=False)
    return None


def check_provider(cfg):
    health = providers.diagnose(cfg)
    if health["status"] == "blocked":
        raise providers.ProviderError(health["detail"])


def obtain(ctrl, policy, turn, retry_unknown, *, context=None):
    state = ctrl.state
    mid = "%s:T%04d" % (state.get_meta("run_id"), turn)
    saved = state.get_optional("Message", mid)
    if saved is not None:
        return saved
    if context is None:
        context = ctrl.context(turn)
    if policy != "live":
        raw = offline_message(ctrl, policy, turn)
        return record_reply(state, mid, context, raw, {"provider": policy, "api_calls": 0})
    check_provider(ctrl.cfg)
    aid = dispatch_id(state, mid, retry_unknown)
    record_dispatch(state, mid, aid, context)
    try:
        raw, usage = providers.complete(context, ctrl.cfg)
    except providers.ProviderError as exc:
        incomplete = record_failure(state, mid, aid, context, exc)
        if incomplete is None:
            raise
        return incomplete
    return record_reply(state, mid, context, raw, usage)


def obtain_parallel(ctrl, requests, retry_unknown):
    """Fetch one wave of independent paid replies. Only this caller touches SQLite."""
    state = ctrl.state
    run_id = state.get_meta("run_id")
    pending = []
    # Resolve ambiguous previous attempts before paying for another request in the wave.
    for turn, context in requests:
        mid = "%s:T%04d" % (run_id, turn)
        if state.get_optional("Message", mid) is None:
            pending.append((mid, dispatch_id(state, mid, retry_unknown), context))
    if not pending:
        return
    if len({mid for mid, _, _ in pending}) != len(pending):
        raise ValueError("A parallel wave must use distinct turn IDs")
    check_provider(ctrl.cfg)
    if not hasattr(ctrl, "provider_window"):
        ctrl.provider_window = providers.RateLimitWindow()
    limit = min(len(pending), ctrl.cfg.get("routed_workers", 1))
    with ThreadPoolExecutor(max_workers=limit) as workers:
        # A resumed wave may exceed a newly lowered worker limit. Journal only the calls
        # we are about to start, so a stopped cohort leaves the next one undispatched.
        for offset in range(0, len(pending), limit):
            futures = {}
            failures = {}
            cohort = pending[offset:offset + limit]
            for mid, aid, context in cohort:
                record_dispatch(state, mid, aid, context)
                future = workers.submit(ctrl.provider_window.complete, context, ctrl.cfg)
                futures[future] = (mid, aid, context)
            for future in as_completed(futures):
                mid, aid, context = futures[future]
                try:
                    raw, usage = future.result()
                except providers.ProviderError as exc:
                    if record_failure(state, mid, aid, context, exc) is None:
                        failures[mid] = exc
                except Exception as exc:
                    # Keep any other paid replies even if a worker raises unexpectedly.
                    failures[mid] = exc
                else:
                    record_reply(state, mid, context, raw, usage)
            for mid, _, _ in cohort:
                if mid in failures:
                    raise failures[mid]


def repaired_locally(raw, parsed):
    if any("error" in p for p in parsed):
        return False
    try:
        json.loads(raw)
    except (ValueError, TypeError):
        return True
    return False


def follow_browsing_with_seed(ctrl, mid):
    """Use one spare local action on a screened seed. Journal the choice before computation."""
    s = ctrl.state
    done = s.get_optional("SeedFollowup", mid)
    if done is not None:
        return done
    plan = s.get_optional("SeedFollowupPlan", mid)
    if plan is None:
        if (len(s.all("SeedFollowupPlan")) >= ctrl.budget.get("max_seed_followups", 3) or
                s.usage()["new_queries"] >= ctrl.budget["max_new_queries"] or
                not ctrl.budget["max_expansion_depth"]):
            return None
        queue = ctrl.pending_seed_reviews()
        if not queue:
            return None
        plan = dict(queue[0], source="controller_seed_review")
        s.put("SeedFollowupPlan", mid, plan)
    cid = plan["candidate_id"]
    try:
        result = ctrl.scientific("expand_candidate", {"candidate_id": cid})
    except (ValueError, KeyError) as exc:
        result = {"status": "blocked", "error": str(exc)[:500]}
    record = {"source": plan["source"], "anchor_id": cid, "peer": plan["peer"],
              "shared_markers": plan["shared_markers"], "views": plan["views"], "result": result}
    s.put("SeedFollowup", mid, record)
    return record


def process_message(ctrl, message, turn):
    s = ctrl.state
    mid = message["message_id"]
    seed_batch = s.get_optional("SeedBatch", mid) if ctrl.exhaustive else None
    pool_batch = s.get_optional("PoolBatch", mid) if ctrl.exhaustive and ctrl.labeled else None
    boundary_batch = s.get_optional("ExpansionBatch", mid) if ctrl.exhaustive and ctrl.adaptive else None
    parser_version = PROTOCOL_VERSION + ("/seed" if seed_batch else "")
    parse_id = fingerprint([parser_version, message["raw"], message.get("complete", True)])
    interpretation = s.get_optional("ParsedMessage", parse_id)
    if interpretation is None:
        parsed = parse_message(message["raw"], origin="agent", complete=message.get("complete", True), seed_batch=bool(seed_batch))
        with s.transaction():
            s.put("ParsedMessage", parse_id, {"id": parse_id, "commands": parsed})
            count = int(repaired_locally(message["raw"], parsed))
            s.set_meta("format_repairs", s.get_meta("format_repairs", 0) + count)
    else:
        parsed = interpretation["commands"]
    limit = ctrl.budget["max_actions_per_turn"]
    closing = () if seed_batch or pool_batch or boundary_batch else ctrl.closing_actions(turn)
    if closing and s.get_meta("decision_start_turn") is None:
        s.set_meta("decision_start_turn", turn)
    if not closing and len(parsed) > limit:
        parsed = parsed[:limit] + [{"error": "Remaining batch actions exceed max_actions_per_turn"}]
    eligible = 0
    receipts = []
    progress = False
    for index, command in enumerate(parsed):
        opid = "%s/%d" % (mid, index)
        if closing and command.get("action") in closing:
            eligible += 1
        old = s.get_optional("Receipt", opid)
        if old is not None:
            receipts.append(old["result"])
            progress |= old["progress"]
            continue
        # A revision/submit and its receipt commit together. Local numeric work
        # has an independent cache-publication journal and can be reattached.
        local_query = command.get("action") in SCIENTIFIC or command.get("action") == "investigate_pool"
        outer = nullcontext() if local_query else s.transaction()
        with outer:
            if s.get_meta("submitted", False):
                result = {"error": "Run already submitted; subsequent action blocked", "status": "blocked"}
            else:
                try:
                    if boundary_batch and command.get("action") and command["action"] != "review_expansion":
                        raise ValueError("Review the measured additions for the assigned pools first")
                    if pool_batch and command.get("action") and command["action"] != "investigate_pool":
                        raise ValueError("Pool investigation is open: choose one test per assigned pool")
                    if seed_batch and command.get("action") and command["action"] != "review_seed":
                        raise ValueError("Written seed review is open: use review_seed for the listed seeds first")
                    if closing and command.get("action") and command["action"] not in closing:
                        raise ValueError("Decision phase accepts only: " + ", ".join(closing))
                    if closing and eligible > limit:
                        raise ValueError("Remaining decision actions exceed max_actions_per_turn")
                    if (closing and command.get("action") == "submit" and "candidate_ids" not in command["args"]
                            and not s.get_meta("selected", [])):
                        raise ValueError("No saved selection. Submit candidate_ids explicitly, including [] to decline all.")
                    result = ctrl.execute(command, mid)
                except (ValueError, KeyError) as exc:
                    clarification = "error" in command or isinstance(exc, ProtocolError)
                    result = {"error": str(exc)[:500], "status": "needs_clarification" if clarification else "blocked",
                              "action": command.get("action")}
            logical = fingerprint([command.get("action"), command.get("args"), result.get("query_key"),
                                   result.get("evidence_id"), result.get("candidate_ids"), result.get("candidate_id")])
            changed = False
            if "error" not in result and result.get("status") != "over_limit":
                if command.get("action") == "propose_revision":
                    # Restating a rationale does not change the measured evidence,
                    # member assignments, labels or statuses.
                    projection = set()
                    for cid in s.get_meta("selected", []):
                        c = ctrl.candidate(cid)
                        observations = sorted(e for e in c["evidence_ids"] if s.get("Evidence", e)["kind"] == "observation")
                        projection.add(fingerprint([c["set_hash"], c["proposed_label"], c["status"], observations]))
                    logical = fingerprint(["selected", sorted(projection)])
                elif command.get("action") == "select_boundary":
                    logical = fingerprint(["draft_choice", result["selected_ids"]])
                elif command.get("action") == "apply_boundary":
                    logical = fingerprint(["boundary", result["evidence_id"], result["membership_hash"]])
                elif command.get("action") in {"get_evidence", "list_candidates", "list_ontology_pools", "preview_submission"}:
                    # Reading rows is distinct from computing their underlying observation.
                    # Equal pages also stay equal when expressed with different paging arguments.
                    page = {k: v for k, v in result.items() if k not in {"offset", "limit"}}
                    logical = fingerprint([command["action"], page])
                elif result.get("query_key"):
                    logical = fingerprint(["query", result["query_key"]])
                changed = s.get_optional("Progress", logical) is None
            with s.transaction():
                if ctrl.adaptive and command.get("action") == "get_evidence" and "sections" in result:
                    args = command["args"]
                    entry = {"evidence_id": result["evidence_id"], "metric": args.get("metric", "default"),
                             "offset": result["offset"], "limit": result["limit"]}
                    pages = s.get_meta("evidence_pages", [])
                    s.set_meta("evidence_pages", ([p for p in pages if p != entry] + [entry])[-5:])
                if changed:
                    s.put("Progress", logical, {"id": logical})
                if "error" in result:
                    metric = "protocol_errors" if result.get("status") == "needs_clarification" else "action_errors"
                    s.set_meta(metric, s.get_meta(metric, 0) + 1)
                    if str(result["error"]).startswith("max_new_queries"):
                        s.set_meta("denied_budget", True)
                s.put("Receipt", opid, {"id": opid, "result": result, "progress": changed,
                                        "metadata": command.get("metadata", {}), "message_id": mid})
                s.event("action", {"message_id": mid, "index": index, "action": command.get("action"),
                                   "status": result.get("status", "ok"), "receipt": opid, "changed": changed})
        receipts.append(result)
        progress |= changed
    if boundary_batch:
        expansion_review.finish_batch(ctrl, mid, boundary_batch)
    if pool_batch:
        pool_followup.finish_batch(ctrl, mid, pool_batch)
    browse_only = parsed and all(c.get("action") == "list_candidates" for c in parsed)
    if (not ctrl.exhaustive and not closing and ctrl.active and s.get_meta("policy") == "live" and browse_only
            and len(parsed) < limit and any("candidates" in r for r in receipts)):
        followup = follow_browsing_with_seed(ctrl, mid)
        if followup:
            receipts.append(dict(followup["result"], source=followup["source"], anchor_id=followup["anchor_id"],
                                 peer=followup["peer"], shared_markers=followup["shared_markers"]))
            progress |= "query_key" in followup["result"]
    with s.transaction():
        if not seed_batch and not pool_batch and not boundary_batch:
            delivered = set(message.get("delivered_receipts", []))
            pending = [rid for rid in s.get_meta("requested_receipts", []) if rid not in delivered]
            focus = s.get_meta("recent_followup_ids", [])
            for i, (command, result) in enumerate(zip(parsed, receipts)):
                if result.get("pool_followup") and command.get("action") == "get_evidence":
                    cid = result["pool_followup"]["candidate_id"]
                    focus = [cid] + [entry for entry in focus if entry != cid]
                if command.get("action") in SCIENTIFIC | {"get_evidence"} and "error" not in result:
                    rid = "%s/%d" % (mid, i)
                    if rid not in pending:
                        pending.append(rid)
            s.set_meta("recent_followup_ids", focus[:8])
            s.set_meta("requested_receipts", pending)
        s.set_meta("receipts", receipts[-12:])
        s.set_meta("finished_turn", turn)
        stall_counter = "seed_review_stall" if seed_batch else "stall"
        first_exposure = seed_batch and seed_batch["includes_new"]
        if pool_batch or boundary_batch:
            s.set_meta("pool_review_turns", s.get_meta("pool_review_turns", 0) + 1)
        else:
            s.set_meta(stall_counter, 0 if progress or first_exposure else s.get_meta(stall_counter, 0) + 1)
        if seed_batch:
            s.set_meta("seed_review_turns", turn)
        s.event("turn", {"message_id": mid, "turn": turn, "progress": progress})


def retain_freeze(run, frozen):
    """An evaluated partial checkpoint remains independently scoreable after resume."""
    archive = run / "history" / hash_file(run / "frozen.json")
    if archive.is_dir():
        return
    archive.parent.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(dir=archive.parent, prefix=".retain-") as tmp:
        stage = Path(tmp) / "export"
        stage.mkdir()
        for name in list(frozen["files"]) + ["frozen.json"]:
            src, dst = run / name, stage / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            try:
                os.link(src, dst)
            except OSError:
                shutil.copyfile(src, dst)
        stage.rename(archive)


def explore(ctrl, policy, pause_after=None, retry_unknown=False):
    """Run bounded discovery turns and publish a saved draft on a terminal stall."""
    s = ctrl.state
    s.set_meta("policy", policy)
    status, detail = "partial", "max_agent_turns"
    review_limit = (2 if ctrl.labeled else 1) * len(ctrl.seed_ids) + ctrl.budget["stall_turns"] if ctrl.exhaustive else 0
    # Each seed can have a named singleton and a conflicting expansion in addition to the CL pools.
    pool_limit = len(ctrl.pool_candidates) + 2 * len(ctrl.seed_ids) if ctrl.exhaustive and ctrl.labeled else 0
    hard_limit = ctrl.budget["max_agent_turns"] + review_limit + pool_limit * (3 if ctrl.adaptive else 1)
    for turn in range(s.get_meta("finished_turn", 0) + 1, hard_limit + 1):
        mid = "%s:T%04d" % (s.get_meta("run_id"), turn)
        queue = []
        review = False
        if ctrl.exhaustive and not s.get_meta("seed_phase_complete", False):
            queue = exhaustive.pending(ctrl)
            # A crash can leave every note committed just before the turn checkpoint.
            review = bool(queue) or s.get_optional("SeedBatch", mid) is not None
            if review and (turn > review_limit or s.get_meta("seed_review_stall", 0) >= ctrl.budget["stall_turns"]):
                s.set_meta("current_issue", "Written seed assessments are still missing after the review allowance.")
                return "partial", "seed_review_incomplete"
            if not review:
                if exhaustive.coverage(ctrl)["measurement_blocked"]:
                    s.set_meta("current_issue", "One or more required seed measurements failed. See seed_reviews.tsv.")
                    return "partial", "seed_measurement_blocked"
                s.set_meta("seed_phase_complete", True)
        # A fresh invocation reconstructs context or replays a persisted full message.
        cap = ctrl.budget.get("max_reported_tokens", 0)
        # Budget gates the next dispatch. A reply checkpointed before an
        # interruption is already paid for and must still be applied.
        if cap and policy == "live" and s.get_optional("Message", mid) is None:
            token_counts = [m["usage"].get("total_tokens") for m in s.all("Message")]
            if any(n is None for n in token_counts) or sum(token_counts) >= cap:
                detail = "token_usage_unreported" if any(n is None for n in token_counts) else "max_reported_tokens"
                break
        try:
            if review:
                exhaustive.prepare_batch(ctrl, mid, queue)
            elif ctrl.exhaustive and ctrl.labeled:
                review = pool_followup.prepare_batch(ctrl, mid)
                if not review and ctrl.adaptive:
                    review = expansion_review.prepare_batch(ctrl, mid)
            if not review and turn > ctrl.turn_limit:
                break
            if not review and ctrl.stall_review():
                closing_turn = s.get_meta("closing_turn", turn)
                if turn > closing_turn:
                    detail = "stalled"
                    break
                s.set_meta("closing_turn", closing_turn)
            message = obtain(ctrl, policy, turn, retry_unknown)
        except providers.ProviderError as exc:
            detail = "provider_blocked" if exc.status == "blocked" else exc.status
            s.set_meta("current_issue", str(exc))
            s.event(detail, {"reason": str(exc), "turn": turn})
            break
        except ValueError as exc:
            detail = "context_blocked"
            s.set_meta("current_issue", str(exc))
            s.event(detail, {"reason": str(exc), "turn": turn})
            break
        closing = not review and ctrl.stall_review()
        process_message(ctrl, message, turn)
        if review:
            if pause_after is not None and turn >= pause_after:
                return "partial", "paused"
            continue
        if s.get_meta("submitted", False):
            status, detail = "complete", "submitted"
            break
        if closing:
            detail = "stalled"
            break
        if s.get_meta("stall", 0) >= ctrl.budget["stall_turns"]:
            if not ctrl.stall_review() or turn == ctrl.turn_limit:
                detail = "stalled"
                break
        if pause_after is not None and turn >= pause_after and (policy != "live" or turn < ctrl.turn_limit):
            detail = "paused"
            break
    # A stalled run may use one remaining turn to choose among held hypotheses.
    # The controller can publish only the actual saved selection when that turn
    # produces no submission, or when the call/token budget precludes the turn.
    terminal_stall = detail == "stalled" or (
        ctrl.stall_review() and detail in {"max_agent_turns", "max_reported_tokens", "token_usage_unreported"})
    budget_stop = detail in {"max_agent_turns", "max_reported_tokens", "token_usage_unreported"}
    if (policy == "live" or ctrl.adaptive) and (terminal_stall or budget_stop) and s.get_meta("selected", []) and not s.get_meta("submitted", False):
        source = "controller_stall" if terminal_stall else "controller_budget"
        with s.transaction():
            outcome = ctrl.decide("submit", {"status": "provisional", "reason":
                "Publish the saved working selection after " + detail + "."}, None, source=source)
            s.set_meta("closeout_trigger", detail)
            s.event("provisional_submission", {"decision_id": outcome["decision_id"], "trigger": detail})
    if s.get_meta("submission_source") == "controller_stall":
        return "partial", "stalled_provisional"
    if s.get_meta("submission_source") == "controller_budget":
        return "partial", detail
    if (policy == "live" and not s.get_meta("submitted", False) and detail in {"stalled", "max_agent_turns"}
            and s.get_meta("decision_start_turn") is not None):
        s.set_meta("closeout_trigger", detail)
        s.set_meta("current_issue", "Decision phase ended without an explicit submission or a saved selection.")
        return "partial", "decision_missing"
    return status, detail


def resume_limits(cfg, spec, previous, freeze_sha, requested):
    """Use absolute ceilings. Only the exhausted resource can reopen a terminal freeze."""
    budget_history = spec.get("budget_extensions", [])
    limits = {key: cfg["budget"][key] for key in ROUTED_BUDGET_STOPS.values()}
    if budget_history:
        limits = budget_history[-1]["limits"]
    if not requested:
        return limits, None
    if any(type(v) is not int or v < limits[k] for k, v in requested.items()):
        raise ValueError("Resume ceilings must be integers at least as high as the current cumulative limits")
    increased = dict(limits, **requested)
    if budget_history and increased == limits:
        return limits, None  # Reissuing the same ceilings after a crash or submission.
    resource = ROUTED_BUDGET_STOPS.get(previous.get("detail"))
    if (not freeze_sha or spec["mode"] != ROUTED_MODE or spec["policy"] != "live"
            or previous.get("submitted") or not resource):
        raise ValueError("Budget extensions require an unsubmitted A3_routed query- or turn-budget stop")
    if increased[resource] <= limits[resource]:
        raise ValueError("Increase the exhausted cumulative ceiling: " + resource)
    return increased, {"source_freeze_sha256": freeze_sha, "previous_limits": limits, "limits": increased}


def _execute(run, spec, cfg, pause_after=None, retry_unknown=False, resumed=False, requested_limits=None,
             routed_workers=None):
    started = time.monotonic()
    snapshot = Path(spec["snapshot"])
    manifest = verify_snapshot(snapshot)
    code = code_fingerprint()
    old_code = spec.get("code_sha256")
    migrate = False
    if old_code != code and resumed and spec["mode"] == ROUTED_MODE and spec["policy"] == "live":
        compatible = read_json(RESOURCES / "routed_checkpoint_upgrade.json")
        migrate = (old_code in [compatible["from"], *compatible.get("also_from", [])] and code == compatible["to"]
                   and (run / "state.sqlite").is_file() and not (run / ".publishing.json").exists())
    if (fingerprint(cfg) != spec["cfg_hash"] or fingerprint(contract()) != spec["contract_hash"]
            or spec["software"] != __version__ or (old_code != code and not migrate)):
        raise ValueError("Checkpoint configuration, contract or method version changed")
    workers = spec.get("routed_workers", cfg.get("routed_workers", 1))
    if routed_workers is not None:
        validate_routed_workers(routed_workers, spec["mode"])
    changed_workers = routed_workers is not None and routed_workers != workers
    if changed_workers:
        workers = routed_workers
    if (run / "manifest.json").exists():
        if read_json(run / "manifest.json")["snapshot_hash"] != manifest["snapshot_hash"]:
            raise ValueError("Prepare snapshot differs from the pinned run")
    else:
        write_json(run / "manifest.json", manifest)
    publishing = run / ".publishing.json"
    previous = {}
    frozen = None
    freeze_sha = None
    if (run / "frozen.json").exists() and not publishing.exists():
        frozen = read_json(run / "frozen.json")
        if any(hash_file(run / f) != h for f, h in frozen["files"].items()):
            raise ValueError("Frozen run checksum mismatch")
        previous = read_json(run / "run_summary.json")
        freeze_sha = hash_file(run / "frozen.json")
    limits, extension = resume_limits(cfg, spec, previous, freeze_sha, requested_limits)
    budget_history = spec.get("budget_extensions", [])
    if frozen is not None:
        retry_blocked_pair = (resumed and spec["mode"] == ROUTED_MODE and spec["policy"] == "live"
                              and previous["detail"] == "routed_action_blocked" and not previous["submitted"]
                              and (previous.get("issue") or "").startswith("Routed task pair/"))
        # The effective extension is checkpointed before executing the first new action.
        # A killed process can therefore restart with plain resume and the same old freeze.
        reopening = extension or (budget_history and budget_history[-1]["source_freeze_sha256"] == freeze_sha)
        if not previous["resumable"] and not retry_blocked_pair and not reopening:
            return previous
        if reopening and not (run / "state.sqlite").is_file():
            raise ValueError("The original state.sqlite is required to continue discovery")
        retain_freeze(run, frozen)
    ann, astatus, annotation_detail = [], "disabled_by_mode", ""
    if spec["mode"] in ANNOTATED_MODES:
        try:
            if not cfg["annotation"]["enabled"] or not spec["annotations"]:
                raise ValueError("Labeled mode requires an enabled, verified OntoAnno directory")
            ann_dir = Path(spec["annotations"])
            ann = verify_annotation(snapshot, ann_dir, _manifest=manifest, allow_partial=True)
            astatus = read_json(ann_dir / "annotation_status.json")["status"]
        except (OSError, KeyError, ValueError) as exc:
            astatus, annotation_detail = "blocked", str(exc)
    if budget_history or extension:
        cfg = dict(cfg, budget=dict(cfg["budget"], **limits))
    if spec["mode"] == ROUTED_MODE:
        cfg = dict(cfg, routed_workers=workers)
    ctrl = Controller(snapshot, run, spec["cache"], cfg, spec["mode"], ann)
    s = ctrl.state
    try:
        if (migrate or extension) and (s.get_meta("run_id") != spec["run_id"] or not s.get_meta("initialized", False)
                        or astatus not in {"success", "partial"} or s.get_meta("annotation_hash") != fingerprint(ann)):
            raise ValueError("Cannot upgrade routed checkpoint: original state or verified annotations differ")
        s.set_meta("run_id", spec["run_id"])
        s.set_meta("replication_id", spec["replication_id"])
        if resumed:
            s.set_meta("resume_count", s.get_meta("resume_count", 0) + 1)
        ctrl.initialize()
        if astatus in {"success", "partial"}:
            # The pending initialisation of a previously blocked arm needs its label projection.
            if not s.all("Annotation") and ann:
                ctrl.state.set_meta("initialized", False)
                ctrl.candidates.clear()
                ctrl.by_set.clear()
                ctrl.initialize()
            annotation_hash = fingerprint(ann)
            if s.get_meta("annotation_hash", annotation_hash) != annotation_hash:
                raise ValueError("Verified annotation contents changed after initialization")
            s.set_meta("annotation_hash", annotation_hash)
        if migrate:
            upgrade = {"from": old_code, "to": code, "after_turn": s.get_meta("finished_turn", 0)}
            history = spec.get("code_upgrades", [spec["code_upgrade"]] if spec.get("code_upgrade") else [])
            spec = dict(spec, code_sha256=code, code_upgrade=upgrade, code_upgrades=[*history, upgrade])
        if extension:
            extension.update(after_turn=s.get_meta("finished_turn", 0), used_queries=s.usage()["new_queries"])
            spec = dict(spec, budget_extensions=[*budget_history, extension])
        if changed_workers:
            spec = dict(spec, routed_workers=workers)
        if migrate or extension or changed_workers:
            write_json(run / "checkpoint.json", spec)
        if ctrl.routed and workers != s.get_meta("routed_workers", 1):
            with s.transaction():
                s.event("routed_workers_changed", {"previous": s.get_meta("routed_workers", 1),
                                                  "workers": workers, "after_turn": s.get_meta("finished_turn", 0)})
                s.set_meta("routed_workers", workers)
        if spec.get("budget_extensions") and s.get_meta("budget_extensions") != spec["budget_extensions"]:
            with s.transaction():
                s.set_meta("budget_extensions", spec["budget_extensions"])
                s.event("budget_extended", spec["budget_extensions"][-1])
        if spec.get("code_upgrade") and s.get_meta("code_upgrade") != spec["code_upgrade"]:
            with s.transaction():
                s.set_meta("code_upgrade", spec["code_upgrade"])
                s.set_meta("code_upgrades", spec.get("code_upgrades", [spec["code_upgrade"]]))
                s.event("checkpoint_upgraded", spec["code_upgrade"])
        if publishing.exists():
            pending = read_json(publishing)
            result = export(ctrl, pending["status"], pending["detail"], spec["policy"], pending["annotation_status"])
            publishing.unlink()
            return result
        status, detail = "partial", "max_agent_turns"
        s.set_meta("current_issue", None)
        if astatus == "blocked":
            s.set_meta("current_issue", annotation_detail)
            detail = "annotation_blocked"
            s.event(detail, {"reason": annotation_detail})
        elif s.get_meta("submission_source") == "controller_stall":
            status, detail = "partial", "stalled_provisional"
        elif s.get_meta("submission_source") == "controller_budget":
            status, detail = "partial", s.get_meta("closeout_trigger")
        elif s.get_meta("submitted", False):
            status, detail = "complete", "submitted"
        else:
            if ctrl.routed:
                status, detail = routed.RoutedDiscovery(ctrl).run(obtain, pause_after, retry_unknown)
            else:
                status, detail = explore(ctrl, spec["policy"], pause_after, retry_unknown)
        if detail == "submitted":
            coverage = pool_followup.coverage(ctrl) if ctrl.exhaustive and ctrl.labeled else {}
            exhausted = ctrl.active and s.usage()["new_queries"] >= ctrl.query_limit
            checked = {r["candidate_id"] for r in s.all("Investigation")}
            pending = any(c["initial"] and c["id"] not in checked for c in ctrl.candidates.values())
            if manifest.get("status") == "partial":
                status, detail = "partial", "submitted_with_blocked_prepare"
            elif astatus == "partial":
                status, detail = "partial", "submitted_with_partial_annotations"
            elif ctrl.routed and s.get_meta("submission_source") == "controller_routed":
                status, detail = "partial", "submitted_with_deferred_pairs"
            elif coverage.get("pending") or coverage.get("conflict_pending") or coverage.get("guided_pending"):
                status, detail = "partial", "submitted_with_pending_pool_followup"
            elif ctrl.exhaustive and exhaustive.coverage(ctrl)["named_seed_assessments"].get("unassessed"):
                status, detail = "partial", "submitted_with_incomplete_name_review"
            elif coverage.get("unreviewed_expansions") or coverage.get("expansion_reviews_pending"):
                status, detail = "partial", "submitted_with_unreviewed_expansions"
            elif coverage.get("blocked") or coverage.get("conflict_blocked") or coverage.get("guided_blocked"):
                status, detail = "partial", "submitted_with_blocked_pool_followup"
            elif s.get_meta("denied_budget", False) or (exhausted and pending and not ctrl.routed):
                status, detail = "partial", "submitted_after_budget_exhaustion"
            elif any(q["status"] != "success" for q in s.all("CanonicalQuery")):
                status, detail = "partial", "submitted_with_blocked_queries"
        s.set_meta("elapsed_seconds", s.get_meta("elapsed_seconds", 0) + time.monotonic() - started)
        s.event("export", {"status": status, "detail": detail})
        write_json(publishing, {"status": status, "detail": detail, "annotation_status": astatus})
        result = export(ctrl, status, detail, spec["policy"], astatus)
        publishing.unlink()
        return result
    finally:
        s.close()


def discover(snapshot, run, cache, cfg, mode=None, policy="live", annotations=None, replication_id=None, pause_after=None):
    mode = mode or cfg["mode"]
    cfg = validate_settings(dict(cfg, mode=mode))
    run = Path(run).resolve()
    spec = new_spec(snapshot, run, cache, annotations, cfg, mode, policy, replication_id)
    run.mkdir(parents=True, exist_ok=True)
    with (run / ".controller.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if any(p.name != ".controller.lock" for p in run.iterdir()):
            raise ValueError("Run directory is in use; use resume or a fresh directory")
        atomic(run / "config_resolved.yaml", yaml.safe_dump(cfg, sort_keys=True))
        atomic(run / "project_contract.yaml", (RESOURCES / "PROJECT_CONTRACT.yaml").read_bytes())
        write_json(run / "checkpoint.json", spec)
        return _execute(run, spec, cfg, pause_after)


def resume(run, retry_unknown=False, *, max_new_queries=None, max_agent_turns=None, routed_workers=None):
    run = Path(run).resolve()
    with (run / ".controller.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        spec = read_json(run / "checkpoint.json")
        cfg = yaml.safe_load((run / "config_resolved.yaml").read_text())
        requested = {k: v for k, v in dict(max_new_queries=max_new_queries, max_agent_turns=max_agent_turns).items()
                     if v is not None}
        if (requested or routed_workers is not None) and (run / ".publishing.json").exists():
            _execute(run, spec, cfg, resumed=True)  # Finish the interrupted export before archiving it.
            spec = read_json(run / "checkpoint.json")
        return _execute(run, spec, cfg, retry_unknown=retry_unknown, resumed=True, requested_limits=requested,
                        routed_workers=routed_workers)
