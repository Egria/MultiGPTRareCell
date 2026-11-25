"""Budgeted model exchanges and projection of immutable nominations into disjoint claims."""
from collections import defaultdict

from . import providers
from .config import fingerprint, json_bytes
from .routed import label, response
from .rule_search import ExpressionSearch
from .single_choice import build_groups, choose, measure, validate_panels, validate_proposal, validate_rulings


class StopReview(Exception):
    pass


def _batches(rows, prompt, limit, size):
    batch = {}
    for key, row in rows.items():
        trial = {**batch, key: row}
        if batch and (len(trial) > size or len(json_bytes(prompt(trial))) > limit):
            yield batch
            trial = {key: row}
        if len(json_bytes(prompt(trial))) > limit:
            raise ValueError("A single forced-choice measurement exceeds context_bytes")
        batch = trial
    if batch:
        yield batch


def adjudicate(ctrl, claims, ask):
    cells = list(ctrl.cells)
    ids = {cell: i for i, cell in enumerate(cells)}
    limit = ctrl.budget["context_bytes"]
    tissue = dict(species=ctrl.manifest.get("species"), tissue=ctrl.manifest.get("tissue"))
    indexed = {c["candidate_id"]: c for c in claims}
    quality = {cid: c.get("label_reliability") for cid, c in indexed.items()}
    unknown = {cid: {k: c[k] for k in ("label", "reason", "support")} for cid, c in indexed.items()
             if label(c["label"]) and quality[cid] not in ("supported", "uncertain")}
    def quality_prompt(batch):
        return dict(tissue, phase="conflict_label_audit", candidates=batch,
                    instruction='Assess the proposed identity of each frozen group using its rationale. Distinguish identity from purity or exact boundary. Return {"reviews":{candidate_id:"supported" or "uncertain"}}. A support number alone is insufficient.')
    for batch in _batches(unknown, quality_prompt, limit, 24):
        def validate_quality(reply):
            reviews = reply.get("reviews")
            if not isinstance(reviews, dict) or reviews.keys() != batch.keys() or any(v not in ("supported", "uncertain") for v in reviews.values()):
                raise ValueError("Assess every supplied group name")
        quality.update(ask(quality_prompt(batch), validate_quality)["reviews"])
    names = sorted({c["label"] for cid, c in indexed.items() if label(c["label"]) and quality[cid] == "supported"})
    def panel_prompt(batch):
        return dict(tissue, phase="conflict_panel_nomination", labels=list(batch), vocabulary=names,
                    instruction='Using only the class names and tissue, propose discriminating coexpression evidence. Return {"labels":{name:{"canonical":exact supplied equivalent label or itself,"panels":[[gene_symbols],...],"min_positive":2,"required_panel":1}}}. Use 2-3 nonoverlapping functional panels per class, with 1-4 genes each. required_panel is optional and 1-based. Use the same canonical for true synonyms; distinguish subtypes and states. For a particular rival, optionally supply versus:{rival_canonical_name:{panels:[[genes],...],min_positive:2}} inside the class proposal. Prioritize mechanisms distinguishing those two names. Each panel is scored in each cell and compared with the opposite flank.')
    specs = {}
    for batch in _batches(dict.fromkeys(names), panel_prompt, limit, 32):
        def validate_batch(reply):
            partial = reply.get("labels")
            if not isinstance(partial, dict) or partial.keys() != batch.keys():
                raise ValueError("Supply marker panels for each supplied name")
            # Validate the per-label shape even when its alias is specified in another batch.
            for name, s in partial.items():
                validate_proposal(name, s, names)
        specs.update(ask(panel_prompt(batch), validate_batch)["labels"])
    validate_panels(names, specs)
    groups = build_groups(claims, quality, specs, ids)
    expr = ExpressionSearch(ctrl)
    pairs, cards = measure(expr.rows, cells, expr.index, groups, specs)
    def rule_prompt(batch):
        panel_names = sorted({c[k] for c in batch.values() for k in ("left", "right")})
        return dict(tissue, phase="conflict_measured_review", panels={n: specs[n] for n in panel_names}, pairs=batch,
                    instruction='Freeze one decision rule per pair. pair_specific_panels in a measurement supersede the default class panels. Profiles encode left/right integer masks of positive panels, with panel 1 as bit 1. Each profile lists [cell count, unilateral support (1 left,-1 right,0 both/neither), median left and right panel scores]. Each small pair additionally gives every cell index, mask and both score vectors. Scores are log1p sums of normalized linear expression. Positive means detected and above the opposing P95. If the opposing flank is empty we use leave-one-out shared references. One sufficiently supported side wins; for other cases specify a default, or override an exact profile or displayed cell. Consider specificity, reference size and competing expression. Return {"rules":{pair_key:{"default":"geometry","profiles":{},"cells":{}}}}. Allowed choices left, right, geometry. Geometry requests anonymous expression proximity as the forced deciding vote.')
    rules = {}
    for batch in _batches(cards, rule_prompt, limit, 8):
        def validate_batch(reply):
            validate_rulings(batch, reply.get("rules"))
        rules.update(ask(rule_prompt(batch), validate_batch)["rules"])
    owners, votes, ties = choose(cells, groups, pairs, rules)
    by_group = defaultdict(list)
    for cell, gid in owners.items():
        by_group[gid].append(cell)
    result = dict(method="routed_single_choice/1", snapshot_hash=ctrl.snapshot_hash,
                frozen_claims={cid: fingerprint(sorted(c["members"])) for cid, c in indexed.items()},
                quality=quality, panels=specs, rules=rules, measurements=cards, pair_votes=votes,
                n_claimed=len(owners), n_wins_tied=ties,
                groups=[{**{k: g[k] for k in ("id", "parents", "label", "reliable")},
                         "n_nominated": len(g["members"]), "members": sorted(by_group[g["id"]])} for g in groups])
    return result


def project(ctrl, result):
    ref = ctrl.state.artifact(result)
    selected = []
    for g in result["groups"]:
        if not g["members"]:
            continue
        parents = g["parents"]
        evidence = sorted({e for i in parents for e in ctrl.candidate(i)["evidence_ids"]})
        basis = dict(method=result["method"], group=g["id"], result_ref=ref)
        c = ctrl.make_revision(g["members"], parents, g["label"], "provisional", evidence, "routed_claim", basis)
        selected.append(c["id"])
    ctrl.state.set_meta("forced_choice_ref", ref)
    return ctrl.decide("submit", dict(candidate_ids=selected, status="provisional",
                       reason="Publish disjoint measured memberships from the original frozen nominations."),
                       None, source="forced_single_choice")


def finish(route, obtain, pause_after, retry_unknown):
    ctrl, s = route.ctrl, route.s
    def ask(context, validate):
        key = fingerprint(context)
        cached = s.get_optional("ForcedReview", key)
        if cached:
            return cached["reply"]
        last = s.get_optional("ForcedReviewError", key)
        turn = s.get_meta("finished_turn", 0) + 1
        if turn > ctrl.budget["max_agent_turns"]:
            raise StopReview("routed_turn_budget")
        cap = ctrl.budget.get("max_reported_tokens")
        if cap:
            used = [m["usage"].get("total_tokens") for m in s.all("Message")]
            if None in used or sum(used) >= cap:
                raise StopReview("max_reported_tokens" if None not in used else "token_usage_unreported")
        sent = dict(context, previous_error=last["error"]) if last else context
        msg = obtain(ctrl, "live", turn, retry_unknown, context=sent)
        mid = msg["message_id"]
        try:
            reply = response(msg["raw"])
            if not msg.get("complete", True) or not isinstance(reply, dict):
                raise ValueError("Expected one complete review object")
            validate(reply)
        except (ValueError, KeyError, TypeError, AttributeError) as exc:
            error = str(exc)[:200]
            with s.transaction():
                # A new turn can retry without overwriting an immutable response.
                if not last:
                    s.put("ForcedReviewError", key, {"error": error})
                route.complete_turn(mid, {"error": error})
            raise StopReview("routed_action_blocked" if last else "retry_forced_review")
        with s.transaction():
            s.put("ForcedReview", key, dict(context_ref=s.artifact(context), reply=reply, message_id=mid))
            route.complete_turn(mid, dict(action=context["phase"]))
        if pause_after is not None and turn >= pause_after:
            raise StopReview("paused")
        return reply
    claims = [dict(r, members=ctrl.candidate(r["candidate_id"])["members"]) for r in s.all("RouteFreeze")
            if r.get("candidate_id") and r.get("kind") == "lineage"]
    while True:
        try:
            result = adjudicate(ctrl, claims, ask)
            with s.transaction():
                project(ctrl, result)
            return "complete", "submitted"
        except StopReview as exc:
            if str(exc) == "retry_forced_review":
                continue
            s.set_meta("current_issue", None if str(exc) == "paused" else str(exc))
            return "partial", str(exc)
        except providers.ProviderError as exc:
            s.set_meta("current_issue", str(exc))
            return "partial", "provider_blocked" if exc.status == "blocked" else exc.status
