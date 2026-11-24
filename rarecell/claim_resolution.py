"""Plan adaptive assignments once, with a checkpointed fallback for competing drafts."""
from collections import Counter, defaultdict

from .finalize import membership_claims


def _known_label(ctrl, cid):
    c = ctrl.candidate(cid)
    return c["proposed_label"] if c["proposed_label"] != "unknown" else ctrl.current_label(cid).get("label", "unknown")


def plan(ctrl, requested, preferred=()):
    ids = sorted(set(requested))
    preferred = list(dict.fromkeys(preferred))
    if not set(preferred) <= set(ids):
        raise ValueError("preferred_candidate_ids must be included in the proposed selection")
    ref = ctrl.state.get_meta("claim_plan_ref")
    prior = ctrl.state.read_artifact(ref) if ref else {}
    if prior.get("requested") == ids and (not preferred or preferred == prior["preferred"]):
        return prior
    if not preferred:
        preferred = [cid for cid in prior.get("preferred", []) if cid in ids]
    prior_owners = prior.get("owners", {})
    prior_provisional = set(prior.get("provisional_cells", []))
    prior_order = prior.get("order", ctrl.state.get_meta("selected", []))
    old = {cid: i for i, cid in enumerate(prior_order)}
    pref = {cid: i for i, cid in enumerate(preferred)}
    sets = {cid: set(ctrl.candidate(cid)["members"]) for cid in ids}
    labels = {cid: _known_label(ctrl, cid) for cid in ids}
    prior_retained = set(prior_owners.values())
    order = sorted(ids, key=lambda cid: (pref.get(cid, len(pref)), cid not in prior_retained,
                                         old.get(cid, len(old)), len(sets[cid]), cid))
    held, keep, alternatives = {}, [], defaultdict(list)
    majority_rivals = defaultdict(set)
    for cid in order:
        c = ctrl.candidate(cid)
        for rep in keep:
            other = ctrl.candidate(rep)
            overlap = len(sets[cid] & sets[rep])
            if not overlap:
                continue
            exact = sets[cid] == sets[rep]
            # Compare each draft to an actual surviving representative, never transitively.
            shared_anchor = (c.get("boundary_rule") and other.get("boundary_rule") and c["parents"] == other["parents"])
            same_label = labels[cid].casefold() == labels[rep].casefold() and labels[cid].casefold() != "unknown"
            nested = overlap == min(len(sets[cid]), len(sets[rep]))
            mutual_majority = 2 * overlap > max(len(sets[cid]), len(sets[rep]))
            if exact or shared_anchor or (same_label and (nested or mutual_majority)):
                held[cid] = {"representative": rep, "intersection": overlap, "draft_n": len(sets[cid]),
                             "representative_n": len(sets[rep]), "rule": "identical" if exact else "alternative_boundary"}
                if exact and labels[cid] != "unknown" and labels[cid] != labels[rep]:
                    alternatives[rep].append(labels[cid])
                break
            if mutual_majority:
                majority_rivals[cid].add(rep)
                majority_rivals[rep].add(cid)
        else:
            keep.append(cid)
    claims = membership_claims(ctrl, keep)
    exclusive = Counter(next(iter(v)) for v in claims.values() if len(v) == 1)
    rank = {cid: i for i, cid in enumerate(keep)}
    owners, unresolved, rescued, provisional = {}, {}, Counter(), []
    for cell, rivals in claims.items():
        if len(rivals) == 1:
            owners[cell] = next(iter(rivals))
            continue
        # An explicit preference or previously viable assignment is already evidence of
        # the agent's boundary choice. Retain exclusive flanks of other hypotheses.
        choices = rivals.intersection(pref)
        if choices:
            owner = min(choices, key=pref.get)
        elif prior_owners.get(cell) in rivals:
            owner = prior_owners[cell]
            if cell in prior_provisional:
                provisional.append(cell)
        else:
            at_risk = [cid for cid in rivals if not exclusive[cid] or majority_rivals[cid].intersection(rivals)]
            if not at_risk:
                unresolved[cell] = sorted(rivals)
                continue
            owner = min(at_risk, key=rank.get)
            provisional.append(cell)
        owners[cell] = owner
        rescued[owner] += 1
    return {"requested": ids, "preferred": preferred, "order": keep, "held": held,
            "alternative_labels": dict(alternatives), "owners": owners, "unresolved": unresolved,
            "resolved_overlaps": dict(rescued), "provisional_cells": sorted(provisional)}


def summary(p):
    retained = Counter(p["owners"].values())
    return {"retained": {cid: retained[cid] for cid in p["requested"]}, "held": p["held"],
            "alternative_labels": p["alternative_labels"], "resolved_overlaps": p["resolved_overlaps"],
            "unresolved_overlap_cells": len(p["unresolved"]), "provisional_overlap_cells": len(p["provisional_cells"])}


def commit(ctrl, requested, preferred=()):
    result = plan(ctrl, requested, preferred)
    ctrl.state.set_meta("claim_plan_ref", ctrl.state.artifact(result))
    return result
