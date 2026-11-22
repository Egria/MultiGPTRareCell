"""Evidence-bound membership decisions for the opt-in adaptive discovery arm."""
from .config import fingerprint

GOAL = (
    " Favor discovery of plausible rare memberships, including unnamed provisional groups. Investigate promising "
    "singleton/outlier seeds across coarse groups early. If a computation was blocked, correct it or record why "
    "you defer the seed. Inspect ranks, gaps and co-detection against a stated reference. For a singleton, compare "
    "positive-only and combined expression before treating zeros as exclusion evidence. You own the boundary: "
    "choose a measured cutoff per candidate, citing core positive support, outside contrast and uncertainty. "
    "Use global and parent-relative <=5%, <=1%, <=0.1% as nested search goals, without quotas. Record the "
    "denominator. Propose biological rarity hypotheses in tissue/species context and test observable programs. "
    "Retain a coherent but uncertain split as provisional when the cost of missing it outweighs a false discovery. "
    "State an alternative explanation and what observation would weaken the split. Low abundance alone is insufficient. "
    "Keep prior knowledge, abundance and boundary support distinct. Reusing the same RNA is selection-dependent. "
    "Preview competing claims before submitting; a broad claim can absorb a rare hypothesis. Explain deferred "
    "promising candidates in your submission reason. Keep the discovery and token budgets."
)
TOOLS = {
    "apply_boundary": "candidate_id (original seed), evidence_ids (one expansion observation), metric "
        "(score|expression|positive_expression), threshold [0,1], reason. Return an anchored provisional candidate "
        "using all scored cells at or above your cutoff, with no hidden expression or structure cutoff. "
        "max_cells bounds additions; an over-budget rule reports its exact size without truncation. "
        "Read the expansion ranking with get_evidence. A new cutoff reuses the same observation.",
    "preview_submission": "candidate_ids (proposed selection; omit for current selection), optional candidate_id "
        "(hypothesis to track, even when omitted from selection), optional reference_id (reference candidate). "
        "Reports exclusive claims, conflicts, unclaimed cells, and local/global fractions before final assignment.",
}


def abundance(count, total):
    bands = [pct for pct in (0.1, 1, 5) if count and 100 * count <= pct * total]
    return {"n": count, "denominator": total, "fraction": count / total if total else None,
            "smallest_band_pct": min(bands) if bands else None}


def coexpression_members(ctrl, eid):
    ctrl.checked_evidence([eid])
    key = ctrl.state.get("Evidence", eid)["object"]
    try:
        observation = ctrl.state.get("QueryObservation", key)
        query = ctrl.state.get("CanonicalQuery", key)
    except KeyError:
        raise ValueError("matched_evidence_id requires a measured coexpress observation") from None
    if observation["evidence_id"] != eid or query["canonical"]["action"] != "coexpress":
        raise ValueError("matched_evidence_id requires a measured coexpress observation")
    return ctrl.state.read_artifact(ctrl.state.result_ref(key))["matched_members"]


def apply_boundary(ctrl, args):
    seed = ctrl.candidate(args["candidate_id"])
    eid, = ctrl.checked_evidence(args["evidence_ids"])
    key = ctrl.state.get("Evidence", eid)["object"]
    # An observation about a different seed (even a similarly sized seed)
    # cannot authorize a new membership for this lineage.
    if not any(row["query_key"] == key and row["candidate_id"] == seed["id"] and row["evidence_id"] == eid
               for row in ctrl.state.all("Investigation")):
        raise ValueError("Use the original candidate_id and its expand_candidate evidence_id")
    ranked = ctrl.state.read_artifact(ctrl.state.result_ref(key)).get("ranking", [])
    metric, threshold = args["metric"], args["threshold"]
    values = [row for row in ranked if row[metric] is not None]
    if not values:
        raise ValueError(f"No scored cells for {metric} in this expansion")
    anchor = set(seed["members"])
    members = anchor.union(row["cell_id"] for row in values if row[metric] >= threshold)
    added = len(members) - len(anchor)
    response = {"action": "apply_boundary", "evidence_id": eid, "anchor_id": seed["id"], "metric": metric,
                "threshold": threshold, "added_count": added, "global_abundance": abundance(len(members), len(ctrl.cells)),
                "membership_hash": fingerprint(sorted(members))}
    if added > ctrl.cfg["search"]["max_cells"]:
        return dict(response, status="over_limit", max_additions=ctrl.cfg["search"]["max_cells"])
    if not added:
        return dict(response, status="no_expansion", candidate_id=seed["id"])
    rule = {"metric": metric, "threshold": threshold, "evidence_id": eid}
    revision = ctrl.make_revision(members, [seed["id"]], "unknown", "provisional", [eid], "boundary", rule)
    return dict(response, status="proposed", candidate_id=revision["id"])


def preview(ctrl, args):
    from .finalize import membership_claims
    selected = args.get("candidate_ids", ctrl.state.get_meta("selected", []))
    claims = membership_claims(ctrl, selected)
    exclusive = {cell for cell, owners in claims.items() if len(owners) == 1}
    focus_ids = [args["candidate_id"]] if "candidate_id" in args else selected
    reference = set(ctrl.candidate(args["reference_id"])["members"]) if "reference_id" in args else None
    report = []
    for cid in focus_ids:
        members = set(ctrl.candidate(cid)["members"])
        unique = members & exclusive
        retained = {cell for cell in unique if cid in claims[cell]}
        row = {"candidate_id": cid, "selected": cid in selected, "retained": len(retained),
               "assigned_to_other": len(unique) - len(retained),
               "conflicts": sum(len(claims.get(cell, ())) > 1 for cell in members),
               "unclaimed": sum(cell not in claims for cell in members),
               "proposed_global": abundance(len(members), len(ctrl.cells)),
               "retained_global": abundance(len(retained), len(ctrl.cells))}
        if reference is not None:
            row.update(reference_id=args["reference_id"], proposed_local=abundance(len(members & reference), len(reference)),
                       retained_local=abundance(len(retained & reference), len(reference)))
        report.append(row)
    return {"action": "preview_submission", "exclusive_cells": len(exclusive),
            "conflicted_cells": len(claims) - len(exclusive), "unclaimed_cells": len(ctrl.cell_set)-len(claims),
            "candidates_reviewed": report}
