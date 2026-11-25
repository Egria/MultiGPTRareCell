"""Evidence-bound membership decisions for the opt-in adaptive discovery arm."""
from .config import fingerprint
from .protocol import resolve_genes

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
    "promising candidates in your submission reason. Keep the discovery and token budgets. "
    "Treat the core and wider membership as competing drafts. Assess the incremental cells' positive program, "
    "their fit to the core versus the outside and the impact on other groups. The wider gap is a search suggestion. "
    "Choose either set, another measured boundary, or neither. Submit once this comparison is sufficient."
)
TOOLS = {
    "apply_boundary": "candidate_id (original seed), evidence_ids (one expansion observation), metric "
        "(score|expression|positive_expression), threshold [0,1], reason. Return an anchored provisional candidate "
        "using all scored cells at or above your cutoff, with no hidden expression or structure cutoff. "
        "max_cells bounds additions; an over-budget rule reports its exact size without truncation. "
        "The first nonconflicting boundary you request per anchor enters the working selection. A wider alternative "
        "is held separately. Subsequent cutoffs are held for comparison. See boundary_reviews for current choices.",
    "select_boundary": "candidate_id from a held boundary review and reason. Replace the working boundary of its "
        "anchor and preserve other disjoint selections. Use propose_revision to resolve cross-anchor overlaps.",
    "preview_submission": "candidate_ids (proposed selection; omit for current selection), optional candidate_id "
        "(hypothesis to track, even when omitted from selection), optional reference_id (reference candidate). "
        "Optional ordered preferred_candidate_ids. Reports resolved, held and unresolved overlapping claims before export.",
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


def unclaimed_coexpression(ctrl):
    """Pair a saved group's core with exact, previously measured outside hits."""
    observations = ctrl.state.all("CoexpressionHits")
    if not observations:
        return []
    selected = set(ctrl.state.get_meta("selected", []))
    occupied = {c for cid in selected for c in ctrl.candidate(cid)["members"]}
    candidates = set(selected)
    for review in ctrl.state.all("BoundaryReview"):
        candidates.update(review[k]["candidate_id"] for k in ("requested", "wider") if k in review)
    rows = []
    for obs in observations:
        hits = set(obs["members"])
        unclaimed = hits - occupied
        if not unclaimed:
            continue
        for cid in sorted(candidates):
            core = set(ctrl.candidate(cid)["members"])
            overlap = len(hits & core)
            extra = sorted(unclaimed - core)
            if overlap and extra:
                rows.append(dict(candidate_id=cid, evidence_id=obs["evidence_id"], genes=obs["genes"],
                                 min_detected=obs["min_detected"], threshold=obs["threshold"],
                                 core_hits=overlap, n_hits=len(hits), n_unclaimed=len(extra),
                                 other_claimed_hits=len(hits & (occupied - core)), examples=extra[:8]))
    return sorted(rows, key=lambda v:(v["n_unclaimed"], v["other_claimed_hits"], -v["core_hits"],
                                       v["evidence_id"], v["candidate_id"]))


def program_members(ctrl, eid, name, k):
    ctrl.checked_evidence([eid])
    key = ctrl.state.get("Evidence", eid)["object"]
    obs = ctrl.state.get_optional("QueryObservation", key)
    if not obs or obs["evidence_id"] != eid:
        raise ValueError("A program boundary requires a measured pool observation")
    result = ctrl.state.read_artifact(ctrl.state.result_ref(key))
    program = result.get("programs", {}).get(name)
    if not program or not 1 <= k <= len(program["genes"]):
        raise ValueError("Unknown program or cutoff for this pool observation")
    return sorted(cell for count, cells in program["members_by_count"].items() if int(count) >= k for cell in cells)


def nb_high_members(ctrl, eid, gene):
    ctrl.checked_evidence([eid])
    key = ctrl.state.get("Evidence", eid)["object"]
    obs = ctrl.state.get_optional("QueryObservation", key)
    if not obs or obs["evidence_id"] != eid or ctrl.state.get("CanonicalQuery", key)["canonical"]["action"] != "compare_groups":
        raise ValueError("nb_high_gene requires a measured contrast")
    gene, = resolve_genes([gene], ctrl.genes, species=ctrl.manifest.get("species"))
    result = ctrl.state.read_artifact(ctrl.state.result_ref(key))
    fit = result["summary"].get("nb_mixture", {})
    split = fit.get("genes", {}).get(gene, {})
    if "minor_members" not in split:
        raise ValueError("No fitted two-component split for " + gene)
    if split["minor_mean_per_10k"] > split["major_mean_per_10k"]:
        return split["minor_members"]
    scope = result["selection_history"]
    return sorted((set(scope["members"]) | set(scope["background"])) - set(fit.get("excluded_zero_depth", []))
                  - set(split["minor_members"]))


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
    def option(cell_ids, cut):
        rule = {"metric": metric, "threshold": cut, "evidence_id": eid}
        candidate = ctrl.make_revision(cell_ids, [seed["id"]], "unknown", "provisional", [eid], "boundary", rule)
        return {"candidate_id": candidate["id"], "threshold": cut, "n": len(cell_ids)}

    requested = option(members, threshold)
    cid = requested["candidate_id"]
    reviews = ctrl.state.all("BoundaryReview")
    review = next((r for r in reviews if r["requested"]["candidate_id"] == cid), None)
    if review is None:
        # Only non-anchor members can change. Enumerate complete tie blocks and
        # compare realized group sizes to the existing additions budget.
        others = sorted((r for r in values if r["cell_id"] not in anchor), key=lambda r: -r[metric])
        best_gap, split = 0.0, None
        for n in range(added + 1, min(len(others) - 1, ctrl.cfg["search"]["max_cells"]) + 1):
            gap = others[n-1][metric] - others[n][metric]
            if gap > max(best_gap, 1e-12):
                best_gap, split = gap, n
        review = {"anchor_id": seed["id"], "evidence_id": eid, "metric": metric, "requested": requested}
        wider = members
        if split is not None:
            wider = anchor.union(r["cell_id"] for r in others[:split])
            review["wider"] = option(wider, others[split-1][metric])
            review["gap_below_wider"] = best_gap
        # Include the nearest excluded cells separately so that a large majority
        # cannot dilute the contrast with the most competitive outside cells.
        nearest = [r for r in others if r["cell_id"] not in wider][:max(1, len(wider)-len(members))]
        # Aggregate the same observed expression, including cells scored on a
        # different component, without treating this comparison as new evidence.
        cohorts = {"core": [], "additional": [], "outside": []}
        for row in ranked:
            cell = row["cell_id"]
            cohort = "core" if cell in members else "additional" if cell in wider else "outside"
            cohorts[cohort].append(row)
        cohorts["nearest_outside"] = nearest
        review["cohorts"] = {}
        for name, rows in cohorts.items():
            metrics = {}
            for field in ("score", "positive_expression", "expression", "structure"):
                observed = [r[field] for r in rows if r[field] is not None]
                if observed:
                    metrics[field] = round(sum(observed) / len(observed), 5)
            review["cohorts"][name] = {"scored": len(rows), "mean": metrics}
        ctrl.state.put("BoundaryReview", cid, review)
        first = not any(r["anchor_id"] == seed["id"] for r in reviews)
        selected = ctrl.state.get_meta("selected", [])
        if first and not any(members.intersection(ctrl.candidate(i)["members"]) for i in selected):
            from .claim_resolution import commit
            commit(ctrl, [*selected, cid])
            ctrl.state.set_meta("selected", sorted([*selected, cid]))
    return dict(response, status="proposed" if added else "no_expansion", candidate_id=cid,
                boundary_review=review_card(ctrl, review))


def other_selections(ctrl, anchor_id):
    # The competing drafts all have this exact anchor as their sole parent.
    return [cid for cid in ctrl.state.get_meta("selected", [])
            if not (ctrl.candidate(cid).get("boundary_rule") and ctrl.candidate(cid)["parents"] == [anchor_id])]


def review_card(ctrl, review):
    from .claim_resolution import plan
    other = other_selections(ctrl, review["anchor_id"])
    selected = ctrl.state.get_meta("selected", [])
    result = dict(review)
    for role in ("requested", "wider"):
        if role not in review:
            continue
        item = review[role]
        cid = item["candidate_id"]
        # select_boundary explicitly prefers the requested option over competing drafts.
        assignment = plan(ctrl, [*other, cid], [cid])
        result[role] = dict(item, selected=cid in selected,
                            exclusive_members=sum(owner == cid for owner in assignment["owners"].values()),
                            partition_conflicts=len(assignment["unresolved"]),
                            partition_unclaimed=len(ctrl.cells)-len(assignment["owners"])-len(assignment["unresolved"]))
    return result


def select_boundary(ctrl, args):
    cid = args["candidate_id"]
    c = ctrl.candidate(cid)
    if not c.get("boundary_rule"):
        raise ValueError("Choose a candidate_id from a boundary review")
    other = other_selections(ctrl, c["parents"][0])
    from .claim_resolution import commit
    chosen = sorted([*other, cid])
    commit(ctrl, chosen, [cid])
    ctrl.state.set_meta("selected", chosen)
    return {"action": "select_boundary", "candidate_id": cid, "selected_ids": chosen, "reason": args["reason"]}


def preview(ctrl, args):
    from .claim_resolution import plan, summary
    selected = args.get("candidate_ids", ctrl.state.get_meta("selected", []))
    resolution = plan(ctrl, selected, args.get("preferred_candidate_ids", ()))
    owners = resolution["owners"]
    conflicts = resolution["unresolved"]
    exclusive = set(owners)
    focus_ids = [args["candidate_id"]] if "candidate_id" in args else selected
    reference = set(ctrl.candidate(args["reference_id"])["members"]) if "reference_id" in args else None
    report = []
    for cid in focus_ids:
        members = set(ctrl.candidate(cid)["members"])
        unique = members & exclusive
        retained = {cell for cell in unique if owners[cell] == cid}
        row = {"candidate_id": cid, "selected": cid in selected, "retained": len(retained),
               "assigned_to_other": len(unique) - len(retained),
               "conflicts": sum(cell in conflicts for cell in members),
               "unclaimed": sum(cell not in owners and cell not in conflicts for cell in members),
               "proposed_global": abundance(len(members), len(ctrl.cells)),
               "retained_global": abundance(len(retained), len(ctrl.cells))}
        if reference is not None:
            row.update(reference_id=args["reference_id"], proposed_local=abundance(len(members & reference), len(reference)),
                       retained_local=abundance(len(retained & reference), len(reference)))
        report.append(row)
    return {"action": "preview_submission", "exclusive_cells": len(exclusive),
            "claim_resolution": summary(resolution), "conflicted_cells": len(conflicts),
            "unclaimed_cells": len(ctrl.cell_set) - len(exclusive) - len(conflicts),
            "candidates_reviewed": report}
