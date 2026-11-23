"""Turn written singleton follow-ups into bounded measurements of exact-CL pools."""
from .config import contract, json_bytes
from .protocol import ProtocolError, resolve_genes
from .state import BudgetExhausted, QueryBlocked


def queue(ctrl):
    notes = {r["candidate_id"]: r for r in ctrl.state.all("SeedAssessment") if r["disposition"] == "follow_up"}
    return [(c, [notes[p] for p in c["parents"] if p in notes]) for c in ctrl.pool_candidates
            if any(p in notes for p in c["parents"])]


def coverage(ctrl):
    planned = queue(ctrl)
    rows = ctrl.state.all("PoolFollowup")
    covered = {r["candidate_id"] for _, notes in planned for r in notes}
    total = sum(r["disposition"] == "follow_up" for r in ctrl.state.all("SeedAssessment"))
    return {"eligible_pools": len(planned), "measured": sum(r["status"] == "ok" for r in rows),
            "blocked": sum(r["status"] != "ok" for r in rows), "pending": len(planned) - len(rows),
            "pooled_follow_up_seeds": len(covered), "unpooled_follow_up_seeds": total - len(covered),
            "handoff_turns": ctrl.state.get_meta("pool_review_turns", 0)}


def prepare_batch(ctrl, mid):
    s = ctrl.state
    if s.get_optional("PoolBatch", mid):
        return True
    if s.get_meta("pool_phase_complete", False):
        return False
    completed = {r["candidate_id"] for r in s.all("PoolFollowup")}
    remaining = ctrl.query_limit - s.usage()["new_queries"]
    size = min(remaining, ctrl.budget["max_actions_per_turn"], ctrl.budget.get("seed_review_batch_size", 4))
    todo = [(c, notes) for c, notes in queue(ctrl) if c["id"] not in completed]
    if not todo or size <= 0:
        s.set_meta("pool_phase_complete", True)
        return False
    ctx = {"phase": "pool_follow_up", "mode": ctrl.mode, "prompt_version": ctrl.prompt, "contract": contract(),
           "species": ctrl.manifest.get("species"), "tissue": ctrl.manifest.get("tissue"),
           "population_size": len(ctrl.cells), "pools": [],
           "actions": {"investigate_pool": "candidate_id, genes (1 to 6), reason; optional min_detected."},
           "turn_instruction": "Return one investigate_pool per listed pool. Use its written questions and marker "
           "panels to choose a discriminating test, including held-out genes when appropriate. Default: compare "
           "this entire seed pool to all other cells on the named genes, including depth-adjusted NB splits when "
           "raw counts exist. For dataset-wide co-detection instead, specify min_detected (how many panel genes "
           "must be detected). Counts and exact memberships will be available in discovery. Labels are hypotheses."}
    if not ctrl.adaptive:
        ctx["turn_instruction"] += " For co-detection in this mode choose one required gene or require all panel genes."
    for c, notes in todo[:size]:
        item = ctrl.pool_card(c)
        item["marker_seed_counts"] = dict(c["ontology_pool"]["marker_counts"][:48])
        item["n_follow_up"] = len(notes)
        item["written_questions"] = [{"seed": r["candidate_id"], "reason": r["reason"][:480]} for r in notes[:4]]
        ctx["pools"].append(item)
        if len(json_bytes(ctx)) > ctrl.budget["context_bytes"]:
            ctx["pools"].pop()
            if not ctx["pools"]:
                raise ValueError("context_bytes cannot fit one ontology pool follow-up")
            break
    s.put("PoolBatch", mid, {"candidate_ids": [c["candidate_id"] for c in ctx["pools"]],
                             "context_ref": s.artifact(ctx)})
    return True


def run_query(ctrl, mid, args, *, fallback=False):
    s = ctrl.state
    cid = args["candidate_id"]
    batch = s.get_optional("PoolBatch", mid)
    if not batch or cid not in batch["candidate_ids"]:
        raise ProtocolError("investigate_pool requires a pool assigned in this message")
    saved = s.get_optional("PoolFollowup", cid)
    if saved:
        return saved
    c = ctrl.candidate(cid)
    if not 1 <= len(args["genes"]) <= 6:
        raise ProtocolError("choose one to six genes for this pool")
    genes = resolve_genes(args["genes"], ctrl.genes, species=ctrl.manifest.get("species"))
    k = args.get("min_detected")
    params = {"genes": genes, "evidence_ids": c["evidence_ids"]}
    if k is None:
        action = "compare_groups"
        params["candidate_id"] = cid
    else:
        action = "coexpress"
        params["scope"] = "all"
        if k == len(genes) or k == 1:
            params["logic"] = "AND" if k == len(genes) else "OR"
        elif ctrl.adaptive and 1 < k < len(genes):
            params["min_detected"] = k
        else:
            raise ProtocolError("min_detected is outside the supported range for this mode and panel")
    row = {"candidate_id": cid, "label": c["proposed_label"], "n_seeds": len(c["members"]),
           "source": "controller_fallback" if fallback else "discovery_agent", "reason": args["reason"],
           "message_id": mid, "action": action, "genes": genes}
    if k is not None:
        row["min_detected"] = k
    try:
        obs = ctrl.scientific(action, params)
    except (QueryBlocked, BudgetExhausted) as exc:
        row.update(status="blocked", error=str(exc)[:300])
    else:
        data = s.read_artifact(obs["artifact_ref"])
        row.update(status=obs["status"], evidence_id=obs["evidence_id"], query_key=obs["query_key"])
        if k is not None:
            hits = data["matched_members"]
            row.update(n_matching=len(hits), matching_seeds=len(set(hits).intersection(c["members"])))
        else:
            row.update(n_background=data["summary"]["n_background"], contrasts={})
            for stat in data.get("statistics", []):
                gene = stat["gene_id"]
                entry = {key: round(stat[key], 4) for key in ("effect", "fraction_query", "fraction_background")}
                diag = data.get("nb_diagnostics", {}).get(gene, {})
                entry["nb"] = diag.get("status", data["summary"].get("nb_mixture", {}).get("status", "unavailable"))
                if "high_foreground" in diag:
                    entry.update(high_seeds=diag["high_foreground"], high_outside=diag["high_background"],
                                 means_per_10k=[round(v, 4) for v in diag["means_per_10k"]])
                row["contrasts"][gene] = entry
    s.put("PoolFollowup", cid, row)
    return row


def finish_batch(ctrl, mid, batch):
    """A missing or malformed nomination still gets a traceable numerical investigation."""
    s = ctrl.state
    for cid in batch["candidate_ids"]:
        if s.get_optional("PoolFollowup", cid):
            continue
        c = ctrl.candidate(cid)
        panel = set()
        for symbol, _ in c["ontology_pool"]["marker_counts"]:
            try:
                panel.update(resolve_genes([symbol], ctrl.genes, species=ctrl.manifest.get("species")))
            except ProtocolError:
                continue
            if len(panel) >= 4:
                break
        if panel:
            run_query(ctrl, mid, {"candidate_id": cid, "genes": sorted(panel), "reason":
                      "No usable nomination. Measure recurrent seed-panel genes against the complement."}, fallback=True)
        else:
            s.put("PoolFollowup", cid, {"candidate_id": cid, "label": c["proposed_label"], "status": "blocked",
                  "message_id": mid, "source": "controller_fallback", "error": "No unambiguous pool marker resolves."})


def card(row):
    return {k: v for k, v in row.items() if k not in {"query_key", "message_id", "reason"}}
