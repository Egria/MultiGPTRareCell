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
           "actions": {"investigate_pool": "candidate_id, reason, seed_genes, optional support_genes and counter_genes (up to six per role)."},
           "turn_instruction": "Choose a discriminating seed_genes program from this pool's observed marker union, "
           "including specific markers seen on different seeds. Name additional positive markers in support_genes "
           "and competing-lineage genes in counter_genes. One local query measures contrast and all k-of-m joint "
           "detection boundaries separately for each program. Inspect pool and outside hits, depth-matched background "
           "when raw counts exist, and competing explanations. Existing genes/min_detected requests are also accepted."}
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
    panels = {}
    k = args.get("min_detected")
    role_fields = ("seed_genes", "support_genes", "counter_genes")
    nominated = any(args.get(field) for field in role_fields)
    if any(len(args.get(field, [])) > 6 for field in (*role_fields, "genes")):
        raise ProtocolError("choose up to six genes per role")
    def resolve(tokens):
        return resolve_genes(tokens, ctrl.genes, species=ctrl.manifest.get("species")) if tokens else []
    base = resolve(args.get("genes", []))
    # A signed request already measures every cutoff, including a supplied hint.
    if nominated:
        k = None
    if k is None:
        observed = dict(c["ontology_pool"]["marker_counts"])
        for field in role_fields:
            tokens = args.get(field, [])
            if not tokens:
                continue
            panels[field[:-6]] = resolve(tokens)
        if "seed" in panels:
            known = set(observed)
            genes_by_id = {r["gene_id"]: r.get("symbol") for r in ctrl.genes}
            if any(g not in known and genes_by_id.get(g) not in known for g in panels["seed"]):
                raise ProtocolError("seed_genes must come from this pool's marker union")
        if not panels.get("seed"):
            seed = []
            for symbol in observed:
                try:
                    seed.extend(resolve([symbol]))
                except ProtocolError:
                    continue
                if len(set(seed)) >= 6:
                    break
            if seed:
                panels["seed"] = sorted(set(seed))
        counter = set(panels.get("counter", []))
        if counter.intersection(panels.get("support", [])) or counter.intersection(panels.get("seed", [])):
            if args.get("seed_genes") or counter.intersection(panels.get("support", [])):
                raise ProtocolError("Counter genes cannot also nominate positive support")
            panels["seed"] = [g for g in panels["seed"] if g not in counter]
        if base and "support" not in panels:
            panels["inferred_support"] = sorted(set(base) - counter)
        genes = sorted(set(base).union(*panels.values()))
        if not genes:
            raise ProtocolError("choose at least one measurable gene")
        action = "compare_groups"
        params = {"candidate_id": cid, "genes": genes, "pool_programs": panels, "evidence_ids": c["evidence_ids"]}
    else:
        genes = base
        params = {"genes": genes, "evidence_ids": c["evidence_ids"]}
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
    if "min_detected" in args:
        row["min_detected"] = args["min_detected"]
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
            row.update(n_background=data["summary"]["n_background"], contrasts={},
                       programs=data["summary"].get("programs", {}), depth_matching=data["summary"].get("depth_matching"),
                       seed_panel_source="agent" if args.get("seed_genes") else "pool_frequency_fallback")
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


def card(row, compact=False):
    item = {k: v for k, v in row.items() if k not in {"query_key", "message_id", "reason"}}
    if compact and row.get("programs"):
        # Keep the full counts in get_evidence, with a small all-cutoffs index in each turn.
        item["programs"] = {name: {"genes": p["genes"], "seed_hits_by_k": [x["seed_hits"] for x in p["curve"]],
                                   "outside_hits_by_k": [x["outside_hits"] for x in p["curve"]]}
                            for name, p in row["programs"].items()}
        item["contrasts"] = {g: {k: v for k, v in stat.items() if k in ("effect", "high_seeds", "high_outside")}
                             for g, stat in row["contrasts"].items()}
    return item
