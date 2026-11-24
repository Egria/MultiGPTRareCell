"""Pack model-requested observations and the measured handoff index before optional context."""
from .config import json_bytes
from .pool_followup import card


def index(row):
    brief = {k: row[k] for k in ("candidate_id", "label", "evidence_id", "status", "n_seeds", "n_background",
                                 "anchor_id", "kind", "n_foreground", "expansion_evidence_id", "matching_seeds", "n_matching") if k in row}
    if row.get("written_questions") and row.get("annotations"):
        brief["follow_up_reason"] = row["written_questions"][0]["reason"][:350]
        brief["global_alternative"] = row["annotations"]["global"]["label"]
    # Put one directly addressable enriched boundary in the directory. The complete
    # per-gene fits and all program cutoffs remain available by candidate ID.
    fits = [(g, f) for g, f in row.get("contrasts", {}).items() if f.get("high_seeds") and f.get("effect", 0) > 0]
    if fits:
        gene, f = max(fits, key=lambda item: item[1]["high_seeds"] /
                      (row["n_seeds"] + item[1]["high_seeds"] + item[1]["high_outside"]))
        brief["nb_high_gene"] = gene
        brief.update({k: f[k] for k in ("high_seeds", "high_outside", "means_per_10k")})
    return brief


def brief_receipt(result):
    if "pool_followup" in result:
        pool = result["pool_followup"]
        return {"candidate_id": pool["candidate_id"], "pool_followup": card(pool, compact=True),
                "offset": result["offset"], "total": result["total"], "members": result["members"]}
    brief = {k: result[k] for k in ("action", "candidate_id", "evidence_id", "query_key", "offset", "limit",
                                     "status", "error", "total", "members", "nb_high_counts") if k in result}
    if "candidate" in result:
        brief["candidate_id"] = result["candidate"]["id"]
    if "summary" in result:
        summary = dict(result["summary"])
        mixture = summary.pop("nb_mixture", None)
        if mixture:
            summary["nb_genes"] = {g: {**{k: v for k, v in fit.items() if k != "minor_members"},
                                       "minor_n": len(fit.get("minor_members", []))}
                                   for g, fit in mixture.get("genes", {}).items()}
        brief["summary"] = summary
    if "sections" in result:
        brief["sections"] = {k: {"total": v["total"], "rows": v["rows"][:2]} for k, v in result["sections"].items()}
    if "evidence" in result:
        brief["evidence"] = [{k: e[k] for k in ("evidence_id", "kind", "object")} for e in result["evidence"]]
    return brief


def requested(ctx, state, limit):
    pending = state.get_meta("requested_receipts", [])
    if not pending:
        return set()
    ctx["requested_evidence"] = []
    ctx["pending_evidence_count"] = len(pending)
    entries = ctx["requested_evidence"]
    space = limit - len(json_bytes(ctx))
    detailed = []
    delivered = set()
    for i, rid in enumerate(pending):
        rec = state.get("Receipt", rid)["result"]
        short = brief_receipt(rec)
        variants = [rec, short]
        if rec.get("pool_followup"):
            variants.append(index(rec["pool_followup"]))
        else:
            # Large member lists can exceed even an otherwise empty context. Preserve
            # the binding and useful numbers, and allow the model to request a small page.
            tiny = {k: short[k] for k in ("action", "candidate_id", "evidence_id", "query_key", "offset", "total", "status")
                    if k in short}
            tiny["abbreviated"] = "Full measurement archived. Read a smaller get_evidence page if needed."
            summary = short.get("summary", {})
            tiny.update({k: summary[k] for k in ("n_query", "n_background", "n_matching", "n_added", "n_passing")
                         if k in summary})
            if short.get("nb_high_counts"):
                gene, counts = max(short["nb_high_counts"].items(), key=lambda v: v[1]["high_foreground"] /
                                   max(1, summary.get("n_query", 0) + sum(v[1].values())))
                tiny.update(nb_high_gene=gene, **counts)
            for name, section in short.get("sections", {}).items():
                tiny.setdefault("sections", {})[name] = {"total": section["total"], "rows": section["rows"][:1]}
            if short.get("evidence"):
                tiny["evidence"] = short["evidence"][:1]
            variants.append(tiny)
        # Fair-share within the turn. If even the compact result cannot fit,
        # keep it pending for the next dispatch and show the remaining count.
        share = space // (len(pending) - i)
        for payload in variants:
            entry = {"receipt_id": rid, "result": payload}
            cost = len(json_bytes(entry)) + 1
            if cost <= share or (payload is variants[-1] and cost <= space):
                entries.append(entry)
                delivered.add(rid)
                space -= cost
                if payload is not rec:
                    detailed.append((entry, rec))
                break
    # Use spare space for the original page, including precisely the requested rows.
    for entry, full in detailed:
        extra = len(json_bytes(full)) - len(json_bytes(entry["result"]))
        if extra <= space:
            entry["result"] = full
            space -= extra
    ctx["pending_evidence_count"] -= len(delivered)
    return delivered


def directory(ctx, state, limit):
    rows = state.all("PoolFollowup")
    if not rows:
        return
    # Give every measured case an entry before spending the remainder on larger cards.
    # Rows with an alternative lineage absent from the global pools lead the directory.
    types = {r.get("ontology_id") for r in rows if not r.get("kind")}
    focus = state.get_meta("recent_followup_ids", [])
    rows.sort(key=lambda r: (focus.index(r["candidate_id"]) if r["candidate_id"] in focus else len(focus),
                             not (r.get("annotations") and r["annotations"]["neighbor"]["ontology_id"] not in types),
                             r.get("n_seeds", 0), r["candidate_id"]))
    available = (limit - len(json_bytes(ctx))) // 2
    entries = ctx["pool_followups"]
    for row in rows:
        entry = index(row)
        cost = len(json_bytes(entry)) + 1
        if cost > available:
            continue
        entries.append(entry)
        available -= cost
    by_id = {r["candidate_id"]: r for r in rows}
    for i, entry in enumerate(entries):
        full = card(by_id[entry["candidate_id"]], compact=True)
        extra = len(json_bytes(full)) - len(json_bytes(entry))
        if extra <= available:
            entries[i] = full
            available -= extra
