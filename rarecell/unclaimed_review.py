"""Offer previously measured boundaries to the agent for still unclaimed cells."""
from collections import defaultdict
import numpy as np

from .config import fingerprint, json_bytes
from .rule_search import positive_otsu


def used_scopes(rule):
    if isinstance(rule, list):
        return set().union(*(used_scopes(item) for item in rule))
    if isinstance(rule, dict):
        if "scope" in rule:
            return {rule["scope"]}
        return set().union(*(used_scopes(v) for v in rule.values()))
    return set()


def measured_rules(route, claims):
    jobs = {j["id"]: j for j in route.jobs()}
    frozen = {r["job"] for r in claims}
    scopes = {}
    records = {}
    for row in route.s.all("RouteExchange"):
        job = row["job"]
        if job not in frozen:
            continue
        exchange = row["exchange"]
        req, resp = exchange["request"], exchange["response"]
        if req["action"] == "views":
            views = zip(req["rules"], resp["views"])
        elif req["action"] == "markers":
            views = [(req["rule"], resp)]
        else:
            continue
        for rule, view in views:
            needed = used_scopes(rule)
            if needed and job not in scopes:
                scopes[job] = route.scope_members(jobs[job])
            entry = dict(job=job, query_key=exchange["query_key"], rule=rule, n=view["n"],
                         scopes={k: scopes[job][k] for k in sorted(needed)},
                         genes=req.get("genes", []) if req["action"] == "views" else
                               [r[0] for r in resp.get("markers", [])[:8]])
            records[fingerprint(entry)] = entry
    return [records[k] for k in sorted(records)]


def review_options(expr, groups, claims, queries, owners, review_batches, limit):
    """Approve an exact former hit set or identified members; defer by default."""
    from_job = {c["job"]: c for c in claims}
    by_parent = {cid: g for g in groups for cid in g["parents"]}
    unclaimed = np.array([c not in owners for c in expr.cells], dtype=bool)
    cards, members, thresholds = {}, {}, {}
    if not unclaimed.any():
        return {}, {}
    for row in queries:
        claim = from_job[row["job"]]
        g = by_parent[claim["candidate_id"]]
        mask = expr.mask(row["rule"], row["scopes"])
        if int(mask.sum()) != row["n"]:
            raise ValueError("Archived boundary size differs from its source query")
        additional = mask & unclaimed
        n_new = int(additional.sum())
        core_hits = int(mask[g["members"]].sum())
        if not n_new or not core_hits:
            continue
        key = fingerprint([g["id"], row["query_key"], row["rule"]])
        expression = {}
        new = expr.cells[additional].tolist()
        cell_high = {c: [] for c in new} if n_new <= 12 else {}
        for gene in row.get("genes", [])[:12]:
            if gene not in expr.index:
                continue
            col = expr.column(gene)
            if gene not in thresholds:
                thresholds[gene] = positive_otsu(col)
            threshold = thresholds[gene]
            if threshold is not None:
                high = col >= threshold
                expression[gene] = [int(high[g["members"]].sum()), int(high[additional].sum())]
                if cell_high:
                    for c in expr.cells[high & additional]:
                        cell_high[c].append(gene)
        cards[key] = dict(core_rationale=claim.get("reason", "")[:320],
                          high_counts_core_and_additions=expression, group=g["id"],
                          label=g["label"], identity_supported=g["reliable"],
                          query_key=row["query_key"], rule=row["rule"], n_hits=row["n"],
                          core_hits=core_hits, core_size=len(g["members"]), n_unclaimed=len(new),
                          other_nominated_hits=row["n"]-core_hits-len(new),
                          examples=new[:12], membership_sha256=fingerprint(sorted(new)))
        if cell_high and expression:
            cards[key]["addition_high_genes"] = cell_high
        members[key] = np.packbits(additional)
    indices = {c:i for i,c in enumerate(expr.cells)}
    def is_eligible(key, cell):
        if not isinstance(cell, str) or cell not in indices:
            return False
        byte, bit = divmod(indices[cell], 8)
        return bool(members[key][byte] & (1 << (7 - bit)))

    def prompt(batch):
        return dict(phase="unclaimed_evidence_review", options=batch,
                    instruction='Review previously measured hits beyond all frozen nominations. Per-gene counts use local positive-expression Otsu high states, in core and additions respectively. Inspect the measured genes, selection dependence, other nominated hits and original identity uncertainty. examples lists up to 12 exact unclaimed IDs; "all" means every unclaimed hit of that rule. Return {"reviews":{option_key:{"add":"all" or [cell IDs] or [], "reason":"..."}}}. Several rules from the same RNA can be correlated. Approve only a supported extension; choose [] to defer. A provisional group can retain an uncertain name.')
    approvals = defaultdict(set)
    decisions = {}
    visible = {}
    for key, card in cards.items():
        if len(json_bytes(prompt({key: card}))) <= limit:
            visible[key] = card
        else:
            decisions[key] = dict(card, decision=dict(add=[], reason="Full measured rule exceeds context_bytes"),
                                  n_accepted=0)
    def validate(batch, reply):
        reviews = reply.get("reviews")
        if not isinstance(reviews, dict) or reviews.keys() != batch.keys():
            raise ValueError("Review every supplied measured boundary")
        for key, decision in reviews.items():
            add = decision["add"]
            if add != "all" and (not isinstance(add, list) or not all(is_eligible(key, c) for c in add)):
                raise ValueError("Select only measured unclaimed hits for this option")
            if not isinstance(decision.get("reason"), str) or not decision["reason"].strip():
                raise ValueError("Explain each unclaimed membership decision")
    for reply in review_batches(visible, prompt, 8, validate):
        for key, decision in reply["reviews"].items():
            if decision["add"] == "all":
                mask = np.unpackbits(members[key], count=len(expr.cells)).astype(bool)
                accepted = expr.cells[mask].tolist()
            else:
                accepted = set(decision["add"])
            approvals[cards[key]["group"]].update(accepted)
            decisions[key] = dict(cards[key], decision=decision, n_accepted=len(accepted))
    return approvals, decisions
