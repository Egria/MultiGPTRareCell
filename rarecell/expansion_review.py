"""Register a bounded set of exact measured alternatives and review their newly added cells."""
from collections import Counter

from . import adaptive
from .config import contract, fingerprint, json_bytes
from .pool_followup import program_curves
from .protocol import ProtocolError


def register(ctrl, row, data):
    cid, eid = row["candidate_id"], row["evidence_id"]
    core = set(ctrl.candidate(cid)["members"])
    outside_n = len(ctrl.cells) - len(core)
    programs = data.get("programs", {})
    bins = {name: {int(k): set(cells) for k, cells in p["members_by_count"].items()}
            for name, p in programs.items()}
    options = {}

    def add(rule, matches):
        matches = set(matches)
        added = matches - core
        if not added or not matches & core:
            return
        key = fingerprint(sorted(matches))
        if key in options:
            options[key]["rules"].append(rule)
        else:
            # Counts are cumulative in k and refer separately to the original core and the new cells.
            added_hits = {}
            for name, counts in bins.items():
                counts_at_k = [len(added & counts[k]) for k in sorted(counts)]
                added_hits[name] = [sum(counts_at_k[i:]) for i in range(len(counts_at_k))]
            c = ctrl.make_revision(matches, [cid], row["label"], "provisional", [eid], "measured_alternative", rule)
            options[key] = dict(candidate_id=c["id"], n=len(matches), retained_core=len(matches & core),
                n_added=len(added), added_hits_by_k=added_hits, rules=[rule])

    def priority(inside, outside):
        # A shortlist score for seed recovery versus expansion size, without reference labels.
        return inside / (len(core) + inside + outside)

    for name, p in programs.items():
        if name == "counter":
            continue
        cutoffs = [v for v in p["curve"] if v["outside_hits"] and
                   v["seed_hits"] * outside_n > v["outside_hits"] * len(core)]
        if cutoffs:
            best = max(cutoffs, key=lambda v: (priority(v["seed_hits"], v["outside_hits"]), v["seed_hits"], v["k"]))
            k = best["k"]
            add(dict(program=name, min_detected=k), (c for n, cells in bins[name].items() if n >= k for c in cells))

    positive = {g for name, p in programs.items() if name != "counter" for g in p["genes"]}
    fits = [(priority(f["high_seeds"], f["high_outside"]), gene) for gene, f in row.get("contrasts", {}).items()
            if (not programs or gene in positive) and f.get("effect", 0) > 0 and f.get("high_outside", 0)
            and f.get("high_seeds", 0) * outside_n > f["high_outside"] * len(core)]
    for _, gene in sorted(fits, key=lambda v: (-v[0], v[1]))[:2]:
        add(dict(nb_high_gene=gene), adaptive.nb_high_members(ctrl, eid, gene))
    if "matched_members" in data:
        add(dict(genes=row["genes"], min_detected=row["min_detected"]), data["matched_members"])
    return list(options.values())


def prepare_batch(ctrl, mid):
    s = ctrl.state
    if s.get_optional("ExpansionBatch", mid):
        return True
    if s.get_meta("expansion_phase_complete", False):
        return False
    attempts = Counter(cid for b in s.all("ExpansionBatch") for cid in b["candidate_ids"])
    reviewed = {r["candidate_id"] for r in s.all("ExpansionReview")}
    pending = [r for r in s.all("PoolFollowup") if r.get("boundaries") and r["candidate_id"] not in reviewed]
    if not pending:
        s.set_meta("expansion_phase_complete", True)
        return False
    pending.sort(key=lambda r: (attempts[r["candidate_id"]], r["n_seeds"], r["candidate_id"]))
    ctx = dict(phase="expansion_review", mode=ctrl.mode, prompt_version=ctrl.prompt, contract=contract(),
        species=ctrl.manifest.get("species"), tissue=ctrl.manifest.get("tissue"), cases=[],
        actions={"review_expansion": "candidate_id (assigned pool), choice (draft ID, core, or defer), reason; optional include_core boolean"},
        turn_instruction="Return one review_expansion per case. Evaluate the proposed name on the original cells. "
        "Use biological knowledge to interpret the measured programs. For each plausible wider rule compare the added "
        "cells with the core using support and counter hits; arrays are cumulative for k=1,2,... . A selecting gene "
        "or reused panel is selection-dependent. Consider a provisional expansion when evidence supports additions. "
        "Choose the draft ID, core alone, or defer. Explain the evidence in the added cells for widening or narrowing, "
        "including missing independent markers. Choosing a draft or core saves it for ordinary discovery and submission. "
        "When a measured rule omits well-supported original cells, include_core=true can retain the full original core. "
        "A different cutoff or subset can be proposed in ordinary discovery using the evidence ID.")
    ids = []
    size = min(ctrl.budget["max_actions_per_turn"], ctrl.budget["seed_review_batch_size"])
    for row in pending[:size]:
        case = {k: row[k] for k in ("candidate_id", "label", "evidence_id", "n_seeds", "boundaries")}
        case["programs"] = program_curves(row)
        genes = {rule["nb_high_gene"] for b in case["boundaries"] for rule in b["rules"] if "nb_high_gene" in rule}
        case["nb_fits"] = {g: row["contrasts"][g] for g in sorted(genes)}
        if row.get("seed_name_assessments"):
            case["seed_name_assessments"] = row["seed_name_assessments"]
        if row.get("written_questions"):
            case["seed_notes"] = row["written_questions"][:2]
        ctx["cases"].append(case)
        if len(json_bytes(ctx)) > ctrl.budget["context_bytes"]:
            ctx["cases"].pop()
            if not ids:
                raise ValueError("context_bytes cannot fit one measured boundary review")
            break
        ids.append(row["candidate_id"])
    s.put("ExpansionBatch", mid, dict(candidate_ids=ids, context_ref=s.artifact(ctx)))
    return True


def record(ctrl, mid, args):
    s = ctrl.state
    batch = s.get_optional("ExpansionBatch", mid)
    cid, choice = args["candidate_id"], args["choice"].strip()
    if not batch or cid not in batch["candidate_ids"]:
        raise ProtocolError("Review an assigned expansion case")
    row = s.get("PoolFollowup", cid)
    if s.get_optional("ExpansionReview", cid):
        raise ProtocolError("This expansion case already has a review")
    review = dict(candidate_id=cid, choice=choice, reason=args["reason"], message_id=mid, source="discovery_agent")
    alternatives = {b["candidate_id"] for b in row["boundaries"]}
    if choice.casefold() in {"core", "defer"}:
        choice = choice.casefold()
    elif choice not in alternatives:
        raise ProtocolError("Choose this case's exact draft ID, core or defer")
    review["choice"] = choice
    include_core = args.get("include_core", False) and choice not in {"core", "defer"}
    review["include_core"] = include_core
    if choice != "defer":
        selected = cid if choice == "core" else choice
        result = ctrl.decide("propose_revision", {"parents": [cid], "reason": args["reason"],
            "evidence_ids": [row["evidence_id"]], "groups": [{"candidate_ids": [selected, cid] if include_core else [selected],
            "proposed_label": row["label"], "status": "provisional"}]}, mid)
        review["working_candidate_id"] = result["revised_candidate_ids"][0]
        review["decision_id"] = result["decision_id"]
    s.put("ExpansionReview", cid, review)
    return dict(review, action="review_expansion", candidate_id=cid, status="recorded")


def finish_batch(ctrl, mid, batch):
    # A missing reply gets a second opportunity without displacing unreviewed cases.
    s = ctrl.state
    attempts = Counter(cid for b in s.all("ExpansionBatch") for cid in b["candidate_ids"])
    for cid in batch["candidate_ids"]:
        if not s.get_optional("ExpansionReview", cid) and attempts[cid] >= 2:
            s.put("ExpansionReview", cid, {"candidate_id": cid,
                "source": "controller", "choice": "unreviewed", "message_id": mid,
                "reason": "No usable review received in two assigned batches."})
