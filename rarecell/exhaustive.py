"""Mandatory singleton measurements and public model assessments for opted-in live runs."""
from collections import Counter

from .contrast_panels import GUIDANCE as PANEL_GUIDANCE
from .config import contract, json_bytes
from .protocol import ProtocolError
from .state import QueryBlocked


def pending(ctrl):
    reviewed = {row[0] for row in ctrl.state.db.execute("SELECT id FROM record WHERE kind='SeedAssessment'")}
    seen = Counter(cid for batch in ctrl.state.all("SeedBatch") for cid in batch["candidate_ids"])
    waiting = [dict(seed, previous_batches=seen[seed["candidate_id"]])
               for seed in ctrl.screened_seeds if seed["candidate_id"] not in reviewed]
    # Preserve the scientific priority within a pass. A missing reply must not hide later seeds.
    return sorted(waiting, key=lambda seed: seed["previous_batches"])


def measure(ctrl, cid):
    s = ctrl.state
    saved = s.get_optional("SeedMeasurement", cid)
    if saved is not None:
        return saved
    try:
        result = ctrl.scientific("expand_candidate", {"candidate_id": cid})
        row = {"candidate_id": cid, "status": result["status"], "evidence_id": result["evidence_id"],
               "query_key": result["query_key"], "summary": result["summary"],
               "proposal_id": result.get("candidate_id")}
    except QueryBlocked as exc:
        row = {"candidate_id": cid, "status": "blocked", "error": str(exc)[:500]}
    s.put("SeedMeasurement", cid, row)
    return row


def prepare_batch(ctrl, mid, queue):
    """Pin exactly the observations sent in this turn, including across interrupted dispatches."""
    s = ctrl.state
    if s.get_optional("SeedBatch", mid) is not None:
        return
    previous = s.get_meta("receipts", [])
    context = {"phase": "exhaustive_seed_review", "mode": ctrl.mode, "prompt_version": ctrl.prompt,
               "contract": contract(), "species": ctrl.manifest.get("species"), "tissue": ctrl.manifest.get("tissue"),
               "population_size": len(ctrl.cells), "total_seeds": len(ctrl.screened_seeds), "pending": len(queue),
               "actions": {"review_seed": "candidate_id, disposition (follow_up, retain, defer or decline), reason; optional evidence_ids."},
               "turn_instruction": "Return one review_seed JSON action per listed seed, in an array. Write a short public "
               "scientific assessment for each: the observed signature or expansion result, what it supports, and what "
               "to test or why to stop. A singleton or no additions can merit follow-up. Specify the candidate ID; the "
               "controller binds the measurement shown. A disposition is a recommendation for later discovery.",
               "previous_errors": [str(r["error"])[:200] for r in previous if "error" in r][:3], "seeds": []}
    if ctrl.labeled:
        context["actions"]["review_seed"] += " For named seeds give guide_label and label_reliability (supported, uncertain, contradicted)."
        context["turn_instruction"] += (" For EACH named seed audit the seed cell itself against its proposed names. "
            "Use tissue/species and your biological knowledge to interpret measured positive and competing markers, "
            "sparsity and plausible alternatives. State which name the cell supports, or why it is ambiguous. "
            "Knowledge can suggest new tests. A supported name guides an all-cell search; justify additions separately. " + PANEL_GUIDANCE)
    size = min(ctrl.budget["seed_review_batch_size"], ctrl.budget["max_actions_per_turn"])
    for seed in queue[:size]:
        cid = seed["candidate_id"]
        anchor = ctrl.candidate(cid)
        measurement = measure(ctrl, cid)
        item = dict(seed, seed_id=ctrl.seed_ids[cid], measurement=measurement, markers=anchor.get("markers", [])[:3])
        label = ctrl.current_label(cid)
        if label:
            item["tentative_label"] = label["label"]
            item["annotation_sources"] = label["sources"]
        context["seeds"].append(item)
        if len(json_bytes(context)) > ctrl.budget["context_bytes"]:
            context["seeds"].pop()
            if not context["seeds"]:
                raise ValueError("context_bytes cannot fit one measured seed for written review")
            break
    ids = [r["candidate_id"] for r in context["seeds"]]
    s.put("SeedBatch", mid, {"message_id": mid, "candidate_ids": ids, "context_ref": s.artifact(context),
                            "includes_new": any(row["previous_batches"] == 0 for row in context["seeds"])})


def record(ctrl, mid, args):
    s = ctrl.state
    batch = s.get_optional("SeedBatch", mid)
    if not batch:
        raise ProtocolError("review_seed requires an assigned review batch")
    cid = args["candidate_id"]
    for assigned in batch["candidate_ids"]:
        if cid == ctrl.seed_ids[assigned]:
            cid = assigned
            break
    if cid not in batch["candidate_ids"]:
        raise ProtocolError("Seed ID is outside the batch in this message")
    obs = s.get("SeedMeasurement", cid)
    eid = obs.get("evidence_id")
    citations = set(args.get("evidence_ids", []))
    label_ids = set(ctrl.current_label(cid).get("evidence_ids", []))
    if citations - ({eid} | label_ids):
        raise ProtocolError("Cite evidence for the assigned seed. Its measurement is bound automatically")
    assessment = {"source": "discovery_agent", "candidate_id": cid, "seed_id": ctrl.seed_ids[cid], "message_id": mid,
                  "evidence_id": eid, "measurement_status": obs["status"], "proposal_id": obs.get("proposal_id"),
                  "disposition": args["disposition"], "reason": args["reason"].strip()}
    names = ctrl.current_label(cid).get("labels", [])
    names = {name.casefold(): name for name in names if name.casefold() != "unknown"}
    if names:
        label = args.get("guide_label", "").strip().casefold()
        name = names.get(label)
        if not label and len(names) == 1:
            name = next(iter(names.values()))
        reliability = args.get("label_reliability", "unassessed")
        if reliability == "unassessed" and sum(cid in b["candidate_ids"] for b in s.all("SeedBatch")) < 2:
            raise ProtocolError("Review the named seed using label_reliability, guide_label and your measured rationale")
        if reliability == "supported" and not name:
            raise ProtocolError("Identify one of the proposed names supported by this seed, or mark it uncertain")
        assessment.update(guide_label=name, label_reliability=reliability)
    if citations & label_ids:
        assessment["annotation_evidence_ids"] = sorted(citations & label_ids)
    old = s.get_optional("SeedAssessment", cid)
    if old and old != assessment:
        raise ProtocolError("This seed already has a written assessment")
    s.put("SeedAssessment", cid, assessment)
    return dict(assessment, action="review_seed", status="recorded")


def page(ctrl, offset, limit):
    order = {"follow_up": 0, "retain": 1, "defer": 2, "decline": 3}
    reviews = sorted(ctrl.state.all("SeedAssessment"), key=lambda r: (order[r["disposition"]], r["candidate_id"]))
    return {"total": len(reviews), "offset": offset, "reviews": reviews[offset:offset+limit]}


def coverage(ctrl):
    s = ctrl.state
    observations = s.all("SeedMeasurement")
    blocked = sum(r["status"] == "blocked" for r in observations)
    assessed = s.all("SeedAssessment")
    completed = sum(r["measurement_status"] != "blocked" for r in assessed)
    total = len(ctrl.seed_ids)
    names = Counter(r["label_reliability"] for r in assessed if "label_reliability" in r)
    return {"named_seed_assessments": dict(names), "total": total, "measured": len(observations) - blocked, "measurement_blocked": blocked,
            "written": len(assessed), "completed": completed, "remaining": total - completed,
            "review_turns": s.get_meta("seed_review_turns", 0), "query_limit": ctrl.query_limit,
            "ordinary_turns": ctrl.budget["max_agent_turns"]}
