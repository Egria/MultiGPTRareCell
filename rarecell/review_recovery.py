"""Recover independently valid adjudication entries from a rejected model batch."""
from functools import partial

from .routed import response


VALIDATION_ERRORS = (ValueError, KeyError, TypeError, AttributeError)


def rejected_reply(message, error):
    parsed = {}
    if message.get("complete", True):
        try:
            parsed = response(message["raw"])
        except VALIDATION_ERRORS:
            pass
    return {"rejected_review": {"message_id": message["message_id"], "error": error,
                                "parsed": parsed if isinstance(parsed, dict) else {}}}


def is_rejected(reply):
    return isinstance(reply, dict) and set(reply) == {"rejected_review"}


def fallback(phase):
    if phase == "conflict_label_audit":
        return "uncertain"
    if phase == "conflict_panel_nomination":
        return None  # No validated biological panel. Use the existing unnamed comparison.
    if phase in {"conflict_measured_review", "unclaimed_conflict_measured_review"}:
        return {"default": "geometry"}
    if phase == "unclaimed_evidence_review":
        return {"add": [], "reason": "Deferred because the individual review could not be validated"}
    raise ValueError("Unexpected review phase: " + phase)


def review_rows(batches, prompt, validate, dispatch, audit):
    if not batches:
        return []
    requests = [(prompt(batch), partial(validate, batch)) for batch in batches]
    replies = dispatch(requests)
    results, retries = [], {}
    for batch, (context, _), reply in zip(batches, requests, replies):
        if not is_rejected(reply):
            results.append(reply)
            continue
        phase = context["phase"]
        field = {"conflict_panel_nomination": "labels", "conflict_measured_review": "rules",
                 "unclaimed_conflict_measured_review": "rules"}.get(phase, "reviews")
        result = {}
        results.append({field: result})
        failure = reply["rejected_review"]
        entries = failure["parsed"].get(field)
        entry = {"phase": phase, "message_id": failure["message_id"], "error": failure["error"],
                 "salvaged": [], "retried": [], "fallback": []}
        audit.append(entry)
        for key, value in batch.items():
            single = {key: value}
            if isinstance(entries, dict) and key in entries:
                try:
                    validate(single, {field: {key: entries[key]}})
                except VALIDATION_ERRORS:
                    pass
                else:
                    result[key] = entries[key]
                    entry["salvaged"].append(key)
                    continue
            if len(batch) > 1:
                entry["retried"].append(key)
                retries[key] = (single, field, result)
            else:
                entry["fallback"].append(key)
                result[key] = fallback(phase)
    if retries:
        singles = [single for single, _, _ in retries.values()]
        repaired = review_rows(singles, prompt, validate, dispatch, audit)
        for key, answer in zip(retries, repaired):
            _, field, destination = retries[key]
            destination[key] = answer[field][key]
    return results


def conservative_relations(proposals):
    """Keep certified direct aliases and strict-parent edges outside cycles."""
    canonical = {name: spec["canonical"] if proposals[spec["canonical"]]["canonical"] == spec["canonical"] else name
                 for name, spec in proposals.items()}
    edges = {name: set() for name in canonical.values()}
    for name, spec in proposals.items():
        edges[canonical[name]].update(canonical[p] for p in spec.get("broader_labels", []))

    def reaches(origin, goal):
        todo = [origin]
        visited = set()
        while todo:
            node = todo.pop()
            if node == goal:
                return True
            if node not in visited:
                visited.add(node)
                todo.extend(edges[node] - visited)
        return False

    return {name: dict(canonical=canonical[name],
                       broader_labels=[p for p in spec.get("broader_labels", [])
                                       if not reaches(canonical[p], canonical[name])])
            for name, spec in proposals.items()}
