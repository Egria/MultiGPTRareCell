"""Bounded, deterministic translation of provider replies into registered requests.

This module has no file, network, state or tool-execution interface. Entity
values survive parsing unchanged. The controller resolves them in its snapshot.
"""
import json
import math
import re
from collections import defaultdict
from decimal import Decimal, InvalidOperation

VERSION = "soft-protocol-1.6"
MAX_CHARS = 512_000
MAX_ACTIONS = 128


class ProtocolError(ValueError):
    """Ambiguous input or an unsupported scientific instruction."""


def _word(value):
    return re.sub(r"[\s_-]+", " ", value.strip().casefold())


def _alias_map(groups):
    result = {}
    for key, aliases in groups.items():
        for spelling in (key, *aliases):
            normalized = _word(spelling)
            if normalized in result and result[normalized] != key:
                raise RuntimeError(f"conflicting built-in alias {spelling}")
            result[normalized] = key
    return result


ACTIONS = _alias_map({
    "review_seed": ("assess_seed", "record_seed_review", "调查种子"),
    "investigate_pool": ("investigate seed pool",),
    "list_seed_reviews": ("list seed assessments",),
    "list_ontology_pools": ("list seed pools",),
    "apply_boundary": (),
    "select_boundary": (),
    "preview_submission": (),
    "list_candidates": ("list candidates", "show candidates", "列出候选", "列出候选群"),
    "get_evidence": ("read evidence", "查看证据", "读取证据"),
    "compare_groups": ("compare groups", "compare populations", "比较群体", "组间比较"),
    "coexpress": ("check co-expression", "co-expression", "check coexpression", "共表达", "检查共表达"),
    "expand_candidate": ("retrieve_cells", "expand candidate", "扩展候选", "检索细胞"),
    "propose_revision": ("revise candidate", "propose revision", "修订候选", "提出修订"),
    "submit": ("submit result", "finalize", "提交", "提交结果"),
})

SEED_DISPOSITIONS = _alias_map({"follow_up": ("investigate", "followup", "深入"), "retain": ("keep", "保留"),
                                "defer": ("uncertain", "待定"), "decline": ("reject", "exclude", "排除")})

FIELDS = _alias_map({
    "disposition": ("recommendation", "outcome", "结论"),
    "seed_id": (),
    "metric": (), "matched_evidence_id": (), "nb_high_gene": (), "min_detected": (), "exclude_candidate_ids": (), "reference_id": (),
    "background_mode": (), "biological_rarity": (), "rarity_reason": (),
    "action": ("operation", "tool", "动作", "操作"),
    "args": ("arguments", "parameters", "params", "参数"),
    "actions": ("requests", "operations", "请求列表", "动作列表"),
    "genes": ("gene_list", "gene_ids", "基因", "基因列表"),
    "candidate_id": ("candidate", "候选", "候选编号"),
    "candidate_ids": ("candidates", "候选列表"),
    "members": ("cell_ids", "cells", "成员", "细胞列表"),
    "scope": ("范围",),
    "background": ("background_cells", "background_members", "背景细胞"),
    "background_id": ("background_candidate", "background_candidate_id", "背景候选"),
    "logic": ("operator", "逻辑"),
    "threshold": ("cutoff", "阈值"),
    "layer": ("expression_layer", "表达层"),
    "evidence_ids": ("evidence_id", "evidence", "证据编号"),
    "offset": ("偏移",),
    "limit": ("page_size", "每页数量"),
    "groups": ("proposals", "群体"),
    "parents": ("parent_ids", "parent_candidates", "父候选"),
    "reason": ("rationale", "理由"),
    "proposed_label": ("label", "暂定标签"),
    "status": ("状态",),
    "run_id": (), "task_id": (), "data_version": (), "contract_version": (), "message_id": (),
    "metadata": ("meta",),
    "note": ("notes", "备注"), "explanation": ("解释",), "comment": ("comments", "说明"),
    "uncertainty": ("不确定性",), "counterevidence": ("反证",), "summary": ("摘要",),
})
EXPLANATIONS = {"note", "explanation", "comment", "uncertainty", "counterevidence", "summary"}
BINDINGS = {"run_id", "task_id", "data_version", "contract_version", "message_id"}
SCOPES = {"candidate_id", "members", "scope"}
BACKGROUNDS = {"background", "background_id"}
ALLOWED = {
    "review_seed": {"candidate_id", "disposition", "reason", "evidence_ids"},
    "investigate_pool": {"candidate_id", "genes", "min_detected", "reason"},
    "list_seed_reviews": {"offset", "limit"},
    "list_ontology_pools": {"offset", "limit"},
    "select_boundary": {"candidate_id", "reason"},
    "apply_boundary": {"candidate_id", "evidence_ids", "metric", "threshold", "reason"},
    "preview_submission": {"candidate_ids", "candidate_id", "reference_id"},
    "list_candidates": {"offset", "limit"},
    "get_evidence": {"candidate_id", "evidence_ids", "metric", "offset", "limit"},
    "compare_groups": SCOPES | BACKGROUNDS | {"genes", "layer", "evidence_ids"},
    "coexpress": SCOPES | {"min_detected", "genes", "logic", "threshold", "layer", "evidence_ids"},
    "expand_candidate": {"candidate_id", "background_mode"} | BACKGROUNDS | {"layer", "evidence_ids"},
    "propose_revision": {"groups", "parents", "reason", "evidence_ids"},
    "submit": {"candidate_ids", "reason", "evidence_ids", "status"},
}


class _Pairs(list):
    """Retain duplicates so a bad object can fail separately inside a batch."""


def _invalid_constant(value):
    # Keep malformed non-finite numbers local to their batch item.
    return float(value)


def _decode(s):
    return json.loads(s, object_pairs_hook=_Pairs, parse_constant=_invalid_constant)


def _unique_pairs(node):
    result = {}
    for key, val in node:
        if key in result:
            raise ProtocolError(f"duplicate JSON field: {key[:80]}")
        result[key] = val
    return result


def _materialize(node):
    if isinstance(node, _Pairs):
        return {k: _materialize(v) for k, v in _unique_pairs(node).items()}
    if isinstance(node, list):
        return [_materialize(v) for v in node]
    if isinstance(node, float) and not math.isfinite(node):
        raise ProtocolError("non-finite value")
    return node


def _field(key):
    name = FIELDS.get(_word(key))
    if name is None:
        raise ProtocolError(f"unsupported field: {key[:80]}")
    return name


def _single(value, key):
    if isinstance(value, list):
        if len(value) != 1:
            raise ProtocolError(f"{key} requires one value")
        return value[0]
    return value


def _string(value, key):
    value = _single(value, key)
    if not isinstance(value, str) or not value.strip():
        raise ProtocolError(f"{key} requires a nonempty string")
    return value


def _ids(value, key, *, allow_empty=False):
    if isinstance(value, str):
        if key == "genes" and re.search(r"[;；]", value):
            value = [token.strip() for token in re.split(r"[;；]", value)]
        else:
            value = [value]
    if not isinstance(value, list) or (not value and not allow_empty):
        raise ProtocolError(f"{key} requires an explicit ID or list")
    if any(not isinstance(item, str) or not item.strip() for item in value):
        raise ProtocolError(f"{key} contains an invalid ID")
    return sorted(set(value))


def _number(value, key, *, integer=False):
    value = _single(value, key)
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ProtocolError(f"{key} requires a number")
    raw = str(value).strip()
    if not re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", raw):
        raise ProtocolError(f"{key} requires a finite number")
    try:
        number = Decimal(raw)
        if number < 0 or not math.isfinite(float(number)):
            raise ProtocolError(f"{key} requires a finite nonnegative number")
        if integer:
            if number != number.to_integral_value() or number > 2**63 - 1:
                raise ProtocolError(f"{key} requires a nonnegative 64-bit integer")
            return int(number)
        result = float(number)
        if result == 0 and number != 0:
            raise ProtocolError(f"{key} is outside the floating-point range")
        return result
    except (OverflowError, InvalidOperation) as exc:
        raise ProtocolError(f"invalid number for {key}") from exc


def _choice(value, key, choices):
    text = _word(_string(value, key))
    if text not in choices:
        raise ProtocolError(f"unsupported {key}: {text[:80]}")
    return choices[text]


def _merge(target, key, value):
    if key in target and target[key] != value:
        raise ProtocolError(f"conflicting values for {key}")
    target[key] = value


def _group(value):
    if not isinstance(value, dict):
        raise ProtocolError("each group requires an object")
    out = {}
    for raw, val in value.items():
        key = _field(raw)
        if key == "evidence_ids":
            clean = _ids(val, key, allow_empty=True)
        elif key in {"members", "candidate_ids", "exclude_candidate_ids"}:
            clean = _ids(val, key)
        elif key in {"matched_evidence_id", "nb_high_gene"}:
            clean = _string(val, key)
        elif key == "candidate_id":
            key, clean = "candidate_ids", [_string(val, key)]
        elif key in EXPLANATIONS or key in {"proposed_label", "rarity_reason", "reason"}:
            clean = _string(val, key)
        elif key == "biological_rarity":
            clean = _choice(val, key, {s: s for s in ("yes", "no", "uncertain")})
        elif key == "status":
            clean = _choice(val, key, {s: s for s in ("confirmed", "provisional", "unresolved")})
        else:
            raise ProtocolError(f"unsupported group field: {key}")
        _merge(out, key, clean)
    if "nb_high_gene" in out and "matched_evidence_id" not in out:
        raise ProtocolError("nb_high_gene requires a matched_evidence_id")
    if not ({"members", "candidate_ids", "matched_evidence_id"} & out.keys()):
        raise ProtocolError("each group requires members, candidate_ids or matched_evidence_id")
    return out


def _value(key, value):
    if key in {"genes", "members", "background", "evidence_ids", "candidate_ids", "parents"}:
        return _ids(value, key, allow_empty=key in {"evidence_ids", "candidate_ids", "parents"})
    if key == "disposition":
        return _choice(value, key, SEED_DISPOSITIONS)
    if key == "metric":
        return _choice(value, key, {_word(name): name for name in ("score", "expression", "positive_expression")})
    if key == "background_mode":
        return _choice(value, key, {name: name for name in ("local", "complement")})
    if key in BINDINGS | {"reason", "candidate_id", "background_id", "reference_id"}:
        return _string(value, key)
    if key == "logic":
        return _choice(value, key, {"and": "AND", "all": "AND", "与": "AND",
                                    "且": "AND", "or": "OR", "any": "OR", "或": "OR"})
    if key == "layer":
        return _choice(value, key, {"log1p": "log1p"})
    if key == "scope":
        return _choice(value, key, {"all": "all", "all cells": "all", "全部细胞": "all"})
    if key in {"limit", "offset", "threshold", "min_detected"}:
        val = _number(value, key, integer=key != "threshold")
        if key in {"limit", "min_detected"} and val == 0:
            raise ProtocolError(f"{key} must be positive")
        return val
    if key == "groups":
        items = [value] if isinstance(value, dict) else value
        if not isinstance(items, list) or not items:
            raise ProtocolError("groups requires a nonempty list of group objects")
        return [_group(item) for item in items]
    if key == "status":
        return _choice(value, key, {s: s for s in ("provisional", "final")})
    raise ProtocolError(f"unsupported parameter: {key}")


def _validate(action, args):
    present = args.keys()
    if action == "investigate_pool" and not {"candidate_id", "genes", "reason"} <= present:
        raise ProtocolError("investigate_pool needs the assigned pool, genes and reason")
    if action == "review_seed" and not {"candidate_id", "disposition", "reason"} <= present:
        raise ProtocolError("review_seed needs a seed candidate_id, disposition, and written reason")
    if len(SCOPES & present) > 1 or len((BACKGROUNDS | {"background_mode"}) & present) > 1:
        raise ProtocolError("choose one foreground scope and one background scope")
    if action in {"compare_groups", "coexpress", "expand_candidate"} and not (SCOPES & present):
        raise ProtocolError("an explicit candidate or cell scope is required")
    if action == "compare_groups" and args.get("scope") == "all":
        raise ProtocolError("compare_groups requires an explicit foreground within the dataset")
    if action == "coexpress":
        if "genes" not in args:
            raise ProtocolError("coexpress requires genes")
        if "min_detected" in args:
            if "logic" in args:
                raise ProtocolError("choose min_detected or logic")
        else:
            args.setdefault("logic", "AND")
        args.setdefault("threshold", 0.0)
    if action == "select_boundary" and not {"candidate_id", "reason"} <= present:
        raise ProtocolError("select_boundary needs a held candidate_id and reason")
    if action == "apply_boundary":
        if not {"candidate_id", "evidence_ids", "metric", "threshold", "reason"} <= present or len(args["evidence_ids"]) != 1:
            raise ProtocolError("apply_boundary needs candidate_id, one expansion evidence_id, metric, threshold and reason")
        if args["threshold"] > 1:
            raise ProtocolError("boundary threshold must be between 0 and 1")
    if action == "get_evidence" and "evidence_ids" in args and "candidate_id" in args:
        raise ProtocolError("choose evidence_ids or candidate_id")
    if action == "propose_revision":
        groups = args.get("groups")
        if not groups:
            raise ProtocolError("revision requires groups")
        for field in ("reason", "evidence_ids"):
            if not args.get(field) and not all(g.get(field) for g in groups):
                raise ProtocolError(f"revision requires {field} on the action or on each group")
        if not args.get("parents") and not (len(groups) == 1 and groups[0].get("candidate_ids")):
            raise ProtocolError("provide explicit parents, or cite candidate_ids for one new group")
    # The one available expression layer was checked before removing the assertion.
    args.pop("layer", None)


def _normalize(obj):
    if not isinstance(obj, dict):
        raise ProtocolError("each action requires an object")
    entries = []
    for key, value in obj.items():
        key = _field(key)
        if key == "action" and isinstance(value, dict):
            # Some providers put the selector and its parameters in an action
            # object. Normalize at this one boundary; unknown keys, duplicate
            # entities and conflicting selectors still pass through all guards.
            if not any(_word(k) in {"type", "name"} or FIELDS.get(_word(k)) == "action" for k in value):
                raise ProtocolError("nested action needs an explicit action selector")
            for inner, item in value.items():
                name = "action" if _word(inner) in {"type", "name"} else _field(inner)
                entries.append((name, item))
        else:
            entries.append((key, value))
    types = set()
    for key, val in entries:
        if key == "action":
            types.add(_choice(val, "action", ACTIONS))
    if not types:
        raise ProtocolError('missing action selector: use {"action":"<registered action>","args":{}} with the actual action and arguments, or an array of such objects')
    if len(types) > 1:
        raise ProtocolError("conflicting action selectors: " + ", ".join(sorted(types)))
    action = types.pop()
    params = []
    metadata = {}
    for key, val in entries:
        if key == "action":
            continue
        if key == "args":
            if not isinstance(val, dict):
                raise ProtocolError("args requires an object")
            params.extend((_field(k), v) for k, v in val.items())
        else:
            params.append((key, val))
    args = {}
    for key, val in params:
        if action == "review_seed" and key == "seed_id":
            key = "candidate_id"
        if action == "review_seed" and key in {"explanation", "summary"} and not any(k == "reason" for k, _ in params):
            key = "reason"
        if key == "metadata":
            if not isinstance(val, dict):
                raise ProtocolError("metadata requires an object")
            for k, v in val.items():
                name = _field(k)
                if name not in EXPLANATIONS:
                    raise ProtocolError(f"{name} cannot be hidden in metadata")
                _merge(metadata, name, _string(v, name))
        elif key in EXPLANATIONS or (key == "reason" and key not in ALLOWED[action]):
            _merge(metadata, key, _string(val, key))
        elif key in ALLOWED[action] | BINDINGS:
            if key == "background" and action == "expand_candidate" and val == "all":
                _merge(args, key, "all")
            else:
                _merge(args, key, _value(key, val))
        else:
            raise ProtocolError(f"unsupported {action} parameter: {key}")
    _validate(action, args)
    for index, group in enumerate(args.get("groups", [])):
        notes = {k: group.pop(k) for k in sorted(group.keys() & EXPLANATIONS)}
        if notes:
            metadata.setdefault("group_notes", []).append({"index": index, **notes})
    return {"action": action, "args": args, "metadata": metadata}


def _json_region(raw):
    # One structured envelope only. Never search inside a malformed parent for
    # an executable child. Brackets after it may begin a truncated second action.
    start = next((i for i, c in enumerate(raw) if c in "[{"), None)
    if start is None:
        raise ProtocolError("supply one JSON action or action list")
    stack = []
    quote = escape = False
    for i in range(start, len(raw)):
        c = raw[i]
        if quote:
            if escape:
                escape = False
            elif c == "\\":
                escape = True
            elif c == '"':
                quote = False
        elif c == '"':
            quote = True
        elif c in "[{":
            stack.append("}" if c == "{" else "]")
            if len(stack) > 32:
                raise ProtocolError("JSON nesting exceeds 32 levels")
        elif c in "]}":
            if not stack or c != stack.pop():
                raise ProtocolError("unbalanced JSON")
            if not stack:
                if any(char in raw[i + 1:] for char in "{}[]"):
                    raise ProtocolError("multiple or incomplete JSON envelopes")
                return raw[start:i + 1]
    raise ProtocolError("incomplete JSON envelope")


def _trim_trailing_commas(text):
    chars = list(text)
    previous = None
    quote = escape = False
    for i, c in enumerate(text):
        if quote:
            if escape:
                escape = False
            elif c == "\\":
                escape = True
            elif c == '"':
                quote = False
                previous = c
        elif c == '"':
            quote = True
        elif c == "," and previous not in (None, "[", "{", ":", ","):
            j = i + 1
            while j < len(text) and text[j] in " \t\r\n":
                j += 1
            if j < len(text) and text[j] in "]}":
                chars[i] = " "
            previous = c
        elif not c.isspace():
            previous = c
    return "".join(chars)


def _line_request(raw):
    """A local one-pass fallback for 'action name' followed by field:value lines."""
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    if not lines or _word(lines[0]) not in ACTIONS:
        raise ProtocolError("supply one JSON action or action list")
    action = _Pairs([("action", lines[0])])
    for line in lines[1:]:
        match = re.fullmatch(r"([^:：]+)[:：]\s*(.+)", line)
        if not match:
            raise ProtocolError("ambiguous field/value line")
        key, value = match.groups()
        if value.startswith(("[", "{", '"')) or value in {"null", "true", "false"}:
            value = _decode(value)
        action.append((key.strip(), value))
    return action


def _error(exc):
    return {"error": str(exc), "status": "needs_clarification"}


def parse_message(raw, *, complete=True, origin="agent", seed_batch=False):
    """Parse a completed, trusted-origin provider reply without executing it."""
    try:
        if origin != "agent" or complete is not True:
            raise ProtocolError("requires a complete message from the discovery agent")
        if not isinstance(raw, str) or not raw.strip() or len(raw) > MAX_CHARS:
            raise ProtocolError(f"message must have 1–{MAX_CHARS} characters")
        first = raw.lstrip().split("\n", 1)[0].strip()
        repaired = None
        if _word(first) in ACTIONS:
            node = _line_request(raw)
            repaired = "field_lines"
        else:
            region = _json_region(raw)
            try:
                node = _decode(region)
            except json.JSONDecodeError:
                cleaned = _trim_trailing_commas(region)
                if cleaned == region:
                    raise
                node = _decode(cleaned)
                repaired = "trailing_comma"
        if isinstance(node, _Pairs):
            # A wrapper can only contain the batch. Parameters on the wrapper
            # cannot silently be ignored or distributed across sibling actions.
            outer = _unique_pairs(node)
            wrapper = _word(next(iter(outer))) if len(outer) == 1 else ""
            if FIELDS.get(wrapper) == "actions" or (seed_batch and wrapper in {"reviews", "assessments"}):
                node = next(iter(outer.values()))
                if isinstance(node, _Pairs) or not isinstance(node, list):
                    raise ProtocolError("actions requires an array")
        batch = [node] if isinstance(node, _Pairs) else node
        if not isinstance(batch, list) or not 1 <= len(batch) <= MAX_ACTIONS:
            raise ProtocolError(f"expected 1–{MAX_ACTIONS} actions")
        result = []
        for item in batch:
            try:
                obj = _materialize(item)
                if seed_batch and isinstance(obj, dict):
                    fields = {FIELDS.get(_word(k)) for k in obj}
                    if ("action" not in fields and "disposition" in fields
                            and fields.intersection({"seed_id", "candidate_id"})):
                        obj["action"] = "review_seed"
                parsed = _normalize(obj)
                if repaired:
                    parsed["metadata"]["format_repair"] = repaired
                result.append(parsed)
            except (ValueError, RecursionError) as exc:
                result.append(_error(exc))
        return result
    except (ValueError, RecursionError) as exc:
        return [_error(exc)]


def resolve_genes(tokens, rows, *, species=None):
    """Resolve IDs or exact symbols using only the caller's versioned manifest.

    An exact stable ID takes precedence in the ID namespace. Symbol collisions
    require the caller to supply a stable ID. No case-folding or edit distance.
    """
    identifiers = set()
    symbols = defaultdict(set)
    for row in rows:
        if species is not None and row.get("species") not in (species, None, ""):
            continue
        gene_id = row["gene_id"]
        identifiers.add(gene_id)
        if row.get("symbol"):
            symbols[row["symbol"]].add(gene_id)
    resolved = set()
    for token in _ids(tokens, "genes"):
        if token in identifiers:
            resolved.add(token)
        elif len(symbols[token]) == 1:
            resolved.update(symbols[token])
        else:
            issue = "ambiguous" if symbols[token] else "unknown"
            raise ProtocolError(f"{issue} gene {token!r}; provide a manifest gene_id")
    return sorted(resolved)
