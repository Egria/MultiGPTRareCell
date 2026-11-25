"""Cell-wise adjudication of immutable routed lineage nominations."""
from collections import Counter, defaultdict
from itertools import combinations
import numpy as np

from .config import fingerprint
from .lineage_review import ModuleVotes, relationships


def _unit(x):
    return x / np.maximum(np.linalg.norm(x, axis=-1, keepdims=True), 1e-12)


def _geometry(x, sides, shared, seed):
    """No class names or marker symbols enter this calculation."""
    flanks = [np.setdiff1d(side, shared, assume_unique=True) for side in sides]
    refs = [flank if len(flank) else side for flank, side in zip(flanks, sides)]
    # Use exclusive cells for feature selection. Nested sets may have only one flank.
    training = [x[f] for f in flanks if len(f)]
    if not training:
        # Identical anonymous memberships contain no information favouring one side.
        return np.full(len(shared), .5), [0, 0]
    mean = sum(np.asarray(a.mean(0)).ravel() for a in training) / len(training)
    moment = sum(np.asarray(a.multiply(a).mean(0)).ravel() for a in training) / len(training)
    var = moment - mean ** 2
    features = np.flatnonzero(var > 0)
    features = features[np.argsort(-var[features], kind="stable")[:1000]]
    if not len(features):
        return np.full(len(shared), .5), list(map(len, flanks))
    q = _unit(x[shared][:, features].toarray())
    rng = np.random.default_rng(seed)
    draws = min(64, *(len(f) for f in refs))
    similarities = []
    for side, flank in zip(refs, flanks):
        matrix = _unit(x[side][:, features].toarray())
        if len(flank):
            samples = rng.integers(len(side), size=(200, draws))
            centers = _unit(np.stack([matrix[s].mean(0) for s in samples]))
            similarities.append(q @ centers.T)
        else:
            # The nested side has no exclusive cells. Remove each query from its centroid.
            center = _unit(matrix.sum(0) - q)
            similarities.append(np.einsum("ij,ij->i", q, center)[:, None])
    diff = similarities[0] - similarities[1]
    votes = np.mean((diff > 1e-12) + .5 * (np.abs(diff) <= 1e-12), axis=1)
    return votes, list(map(len, flanks))


def build_groups(claims, supports, specs, ids):
    buckets = defaultdict(list)
    for c in sorted(claims, key=lambda v: v["candidate_id"]):
        cid = c["candidate_id"]
        named = supports.get(cid) == "supported" and c["label"] in specs
        key = (True, specs[c["label"]]["canonical"]) if named else (False, cid)
        rows = frozenset(ids[m] for m in c["members"])
        # Same supported identity can unite overlapping drafts transitively.
        parts = buckets[key]
        sources = [cid]
        for i in reversed(range(len(parts))):
            if not rows.isdisjoint(parts[i][0]):
                existing, other = parts.pop(i)
                rows = rows | existing
                sources.extend(other)
        parts.append((rows, sources))
    labels = {c["candidate_id"]: c["label"] for c in claims}
    groups = []
    for (named, key), parts in sorted(buckets.items()):
        for members, sources in parts:
            sources.sort()
            groups.append(dict(id="G" + fingerprint(sources)[:24], parents=sources,
                               label=key if named else labels[sources[0]], reliable=named,
                               members=np.array(sorted(members), dtype=int)))
    return sorted(groups, key=lambda g:g["id"])


def validate_panels(names, specs):
    if not isinstance(specs, dict) or set(specs) != set(names):
        raise ValueError("One panel proposal is required for each requested label")
    for name, s in specs.items():
        canonical = s["canonical"]
        if canonical not in specs or specs[canonical]["canonical"] != canonical:
            raise ValueError("Canonical names must refer directly to a supplied canonical label")
        validate_proposal(name, s, names)
    relationships(specs)


def validate_proposal(name, proposal, names):
    if proposal["canonical"] not in names:
        raise ValueError("Use a supplied canonical name")
    parents = proposal.get("broader_labels", [])
    if not isinstance(parents, list) or any(not isinstance(p, str) or p not in names or p == name for p in parents):
        raise ValueError("broader_labels must list supplied, strictly broader cell types")
    family = proposal.get("lineage_family", "")
    if not isinstance(family, str) or len(family) > 160:
        raise ValueError("lineage_family must be a short shared name")
    contrasts = proposal.get("versus", {})
    if not isinstance(contrasts, dict) or set(contrasts) - (set(names) - {name}):
        raise ValueError("A pair-specific panel must name another supplied class")
    for override in contrasts.values():
        if not isinstance(override, dict) or set(override) - {"panels", "required_panel", "min_positive"}:
            raise ValueError("Invalid pair-specific panel fields")
    for spec in (proposal, *contrasts.values()):
        panels = spec["panels"]
        if not isinstance(panels, list) or not 1 <= len(panels) <= 4:
            raise ValueError("Supply 1 to 4 marker panels")
        for panel in panels:
            if not isinstance(panel, list) or not 1 <= len(panel) <= 8 or not all(isinstance(g, str) and g.strip() for g in panel):
                raise ValueError("Supply nonempty lists of gene symbols")
        genes = [gene for panel in panels for gene in panel]
        if len(set(genes)) != len(genes):
            raise ValueError("Use nonoverlapping panels within each label")
        for field in ("min_positive", "required_panel"):
            if field in spec and (type(spec[field]) is not int or not 1 <= spec[field] <= len(panels)):
                raise ValueError("Panel index or threshold is outside the range of panels")


def _cutoff(values, other, shared):
    flank = np.setdiff1d(other, shared, assume_unique=True)
    if len(flank):
        return np.quantile(values[flank], .95, axis=0)
    # All opposite members are disputed. Use a leave-one-out P95 for each query.
    if len(shared) == 1:
        return np.full(values.shape[1], np.inf)
    data = values[shared]
    order = np.argsort(data, axis=0, kind="stable")
    rank = np.argsort(order, axis=0)
    sorted_data = np.take_along_axis(data, order, axis=0)
    pos = .95 * (len(shared) - 2)
    lo, hi = int(np.floor(pos)), int(np.ceil(pos))
    cols = np.arange(data.shape[1])
    return ((1 - (pos - lo)) * sorted_data[lo + (rank <= lo), cols]
            +(pos - lo) * sorted_data[hi + (rank <= hi), cols])


def _panel_values(x, gene_index, spec):
    scores = []
    missing = []
    for names in spec["panels"]:
        available = [gene_index[g] for g in names if g in gene_index]
        missing.extend(g for g in names if g not in gene_index)
        linear = np.expm1(x[:, available].toarray())
        scores.append(np.log1p(linear.sum(1)))
    return np.stack(scores, axis=1), missing


def measure(x, cells, gene_index, groups, specs, focus=None):
    """Return pair-wise numerical evidence plus compact feedback for isolated biological review."""
    measured_panels = {}
    pairs = []
    cards = {}
    ancestors = relationships(specs)
    modules = None
    for left, right in combinations(groups, 2):
        sides = [left["members"], right["members"]]
        shared = np.intersect1d(*sides, assume_unique=True)
        if not len(shared) or (focus is not None and not focus[shared].any()):
            continue
        key = left["id"] + "/" + right["id"]
        geometry, flanks = _geometry(x, sides, shared, 20261125)
        bio = left["reliable"] and right["reliable"]
        keep = np.arange(len(shared)) if focus is None else np.flatnonzero(focus[shared])
        pair = dict(key=key, groups=(left["id"], right["id"]), rows=shared[keep], geometry=geometry[keep], named=bio)
        pairs.append(pair)
        if not bio:
            continue
        masks = []
        supported = []
        magnitudes = []
        missing = []
        overrides = {}
        effective = []
        for group, other in [(left, right), (right, left)]:
            name = group["label"]
            spec = specs[name]
            override = spec.get("versus", {}).get(other["label"])
            if override:
                overrides[name] = override
                spec = override
            effective.append(spec)
            panel_key = tuple(tuple(panel) for panel in spec["panels"])
            if panel_key not in measured_panels:
                measured_panels[panel_key] = _panel_values(x, gene_index, spec)
            score, absent = measured_panels[panel_key]
            missing.append(absent)
            threshold = _cutoff(score, other["members"], shared)
            positive = (score[shared] > threshold) & (score[shared] > 0)
            sufficient = positive.sum(1) >= spec.get("min_positive", min(2, score.shape[1]))
            if "required_panel" in spec:
                sufficient &= positive[:, spec["required_panel"] - 1]
            positive = positive[keep]
            masks.append(positive @ (1 << np.arange(score.shape[1])))
            supported.append(sufficient[keep])
            magnitudes.append(score[pair["rows"]])
        a, b = left["label"], right["label"]
        fa, fb = (specs[n].get("lineage_family", "").strip().casefold() for n in (a, b))
        if fa and fb and fa != fb and a not in ancestors[b] and b not in ancestors[a]:
            if modules is None:
                modules = ModuleVotes(x, cells, gene_index)
            refs = [g.get("reference", g["members"]) for g in (left, right)]
            pair["module_contrast"] = modules.contrast(pair["rows"], refs, effective)
        patterns = [f"{a}/{b}" for a, b in zip(*masks)]
        pair.update(patterns=patterns, unilateral=supported[0].astype(int) - supported[1].astype(int))
        bins = defaultdict(list)
        for j, pattern in enumerate(patterns):
            bins[pattern].append(j)
        # Agent sees measured per-cell programs grouped by identical detection profiles.
        # Pairs with few shared cells also carry the individual measurements.
        def rounded(row):
            return np.round(row, 2).tolist()
        card = dict(left=left["label"], right=right["label"], exclusive=flanks,
                  shared=len(keep), missing=missing,
                  profiles={pat: [len(ix), int(pair["unilateral"][ix[0]]),
                                 rounded(np.median(magnitudes[0][ix], axis=0)),
                                 rounded(np.median(magnitudes[1][ix], axis=0))]
                            for pat, ix in sorted(bins.items())})
        if overrides:
            card["pair_specific_panels"] = overrides
        if len(keep) <= 12:
            card["cells"] = {str(j): [patterns[j], rounded(magnitudes[0][j]), rounded(magnitudes[1][j])]
                           for j in range(len(keep))}
        cards[key] = card
    return pairs, cards


def validate_rulings(cards, rules):
    if not isinstance(rules, dict) or set(rules) != set(cards):
        raise ValueError("Review every requested pair")
    for key, rule in rules.items():
        if not isinstance(rule, dict) or set(rule) - {"default", "profiles", "cells"}:
            raise ValueError("Invalid biological review fields")
        if rule.get("default") not in ("geometry", "left", "right"):
            raise ValueError("Supply an ambiguous-profile fallback for each pair")
        for field in ("profiles", "cells"):
            overrides = rule.get(field, {})
            if (not isinstance(overrides, dict) or set(overrides) - set(cards[key].get(field, {}))
                    or any(v not in ("geometry", "left", "right") for v in overrides.values())):
                raise ValueError("Decision refers to unmeasured cells or profiles")


def choose(cells, groups, pairs, rules, ancestors=None, fallback_log=None, focus=None):
    """All rivals vote on the original, immutable memberships before one owner is chosen."""
    tally = [{} for _ in cells]
    for group in groups:
        for i in group["members"]:
            tally[i][group["id"]] = [0, 0.]
    broader = {}
    if ancestors:
        for group in groups:
            if group["reliable"]:
                broader[group["id"]] = {p["id"] for p in groups if p["reliable"]
                                        and p["label"] in ancestors.get(group["label"], ())}
    shadowed = {i: {parent for child in rivals for parent in broader.get(child, ()) if parent in rivals}
                for i, rivals in enumerate(tally) if broader and len(rivals) > 1}
    redo = defaultdict(list)
    audit = []
    for p in pairs:
        a, b = p["groups"]
        rule = rules.get(p["key"], {})
        counts = Counter()
        outcomes = []
        module_votes = []
        for j, i in enumerate(p["rows"]):
            geo = p["geometry"][j]
            action = "geometry"
            if p["named"]:
                solo = p["unilateral"][j]
                base = "left" if solo > 0 else "right" if solo < 0 else rule["default"]
                action = rule.get("cells", {}).get(str(j), rule.get("profiles", {}).get(p["patterns"][j], base))
            gap = p["module_contrast"][j] if "module_contrast" in p else 0.
            use_module = abs(gap) > np.log(3.)
            if use_module:
                module_votes.append(dict(cell_id=cells[i], log_ratio_difference=float(gap), prior_action=action))
                action = "left" if gap > 0 else "right"
            tie = abs(geo - .5) <= 1e-12
            left = action == "left" or (action == "geometry" and (geo > .5 or (tie and a < b)))
            winner = a if left else b
            tally[i][winner][0] += 1
            tally[i][a][1] += geo
            tally[i][b][1] += 1 - geo
            mode = "biology" if action != "geometry" else ("stable_id_tie" if tie else "geometry")
            if shadowed.get(i):
                redo[i].append((a, b, winner, geo))
            record = [cells[i], winner, "cross_fitted_modules" if use_module else mode, round(float(geo), 4)]
            if use_module:
                mode = "cross_fitted_modules"
            outcomes.append(record)
            counts[mode] += 1
        record = dict(pair=p["key"], route="named" if p["named"] else "blind", counts=dict(counts), outcomes=outcomes)
        if module_votes:
            record["module_votes"] = module_votes
        audit.append(record)
    owners = {}
    ties = 0
    for i, rivals in enumerate(tally):
        if rivals and (focus is None or focus[i]):
            order = sorted(rivals, key=lambda g:(-rivals[g][0], -rivals[g][1], g))
            initial = order[0]
            removed = []
            while order[0] in shadowed.get(i, ()):
                removed.append(order[0])
                remaining = {g: [0, 0.] for g in rivals if g != order[0]}
                for a, b, winner, geo in redo[i]:
                    if a in remaining and b in remaining:
                        remaining[winner][0] += 1
                        remaining[a][1] += geo
                        remaining[b][1] += 1 - geo
                order = sorted(remaining, key=lambda g: (-remaining[g][0], -remaining[g][1], g))
                rivals = remaining
            if removed and fallback_log is not None:
                fallback_log.append([cells[i], initial, order[0], removed])
            ties += int(len(order) > 1 and rivals[order[0]][0] == rivals[order[1]][0])
            owners[cells[i]] = order[0]
    return owners, audit, ties
