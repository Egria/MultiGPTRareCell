"""Measured expression boundaries for independent named and anchored discovery sessions."""

from functools import lru_cache
import math

import numpy as np
from scipy import sparse


def positive_otsu(x):
    values, counts = np.unique(x[x > 0], return_counts=True)
    if len(values) < 2:
        return float(values[0]) if len(values) else None
    left = np.cumsum(counts)[:-1]
    right = counts.sum() - left
    sums = np.cumsum(values * counts)[:-1]
    difference = sums / left - (values @ counts - sums) / right
    k = int(np.argmax(left * right * difference**2))
    return float((values[k] + values[k + 1]) / 2)


class ExpressionSearch:
    def __init__(self, ctrl):
        self.cells = np.asarray(ctrl.cells)
        self.rows = sparse.load_npz(ctrl.snapshot / "data.npz").tocsr()
        self.rows.sum_duplicates()
        self.rows.eliminate_zeros()
        self.cols = self.rows.tocsc()
        self.names = [r["gene_id"] for r in ctrl.genes]
        self.index = {g: i for i, g in enumerate(self.names)}
        # A duplicated symbol must be addressed by feature ID.
        aliases = {}
        for i, row in enumerate(ctrl.genes):
            sym = row.get("symbol")
            if sym:
                aliases[sym] = None if sym in aliases else i
        self.index.update({s: i for s, i in aliases.items() if s not in self.index and i is not None})
        self.totals = np.asarray(self.rows.sum(axis=0)).ravel()
        self.detected = np.diff(self.cols.indptr)
        self.column = lru_cache(maxsize=32)(self._column)

    def _column(self, gene):
        if not isinstance(gene, str) or gene not in self.index:
            raise ValueError(f"Unknown or ambiguous gene {gene!r}")
        return self.cols[:, self.index[gene]].toarray().ravel()

    def mask(self, rule, scopes, depth=0):
        if depth > 12:
            raise ValueError("Expression rule nesting exceeds 12")
        if isinstance(rule, str):
            return self.column(rule) > 0
        if isinstance(rule, list) and len(rule) == 2:
            gene, cutoff = rule
            if not isinstance(gene, str) or type(cutoff) not in (float, int) or not math.isfinite(cutoff):
                raise ValueError("Expected [gene, finite cutoff]")
            return self.column(gene) >= cutoff
        if not isinstance(rule, dict) or len(rule) != 1:
            raise ValueError("Expected one expression rule operator")
        key, value = next(iter(rule.items()))
        if key == "scope" and isinstance(value, str) and value in scopes:
            return np.isin(self.cells, scopes[value])
        if key == "otsu" and isinstance(value, str):
            x = self.column(value)
            cut = positive_otsu(x)
            return x >= cut if cut is not None else np.zeros(len(x), bool)
        if key == "not":
            return ~self.mask(value, scopes, depth + 1)
        if key not in ("k", "all", "any") or not isinstance(value, list):
            raise ValueError("Unknown rule or unavailable scope")
        if key == "k":
            if len(value) != 2 or type(value[0]) is not int or not isinstance(value[1], list):
                raise ValueError("Expected {k:[integer,[rules]]}")
            k, children = value
        else:
            children = value
            k = len(children) if key == "all" else 1
        if not 1 <= k <= len(children) <= 64:
            raise ValueError("Invalid k or number of child rules")
        if key == "k" and any(child in children[:i] for i, child in enumerate(children)):
            raise ValueError("k-of-m children must be distinct")
        hits = np.zeros(len(self.cells), dtype=np.int32)
        for child in children:
            hits += self.mask(child, scopes, depth + 1)
        return hits >= k

    def measure(self, request, scopes):
        op = request["action"]
        genes = request.get("genes", [])
        if not isinstance(genes, list) or len(genes) > 64 or not all(isinstance(g, str) for g in genes):
            raise ValueError("Expected at most 64 gene names")
        if op == "survey":
            result = {}
            for g in genes:
                if g not in self.index:
                    result[g] = "missing"
                    continue
                x = self.column(g)
                pos = x[x > 0]
                cut = positive_otsu(x)
                result[g] = dict(
                    detected=len(pos),
                    positive_quantiles=np.quantile(pos, [0.1, 0.5, 0.75, 0.9, 0.95, 0.99, 1]).tolist()
                    if len(pos)
                    else [],
                    otsu=cut,
                    n_high=int((x >= cut).sum()) if cut is not None else 0,
                )
            return {"genes": result}
        if op == "views":
            rules = request.get("rules")
            if not isinstance(rules, list) or not 1 <= len(rules) <= 12:
                raise ValueError("Supply 1 to 12 rules")
            rows = []
            for rule in rules:
                mask = self.mask(rule, scopes)
                row = {"n": int(mask.sum()), "genes": {}}
                for g in genes:
                    if g not in self.index:
                        row["genes"][g] = "missing"
                        continue
                    x = self.column(g)
                    inside, outside = x[mask], x[~mask]
                    row["genes"][g] = [
                        int(np.count_nonzero(inside)),
                        round(float(inside.mean()), 5) if inside.size else None,
                        int(np.count_nonzero(outside)),
                        round(float(outside.mean()), 5) if outside.size else None,
                    ]
                rows.append(row)
            return {"views": rows}
        mask = self.mask(request["rule"], scopes)
        if op == "freeze":
            return {"members": self.cells[mask].tolist()}
        if op != "markers":
            raise ValueError("Unknown scientific action")
        n, outside = int(mask.sum()), int((~mask).sum())
        if not n or not outside:
            return {"n": n, "markers": []}
        subset = self.rows[mask]
        sums = np.asarray(subset.sum(axis=0)).ravel()
        diff = sums / n - (self.totals - sums) / outside
        detected = subset.getnnz(axis=0)
        excluded = request.get("exclude", [])
        limit = request.get("limit", 16)
        if (
            not isinstance(excluded, list)
            or not all(isinstance(g, str) for g in excluded)
            or type(limit) is not int
            or not 1 <= limit <= 40
        ):
            raise ValueError("Expected exclude gene list and limit 1..40")
        skip = {self.index[g] for g in excluded if g in self.index}
        order = [i for i in np.argsort(-diff, kind="stable") if i not in skip][:limit]
        return {
            "n": n,
            "markers": [
                [
                    self.names[i],
                    round(float(diff[i]), 5),
                    round(float(detected[i] / n), 5),
                    round(float((self.detected[i] - detected[i]) / outside), 5),
                ]
                for i in order
            ],
        }
