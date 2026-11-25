"""Relationships between supported names and cross-fitted marker-module votes."""
import hashlib
import numpy as np
from .rule_search import positive_otsu


def relationships(specs):
    """Return the transitive broader names of each canonical name."""
    direct = {k: set() for k, v in specs.items() if k == v["canonical"]}
    for name, spec in specs.items():
        child = spec["canonical"]
        direct[child].update(specs[p]["canonical"] for p in spec.get("broader_labels", []))
    ancestors = {}
    def visit(name, path):
        if name in path:
            raise ValueError("Broader labels must describe an acyclic strict subtype relationship")
        if name not in ancestors:
            parents = direct[name]
            ancestors[name] = parents.union(*(visit(p, path | {name}) for p in sorted(parents)))
        return ancestors[name]
    for name in direct:
        visit(name, set())
    return ancestors


class ModuleVotes:
    def __init__(self, x, cells, index):
        self.x = x.tocsc()
        self.index = index
        self.cache = {}
        self.fold = np.array([int.from_bytes(hashlib.sha256(c.encode()).digest()[:4], "big") % 5 for c in cells])

    def hits(self, spec):
        cols = []
        for panel in spec["panels"]:
            key = tuple(sorted(panel))
            if key not in self.cache:
                hit = np.zeros(self.x.shape[0], dtype=bool)
                for gene in key:
                    if gene in self.index:
                        values = self.x[:, self.index[gene]].toarray().ravel()
                        threshold = positive_otsu(values)
                        if threshold is not None:
                            hit |= values >= threshold
                self.cache[key] = hit
            cols.append(self.cache[key])
        return np.column_stack(cols)

    @staticmethod
    def score(matrix, spec, query, own, other):
        def rates(ref):
            return (matrix[ref].sum(0) + .5) / (len(ref) + 1.)
        weight = np.maximum(0., np.log(rates(own) / rates(other)))
        hits = matrix[query]
        ok = hits.sum(1) >= spec.get("min_positive", min(2, matrix.shape[1]))
        if "required_panel" in spec:
            ok &= hits[:, spec["required_panel"] - 1]
        return np.where(ok, np.sort(hits * weight, axis=1)[:, -2:].mean(1), 0.)

    def contrast(self, rows, sides, specs):
        hits = [self.hits(s) for s in specs]
        diff = np.zeros(len(rows))
        for fold in range(5):
            idx = np.flatnonzero(self.fold[rows] == fold)
            a, b = [s[self.fold[s] != fold] for s in sides]
            if len(idx) and min(len(a), len(b)) >= 3:
                diff[idx] = (self.score(hits[0], specs[0], rows[idx], a, b)
                             - self.score(hits[1], specs[1], rows[idx], b, a))
        return diff
