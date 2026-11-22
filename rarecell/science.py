"""Sparse, annotation-blind scientific tools for one discovery snapshot."""
import csv
import json
import math
from dataclasses import asdict, dataclass, fields
from numbers import Real
from pathlib import Path

import numpy as np
from scipy import sparse


METHOD_VERSION = "local-contrast/1"
VIEW_NAMES = frozenset(("V_HVG", "V_GINI", "V_LINEAR", "V_TENSOR2"))


@dataclass(frozen=True)
class SearchConfig:
    neighbor_k: int = 32
    background_size: int = 24
    min_background: int = 3
    signature_per_direction: int = 6
    min_effect: float = 0.5
    structure_weight: float = 0.25
    min_expression_score: float = 0.65
    min_score: float = 0.72
    max_cells: int = 50
    block_size: int = 1024

    def __post_init__(self):
        for field in fields(self):
            value = getattr(self, field.name)
            if field.type is int:
                minimum = 0 if field.name == "max_cells" else 1
                if type(value) is not int or value < minimum:
                    raise ValueError(f"{field.name} must be an integer >= {minimum}")
            elif isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
                raise ValueError(f"{field.name} requires a finite number")
            elif value < 0 or (field.name != "min_effect" and value > 1):
                raise ValueError(f"{field.name} outside the allowed range")


def _column(path, name):
    with path.open(encoding="utf-8", newline="") as source:
        return [row[name] for row in csv.DictReader(source, delimiter="\t")]


def _canonical_set(values, index, kind):
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError(f"{kind} must be an explicit nonempty list")
    if any(not isinstance(v, str) or v not in index for v in values):
        raise ValueError(f"unknown {kind} ID")
    return sorted(set(values))


def _rank(indices, scores, barcodes, ids):
    """Tied scores use original barcode, then unique cell ID across samples."""
    keys = (ids[indices], barcodes[indices], -np.round(scores[indices], 12))
    return indices[np.lexsort(keys)]


def _row_units(block):
    block = np.asarray(block, dtype=np.float64)
    if not np.isfinite(block).all():
        raise ValueError("nonfinite view embedding")
    lengths = np.linalg.norm(block, axis=1)
    valid = lengths > 0
    units = np.divide(block, lengths[:, None], out=np.zeros_like(block), where=valid[:, None])
    return units, valid


class Science:
    """Read only declared matrix/manifest fields. Never open evaluation or annotation data."""

    def __init__(self, snapshot: Path, cfg: dict, *, adaptive=False):
        self.adaptive = adaptive
        self.snapshot = Path(snapshot)
        self.config = SearchConfig(**cfg.get("search", {}))
        manifest = json.loads((self.snapshot / "manifest.json").read_text(encoding="utf-8"))
        with (self.snapshot / "cells.tsv").open(encoding="utf-8", newline="") as manifest_cells:
            keys = [(r["cell_id"], r.get("barcode") or r["cell_id"])
                    for r in csv.DictReader(manifest_cells, delimiter="\t")]
        self.cells = np.asarray([c for c, _ in keys], dtype=str)
        self.barcodes = np.asarray([b for _, b in keys], dtype=str)
        self.genes = np.asarray(_column(self.snapshot / "genes.tsv", "gene_id"), dtype=str)
        self.cell_index = {cell: i for i, cell in enumerate(self.cells)}
        self.gene_index = {gene: i for i, gene in enumerate(self.genes)}
        self.data = sparse.load_npz(self.snapshot / "data.npz").tocsr()
        self.data.sum_duplicates()
        self.data.eliminate_zeros()
        self.rna_valid = np.diff(self.data.indptr) > 0
        if self.data.shape != (len(self.cell_index), len(self.gene_index)):
            raise ValueError("snapshot IDs and matrix dimensions disagree")
        self.views = {}
        for name in sorted(manifest["views"]):
            if name not in VIEW_NAMES:
                raise ValueError("unregistered snapshot view")
            embedding = np.load(self.snapshot / "views" / (name + ".npy"), mmap_mode="r", allow_pickle=False)
            if embedding.ndim != 2 or embedding.shape[0] != len(self.cells):
                raise ValueError("view and manifest dimensions disagree")
            self.views[name] = embedding

    def execute(self, action, args):
        """Compute a registered action from canonical (fully resolved) arguments."""
        if action == "retrieve_cells":
            action = "expand_candidate"
        allowed = {
            "coexpress": {"members", "genes", "logic", "threshold", "layer"},
            "compare_groups": {"members", "background", "genes", "layer"},
            "expand_candidate": {"members", "background", "layer"},
        }
        if self.adaptive:
            allowed["coexpress"].add("min_detected")
        if action not in allowed or set(args) - allowed[action]:
            raise ValueError("unknown action or scientific parameter")
        if args.get("layer", "log1p") != "log1p":
            raise ValueError("requested layer is unavailable")
        members = _canonical_set(args.get("members"), self.cell_index, "cell")
        group = self._indices(members)
        if action == "coexpress":
            genes = _canonical_set(args.get("genes"), self.gene_index, "gene")
            return self._coexpress(group, genes, args.get("logic", "AND"), args.get("threshold", 0), args.get("min_detected"))
        background = None
        if "background" in args:
            background = _canonical_set(args["background"], self.cell_index, "cell")
            if set(background).intersection(members):
                raise ValueError("query members and background overlap")
        if action == "expand_candidate":
            return self._expand(group, background)
        genes = _canonical_set(args["genes"], self.gene_index, "gene") if "genes" in args else None
        if background is None:
            background = sorted(self.cell_index.keys() - set(members))
        if not background:
            return {"status": "insufficient_background", "summary": {"n_query": len(group), "n_background": 0},
                    "selection_history": {"members": members, "background": [], "layer": "log1p"}}
        statistics = self._statistics(group, self._indices(background), genes)
        if genes is None:
            selected = self._signature_indices(statistics)
        else:
            selected = list(range(len(genes)))
        records = self._gene_records(statistics, selected)
        return {
            "status": "ok", "statistics": records, "purpose": "exploration",
            "signature_kind": "exploratory_singleton" if len(group) == 1 else "descriptive_contrast",
            "selection_history": {"members": members, "background": background,
                                  "genes": [r["gene_id"] for r in records], "layer": "log1p"},
            "summary": {"n_query": len(group), "n_background": len(background),
                        "n_genes_reported": len(records), "leading_genes": records[:6]},
        }

    def _indices(self, names):
        return np.array([self.cell_index[name] for name in names], dtype=np.intp)

    def _chunks(self, n):
        step = self.config.block_size
        for start in range(0, n, step):
            yield slice(start, min(n, start + step))

    def _coexpress(self, group, genes, logic, threshold, min_detected=None):
        if logic not in ("AND", "OR"):
            raise ValueError("logic requires AND or OR")
        if isinstance(threshold, bool) or not isinstance(threshold, Real) or not math.isfinite(threshold) or threshold < 0:
            raise ValueError("threshold requires a finite nonnegative number")
        required = len(genes) if logic == "AND" else 1
        if min_detected is not None:
            if type(min_detected) is not int or not 1 <= min_detected <= len(genes):
                raise ValueError("min_detected must be within the number of distinct resolved genes")
            required = min_detected
        cols = [self.gene_index[g] for g in genes]
        detected = np.zeros(len(genes), dtype=np.int64)
        hits = []
        for part in self._chunks(len(group)):
            block = self.data[group[part]][:, cols].copy()
            block.data = (block.data > threshold).astype(np.int8)
            block.eliminate_zeros()
            detected += np.bincount(block.indices, minlength=len(genes))
            counts = np.diff(block.indptr)
            selected = counts >= required
            hits.extend(self.cells[group[part][selected]].tolist())
        summary = {"n_query": len(group), "n_matching": len(hits), "fraction": len(hits) / len(group),
                   "logic": logic, "threshold": threshold, "layer": "log1p"}
        if min_detected is not None:
            summary.update(logic="k_of_m", min_detected=required)
        return {"status": "ok", "purpose": "exploration", "matched_members": hits, "summary": summary,
                "gene_counts": dict(zip(genes, detected.tolist())),
                "selection_history": {"members": self.cells[group].tolist(), "genes": genes,
                                      "logic": summary["logic"], "threshold": threshold, "layer": "log1p",
                                      **({"min_detected": required} if min_detected is not None else {})}}

    def _moments(self, indices, cols):
        n_genes = self.data.shape[1] if cols is None else len(cols)
        sums = np.zeros(n_genes)
        detected = np.zeros(n_genes, dtype=np.int64)
        for part in self._chunks(len(indices)):
            block = self.data[indices[part]]
            if cols is not None:
                block = block[:, cols]
            sums += np.bincount(block.indices, weights=block.data, minlength=n_genes)
            detected += np.bincount(block.indices, minlength=n_genes)
        return sums / len(indices), detected / len(indices)

    def _statistics(self, group, background, genes=None):
        cols = None if genes is None else [self.gene_index[g] for g in genes]
        a, pa = self._moments(group, cols)
        b, pb = self._moments(background, cols)
        return {"gene_id": self.genes if genes is None else np.asarray(genes),
                "mean_query": a, "mean_background": b, "fraction_query": pa, "fraction_background": pb,
                "effect": a - b, "detection_difference": pa - pb}

    def _signature_indices(self, statistics):
        effect = statistics["effect"]
        weight = abs(effect) * (1 + abs(statistics["detection_difference"])) / 2
        chosen = []
        for sign in (1, -1):
            rows = np.flatnonzero((sign * effect > 0) & (abs(effect) >= self.config.min_effect))
            order = np.lexsort((statistics["gene_id"][rows], -weight[rows]))
            chosen.extend(rows[order[:self.config.signature_per_direction]].tolist())
        return chosen

    @staticmethod
    def _gene_records(statistics, indices):
        records = []
        for i in indices:
            record = {key: value[i].item() for key, value in statistics.items()}
            record["direction"] = "up" if record["effect"] > 0 else "down" if record["effect"] < 0 else "equal"
            records.append(record)
        return records

    def _query_cosines(self, view, group):
        center = np.zeros(view.shape[1])
        for part in self._chunks(len(group)):
            indices = group[part]
            units, _ = _row_units(view[indices])
            center += units[self.rna_valid[indices]].sum(axis=0)
        length = np.linalg.norm(center)
        if length <= np.finfo(float).eps * len(group):
            return None
        center /= length
        scores = np.full(len(self.cells), np.nan)
        for part in self._chunks(len(self.cells)):
            rows, valid = _row_units(view[part])
            valid &= self.rna_valid[part]
            cosine = np.clip(rows @ center, -1, 1)
            scores[part] = np.where(valid, cosine, np.nan)
        return scores

    def _structure(self, group):
        n = len(self.cells)
        score_sum = np.zeros(n)
        counts = np.zeros(n, dtype=np.int32)
        nominated = np.full(n, -np.inf)
        member_mask = np.zeros(n, dtype=bool)
        member_mask[group] = True
        used = []
        for name, view in self.views.items():
            cosine = self._query_cosines(view, group)
            if cosine is None:
                continue
            used.append(name)
            valid = np.isfinite(cosine)
            score_sum[valid] += np.maximum(0, cosine[valid])
            counts += valid
            eligible = np.flatnonzero(valid & ~member_mask)
            nearest = _rank(eligible, cosine, self.barcodes, self.cells)[:self.config.neighbor_k]
            nominated[nearest] = np.maximum(nominated[nearest], cosine[nearest])
        structure = np.divide(score_sum, counts, out=np.full(n, np.nan), where=counts > 0)
        pool = np.flatnonzero(np.isfinite(nominated))
        local = _rank(pool, nominated, self.barcodes, self.cells)[:self.config.background_size]
        return structure, sorted(self.cells[local].tolist()), used

    def _expression(self, statistics, selected):
        effect = statistics["effect"][selected]
        background = statistics["mean_background"][selected]
        weights = effect ** 2
        positive = effect > 0
        expression = np.empty(len(self.cells))
        positive_score = np.empty(len(self.cells)) if positive.any() else None
        for part in self._chunks(len(self.cells)):
            # Dense block has at most 2 * signature_per_direction columns.
            block = self.data[part, :][:, selected].toarray()
            values = np.clip((block - background) / effect, 0, 1)
            expression[part] = values @ weights / weights.sum()
            if positive_score is not None:
                positive_score[part] = values[:, positive] @ weights[positive] / weights[positive].sum()
        return expression, positive_score

    def _expand(self, group, specified_background):
        members = self.cells[group].tolist()
        structure, local, used = self._structure(group)
        background = local if specified_background is None else specified_background
        history = {"method_version": METHOD_VERSION, "layer": "log1p", "members": members,
                   "background": background, "background_rule": "multiview_local" if specified_background is None else "explicit",
                   "views": used, "parameters": asdict(self.config), "signature": [], "added_members": [],
                   "dependency": "selection-dependent"}
        result = {
            "purpose": "exploration", "signature_kind": "exploratory_singleton" if len(group) == 1 else "descriptive_contrast",
            "proposed_members": members, "added_members": [], "scores": [], "selection_history": history,
            "summary": {"n_query": len(group), "n_background": len(background), "n_added": 0, "views": used},
        }
        if not used:
            result["status"] = "no_valid_structure"
            return result
        if len(background) < self.config.min_background:
            result["status"] = "insufficient_background"
            return result
        statistics = self._statistics(group, self._indices(background))
        selected = self._signature_indices(statistics)
        if not selected:
            result["status"] = "no_signature"
            return result
        signature = self._gene_records(statistics, selected)
        history["signature"] = signature
        expression, pos_score = self._expression(statistics, selected)
        weight = self.config.structure_weight
        scores = weight * structure + (1 - weight) * expression
        eligible = (self.rna_valid & np.isfinite(structure) & (expression >= self.config.min_expression_score)
                    & (scores >= self.config.min_score))
        if pos_score is not None:
            eligible &= pos_score >= self.config.min_expression_score
        eligible[group] = False
        ranked = _rank(np.flatnonzero(eligible), scores, self.barcodes, self.cells)
        returned = ranked[:self.config.max_cells]
        added = self.cells[returned].tolist()
        for i in returned:
            result["scores"].append({"cell_id": str(self.cells[i]), "score": float(scores[i]),
                                     "expression": float(expression[i]),
                                     "positive_expression": float(pos_score[i]) if pos_score is not None else None,
                                     "structure": float(structure[i])})
        result["added_members"] = added
        result["proposed_members"] = sorted(members + added)
        result["status"] = "incomplete" if len(ranked) > len(returned) else "expanded" if added else "no_expansion"
        result["summary"].update(n_added=len(added), n_passing=len(ranked),
                                 signature_genes=[r["gene_id"] for r in signature])
        if self.adaptive:
            valid = self.rna_valid & np.isfinite(structure)
            # Floating-point dot products can exceed 1 by an ulp. Every
            # reported component must be a valid threshold for apply_boundary.
            for component in (scores, expression, structure, pos_score):
                if component is not None:
                    np.clip(component, 0.0, 1.0, out=component)
            ordered = _rank(np.flatnonzero(valid), scores, self.barcodes, self.cells)
            result["ranking"] = [{"cell_id": str(self.cells[i]), "score": float(scores[i]),
                                  "expression": float(expression[i]), "structure": float(structure[i]),
                                  "positive_expression": float(pos_score[i]) if pos_score is not None else None}
                                 for i in ordered]
            gaps = scores[ordered[:-1]] - scores[ordered[1:]]
            largest = np.argsort(-gaps, kind="stable")[:5]
            result["summary"].update(n_scored=len(ordered), n_unscored=len(self.cells)-len(ordered),
                                     rank_gaps=[{"after_rank": int(i+1), "upper": float(scores[ordered[i]]),
                                                 "lower": float(scores[ordered[i+1]])} for i in largest if gaps[i] > 1e-12],
                                     positive_genes=[r["gene_id"] for r in signature if r["direction"] == "up"])
        history.update(added_members=added, n_passing=len(ranked), status=result["status"])
        return result
