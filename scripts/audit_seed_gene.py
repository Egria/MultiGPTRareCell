"""Inspect one gene after discovery using the frozen matrix, backgrounds and reference labels."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import sparse

from rarecell.config import read_lines
from rarecell.prepare import verify_snapshot


def table(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def audit(snapshot, gene, gt, label):
    manifest = verify_snapshot(snapshot)
    cells = {r["cell_id"]: i for i, r in enumerate(table(snapshot / "cells.tsv"))}
    genes = [r["gene_id"] for r in table(snapshot / "genes.tsv")]
    if gene not in genes:
        raise ValueError("gene ID missing from snapshot")
    refs = table(gt)
    labels = {r["cell_id"]: r["label"] for r in refs}
    if len(labels) != len(refs) or set(labels) != set(cells) or label not in labels.values():
        raise ValueError("expected one reference label per snapshot cell and a known target label")
    x = sparse.load_npz(snapshot / "data.npz").tocsr()
    j = genes.index(gene)
    y = x[:, j].toarray().ravel()
    classes = defaultdict(list)
    for c, cl in labels.items():
        classes[cl].append(y[cells[c]])
    distribution = {}
    for cl, values in sorted(classes.items()):
        positives = [float(v) for v in values if v > 0]
        distribution[cl] = dict(cells=len(values), detected=len(positives),
                                positive_median=float(np.median(positives)) if positives else None)
    hvg_file = snapshot / "seurat" / "hvg.tsv"
    hvg = next((i for i, r in enumerate(table(hvg_file), 1) if r["gene_id"] == gene), None) if hvg_file.is_file() else None
    gini_file = snapshot / "views" / "gini_features.tsv"
    gini = next((r for r in table(gini_file) if r["gene_id"] == gene), None) if gini_file.is_file() else None
    seeds = {r["candidate_id"]: r["cell_id"] for r in read_lines(snapshot / "seeds.jsonl")}
    totals = None
    target_seeds = []
    for row in read_lines(snapshot / "initial.jsonl"):
        cid = seeds.get(row["id"])
        if cid is None or labels[cid] != label:
            continue
        i = cells[cid]
        bg = row["marker_background"]
        inside = x[i].toarray().ravel()
        if bg["type"] == "local_neighbours":
            indices = [cells[c] for c in bg["members"]]
            n = len(indices)
            sums = np.asarray(x[indices].sum(axis=0)).ravel()
            positive = int(np.count_nonzero(y[indices]))
        elif bg["type"] == "complement":
            if totals is None:
                totals = np.asarray(x.sum(axis=0)).ravel()
            n = x.shape[0] - 1
            sums = totals - inside
            positive = int(np.count_nonzero(y)) - int(y[i] > 0)
        else:
            raise ValueError("unknown singleton background")
        result = dict(cell_id=cid, candidate_id=row["id"], expression=float(y[i]),
                      background=bg["type"], n_background=n, n_background_detected=positive)
        if n:
            effect = inside - sums / n
            def rank(values):
                return 1 + int((values > values[j]).sum()) + int((values[:j] == values[j]).sum())
            result.update(effect=float(effect[j]), absolute_rank=rank(abs(effect)) if abs(effect[j]) > 1e-12 else None,
                          positive_rank=rank(effect) if effect[j] > 1e-12 else None)
        target_seeds.append(result)
    return dict(snapshot_hash=manifest["snapshot_hash"], gene_id=gene, layer="log1p",
                hvg_table_available=hvg_file.is_file(), hvg_rank=hvg, gini=gini,
                distribution=distribution, target_label=label, target_seeds=target_seeds)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--gt", type=Path, required=True)
    parser.add_argument("--gene", required=True, help="exact frozen gene ID")
    parser.add_argument("--label", required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.snapshot, args.gene, args.gt, args.label), indent=2))
