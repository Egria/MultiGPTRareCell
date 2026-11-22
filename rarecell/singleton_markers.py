"""Positive expression contrasts for OntoAnno on frozen singleton candidates."""
import csv

import numpy as np
from scipy import sparse

from .config import fingerprint, read_lines


SOURCE = "singleton_positive_contrast"
PROMPT = "ontoanno-singleton-cl-v1"


def panels(snapshot, top):
    seeds = {r["candidate_id"]: r for r in read_lines(snapshot / "seeds.jsonl")}
    if not seeds:
        return
    with (snapshot / "cells.tsv").open(newline="") as f:
        cells = {r["cell_id"]: i for i, r in enumerate(csv.DictReader(f, delimiter="\t"))}
    with (snapshot / "genes.tsv").open(newline="") as f:
        genes = list(csv.DictReader(f, delimiter="\t"))
    x = sparse.load_npz(snapshot / "data.npz").tocsr()
    totals = detected = None
    for row in read_lines(snapshot / "initial.jsonl"):
        seed = seeds.pop(row["id"], None)
        if seed is None:
            continue
        members = row["members"]
        if members != [seed["cell_id"]]:
            raise ValueError("singleton seed membership mismatch")
        bg = row["marker_background"]
        i = cells[seed["cell_id"]]
        cell = x.getrow(i)
        indices = cell.indices
        if bg["type"] == "local_neighbours":
            neighbours = [cells[c] for c in bg["members"]]
            if i in neighbours or len(set(neighbours)) != len(neighbours):
                raise ValueError("invalid singleton background")
            n = len(neighbours)
            other = x[neighbours][:, indices]
            sums = np.asarray(other.sum(axis=0)).ravel()
            counts = other.getnnz(axis=0)
        elif bg["type"] == "complement":
            if totals is None:
                totals = np.asarray(x.sum(axis=0)).ravel()
                detected = x.getnnz(axis=0)
            n = x.shape[0] - 1
            sums = totals[indices] - cell.data
            counts = detected[indices] - 1
        else:
            raise ValueError("unknown singleton background")
        markers = []
        symbols = set()
        if n and len(indices):
            effects = cell.data - sums / n
            for k in np.lexsort((indices, -effects)):
                if effects[k] <= 1e-12:
                    break
                gene = genes[indices[k]]
                symbol = gene.get("symbol") or gene["gene_id"]
                if symbol in symbols:
                    continue
                symbols.add(symbol)
                markers.append(dict(gene_id=gene["gene_id"], symbol=symbol, effect=float(effects[k]),
                                    pct_in=1.0, pct_out=float(counts[k] / n)))
                if len(markers) == top:
                    break
        yield dict(candidate_id=row["id"], seed_id=seed["id"], n_members=1, resolution=None,
                   partition_id=None, marker_source=SOURCE, cell_set_hash=fingerprint(members),
                   background_type=bg["type"], background_hash=fingerprint(bg), n_background=n,
                   marker_hash=fingerprint(markers), markers=markers, genes=[m["symbol"] for m in markers],
                   **({} if markers else {"reason": "no_expression" if not len(indices) else
                                          "no_background" if not n else "no_positive_markers"}))
    if seeds:
        raise ValueError("singleton seed missing from frozen candidates")
