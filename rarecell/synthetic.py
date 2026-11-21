"""Build a sparse count fixture; the only return value is a discovery manifest path."""
import csv
import hashlib
from pathlib import Path

import numpy as np
from scipy import sparse
import yaml


def write_tsv(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({c: r[c] for c in columns} for r in rows)


def generate(root: Path, seed: int = 1009) -> Path:
    root = Path(root)
    dest, fixture = root / "input", root / "evaluation_fixture"
    if any(d.exists() and any(d.iterdir()) for d in (dest, fixture)):
        raise ValueError("Synthetic generation requires empty input and evaluation_fixture directories")
    rng = np.random.default_rng(seed)
    cells, genes, rare_count = 160, 96, 8
    barcode = ["C" + hashlib.sha256(("cell/%s/%s" % (seed, i)).encode()).hexdigest()[:14]
               for i in range(cells)]
    bg = np.repeat([0, 1], cells // 2)
    injected = np.zeros(cells, dtype=bool)
    for g in (0, 1):
        injected[rng.choice(np.flatnonzero(bg == g), rare_count // 2, replace=False)] = True
    depth = rng.lognormal(0, 0.30, cells)
    rates = np.full((cells, genes), 0.12)
    rates[:, :10] += 1.6
    rates[bg == 0, 10:23] += 3.8
    rates[bg == 1, 23:36] += 3.8
    rates[injected, 36:43] += 8
    counts = sparse.csr_matrix(rng.poisson(rates * depth[:, None]), dtype=np.int32)
    labels = np.where(injected, "synthetic_rare", np.where(bg == 0, "common_A", "common_B"))
    order = rng.permutation(cells)
    dest.mkdir(exist_ok=True, parents=True)
    fixture.mkdir(exist_ok=True, parents=True)
    sparse.save_npz(dest / "counts.npz", counts[order])
    write_tsv(dest / "cells.tsv", ({"cell_id": barcode[i]} for i in order), ["cell_id"])
    write_tsv(dest / "genes.tsv", ({"gene_id": "G%04d" % j} for j in range(genes)), ["gene_id"])
    write_tsv(fixture / "gt.tsv", ({"cell_id": barcode[i], "label": str(labels[i])} for i in order),
              ["cell_id", "label"])
    manifest = {"counts": "counts.npz", "cells": "cells.tsv", "genes": "genes.tsv",
                "orientation": "cells_by_genes", "layer": "counts", "species": "synthetic",
                "tissue": "synthetic", "allowed_metadata": []}
    path = dest / "dataset.yaml"
    path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return path
