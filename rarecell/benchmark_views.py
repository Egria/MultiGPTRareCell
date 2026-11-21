"""Evaluation-side abundance stress series, with no GT-bearing discovery metadata."""
import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

from scipy import sparse
import yaml

from .data import read_dataset
from .evaluation import read_truth
from .synthetic import write_tsv


def generate_abundance_series(dataset, gt_file, out, target="MAIT", sizes=(144, 72, 36, 18, 9), seed=1009):
    dataset, gt_file, out = [Path(p).resolve() for p in (dataset, gt_file, out)]
    if not sizes or any(type(k) is not int or k < 1 for k in sizes):
        raise ValueError("Sizes must be positive integers")
    if any(a <= b for a, b in zip(sizes, sizes[1:])):
        raise ValueError("Specify strictly decreasing stress sizes")
    if out.exists():
        raise ValueError("Use a fresh output directory")
    src = yaml.safe_load(dataset.read_text(encoding="utf-8"))
    ds = read_dataset(dataset)
    gt = read_truth(gt_file)
    cell_ids = ds.cell_ids
    if set(cell_ids) != set(gt):
        raise ValueError("Input and GT cell identities differ")
    def order(cell):
        return hashlib.sha256(json.dumps([seed, cell], ensure_ascii=False).encode()).digest(), cell
    eligible = sorted((cell for cell, label in gt.items() if label == target), key=order)
    if len(eligible) < max(sizes):
        raise ValueError("GT has fewer target cells than the largest requested size")
    mat = ds.raw if ds.raw is not None else ds.x
    if ds.raw is None and src["layer"] != "log1p_normalized":
        raise ValueError("Cannot reconstruct source expression layer")
    out.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".stress-", dir=out.parent))
    try:
        for i, size in enumerate(sizes):
            kept = set(eligible[:size])
            view = "view" + str(i)
            dest, fixture = work / "discovery" / view, work / "evaluation_fixture" / view
            dest.mkdir(parents=True)
            fixture.mkdir(parents=True)
            selected = [j for j, cell in enumerate(cell_ids) if gt[cell] != target or cell in kept]
            sparse.save_npz(dest / "counts.npz", mat[selected].tocsr())
            # A stress view intentionally strips optional sample metadata as well.
            write_tsv(dest / "cells.tsv", ({"cell_id": cell_ids[j]} for j in selected), ["cell_id"])
            gene_fields = [c for c in ("gene_id", "symbol", "species") if c in ds.genes[0]]
            write_tsv(dest / "genes.tsv", ds.genes, gene_fields)
            manifest = {"counts": "counts.npz", "cells": "cells.tsv", "genes": "genes.tsv",
                        "orientation": "cells_by_genes", "layer": src["layer"],
                        "species": src["species"], "tissue": src["tissue"], "allowed_metadata": []}
            (dest / "dataset.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
            write_tsv(fixture / "gt.tsv",
                      ({"cell_id": cell_ids[j], "label": gt[cell_ids[j]]} for j in selected),
                      ["cell_id", "label"])
            write_tsv(fixture / "retained_target.tsv", ({"cell_id": c} for c in sorted(kept)), ["cell_id"])
            provenance = {"target": target, "size": size, "source_n_target": len(eligible),
                          "seed": seed, "selection": "nested sha256([seed,barcode]) v1"}
            (fixture / "selection.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
        os.rename(work, out)
    finally:
        if work.exists():
            shutil.rmtree(work)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--gt", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--target", default="MAIT")
    parser.add_argument("--seed", type=int, default=1009)
    parser.add_argument("--sizes", type=int, nargs="+", default=[144, 72, 36, 18, 9])
    a = parser.parse_args()
    generate_abundance_series(a.dataset, a.gt, a.out, a.target, a.sizes, a.seed)


if __name__ == "__main__":
    main()
