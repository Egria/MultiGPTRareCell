"""Strict discovery input boundary and deterministic sparse preprocessing."""
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse


METADATA = {"barcode", "sample_id", "batch", "donor_id", "library_id", "chemistry", "timepoint", "assay"}
SCHEMA = {"counts", "cells", "genes", "orientation", "layer", "species", "tissue", "allowed_metadata"}


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def manifest_path(root, value):
    if not isinstance(value, str) or Path(value).is_absolute():
        raise ValueError("input manifest needs relative paths")
    path = (root / value).resolve()
    if root not in path.parents or not path.is_file():
        raise ValueError("input artifact must exist below the dataset directory")
    return path


def table(path, key, allowed, required):
    with open(path, encoding="utf8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("missing or duplicate TSV header")
        if not set(required).issubset(reader.fieldnames):
            raise ValueError(f"missing required fields in {path.name}: {sorted(required)}")
        cols = [k for k in [key] + list(allowed) if k in reader.fieldnames]
        rows = []
        for row in reader:
            if None in row or any(row[k] is None for k in cols):
                raise ValueError(f"malformed TSV in {path.name}")
            item = {k: row[k] for k in cols}
            if not item[key] or any(any(c in s for c in "\t\r\n\x00") for s in item.values()):
                raise ValueError("empty identifier or control character in input metadata")
            rows.append(item)
    if not rows or len(rows) != len({r[key] for r in rows}):
        raise ValueError(f"missing or duplicate {key}")
    return rows


@dataclass
class Dataset:
    raw: object
    x: sparse.csr_matrix
    cells: list
    genes: list
    info: dict

    @property
    def cell_ids(self):
        return [r["cell_id"] for r in self.cells]

    @property
    def gene_ids(self):
        return [r["gene_id"] for r in self.genes]


def combine_duplicate_entries(a):
    """Sum duplicate coordinates before sparse fancy indexing can coalesce them."""
    if a.has_canonical_format:
        return a
    if a.dtype.kind in "iu":
        dtype = np.uint64 if a.dtype == np.uint64 else np.int64
        bound = int(np.iinfo(dtype).max)
        a = sparse.csr_matrix(
            (a.data.astype(dtype), a.indices.copy(), a.indptr.copy()), shape=a.shape)
        a.sort_indices()
        largest_row = int(np.diff(a.indptr).max())
        if a.nnz and int(a.data.max()) * largest_row > bound:
            for start, stop in zip(a.indptr[:-1], a.indptr[1:]):
                previous, total = -1, 0
                for col, count in zip(a.indices[start:stop], a.data[start:stop]):
                    total = int(count) + (total if col == previous else 0)
                    if total > bound:
                        raise ValueError("Duplicate count exceeds the exact integer storage range")
                    previous = col
    else:
        a = a.astype(np.float64, copy=True)
    a.sum_duplicates()
    if not np.isfinite(a.data).all():
        raise ValueError("Duplicate expression sum is not finite")
    return a


def read_dataset(path, target_sum=10000.0):
    path = Path(path).resolve()
    if path.is_dir():
        path = path / "dataset.yaml"
    spec = yaml.safe_load(path.read_text(encoding="utf8"))
    if not isinstance(spec, dict) or set(spec) - SCHEMA or SCHEMA - {"allowed_metadata"} - set(spec):
        raise ValueError("unexpected or missing dataset.yaml keys")
    if spec["orientation"] != "cells_by_genes" or spec["layer"] not in {"counts", "linear_expression", "log1p_normalized"}:
        raise ValueError("expected cells_by_genes and counts, linear_expression or log1p_normalized")
    if not all(isinstance(spec[key], str) and spec[key].strip() for key in ("species", "tissue")):
        raise ValueError("species and tissue must be nonempty text")
    allowed = spec.get("allowed_metadata", [])
    if not isinstance(allowed, list) or set(allowed) - METADATA or len(allowed) != len(set(allowed)):
        raise ValueError("unapproved discovery metadata")
    if not np.isfinite(target_sum) or target_sum <= 0:
        raise ValueError("target_sum must be finite and positive")
    root = path.parent
    cells = table(manifest_path(root, spec["cells"]), "cell_id", sorted(set(allowed) | {"barcode"}), {"cell_id", *allowed})
    genes = table(manifest_path(root, spec["genes"]), "gene_id", ["symbol", "species"], {"gene_id"})
    a = sparse.load_npz(manifest_path(root, spec["counts"]))
    if a.format != "csr" or a.shape != (len(cells), len(genes)):
        raise ValueError("CSR shape disagrees with cells.tsv / genes.tsv")
    if a.dtype.kind not in "iuf" or not np.all(np.isfinite(a.data)) or np.any(a.data < 0):
        raise ValueError("expression must be finite nonnegative real values")
    a.check_format(full_check=True)
    a = combine_duplicate_entries(a)
    # Canonicalise axes before any stochastic step, graph construction or tie-break.
    ci = sorted(range(len(cells)), key=lambda i: (cells[i].get("barcode") or cells[i]["cell_id"], cells[i]["cell_id"]))
    gi = sorted(range(len(genes)), key=lambda i: genes[i]["gene_id"])
    a = sparse.csr_matrix(a)[ci, :][:, gi]
    a.sum_duplicates()
    a.eliminate_zeros()
    a.sort_indices()
    x = a.astype(np.float64, copy=True)
    totals = np.asarray(x.sum(axis=1)).ravel()
    if not np.all(np.isfinite(totals)):
        raise ValueError("expression totals overflow float64")
    raw = a if spec["layer"] == "counts" else None
    if spec["layer"] in {"counts", "linear_expression"}:
        # Divide by the row total first; multiplying very small libraries by
        # an explicit inverse can overflow before the product is formed.
        x.data /= np.repeat(totals, np.diff(x.indptr))
        np.log1p(x.data * target_sum, out=x.data)
    return Dataset(raw, x, [cells[i] for i in ci], [genes[i] for i in gi],
                   {k: spec[k] for k in ("layer", "species", "tissue")})
