"""Export locked single-agent benchmark releases as GUI-ready RareCell datasets."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy import sparse
import yaml

from .config import hash_file, read_json, write_json
from .data import write_tsv


DATASETS = {
    "pbmc_tcr": ("human", "peripheral blood", "counts"),
    "cellsius": ("human", "cell line mixture", "counts"),
    "cellsius_spheroid": ("human", "cortical spheroid", "counts"),
    "villani": ("human", "peripheral blood", "linear_expression"),
    "baron": ("human", "pancreas", "counts"),
    "gse156490": ("human", "pancreas", "linear_expression"),
    "human_pgc": ("human", "gastrula", "linear_expression"),
    "airway": ("mouse", "tracheal epithelium", "counts"),
    "pbmc68k": ("human", "peripheral blood", "counts"),
}
ALIASES = {"pgc": "human_pgc", "spheroid": "cellsius_spheroid"}
MAIT_SIZES = dict(zip(("view" + str(i) for i in range(5)), (144, 72, 36, 18, 9)))


def rows(path):
    with path.open(encoding="utf8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def ids(path, key):
    result = [row[key] for row in rows(path)]
    if not result or not all(result) or len(set(result)) != len(result):
        raise ValueError(f"Missing or duplicate {key} in {path}")
    return result


def default_views(name, table):
    selected = [r["view_id"] for r in table if r["role"] == "main"]
    if name == "pbmc_tcr":
        selected = list(MAIT_SIZES)
    if name == "cellsius" and len(selected) != 6:
        raise ValueError("Expected six main CellSIUS conditions")
    if not selected or not set(selected) <= {r["view_id"] for r in table}:
        raise ValueError(f"Required {name} views missing. Rebuild with the pinned single-agent preparer in a fresh project.")
    return selected


def export_release(project, out, name, all_views=False):
    name = ALIASES.get(name, name)
    project, out = Path(project).resolve(), Path(out).resolve()
    if name not in DATASETS:
        raise ValueError(f"Unknown dataset {name}")
    data, bench = project / "data" / name, project / "benchmark" / name
    lockfile = bench / "release.lock.json"
    lock = read_json(lockfile)["artifacts"]
    checked = set()

    def verified(path):
        rel = str(path.relative_to(project))
        if path not in checked:
            if rel not in lock or hash_file(path) != lock[rel]:
                raise ValueError(f"Single-agent release checksum mismatch: {rel}")
            checked.add(path)
        return path

    table = list(rows(verified(data / "views.tsv")))
    sizes = {r["view_id"]: int(r["n_cells"]) for r in table}
    views = default_views(name, table)
    if all_views:
        views = list(sizes)
    if len(sizes) != len(table) or any(not v or Path(v).name != v or v in {".", ".."} for v in views):
        raise ValueError("Invalid source view IDs")
    manifest = read_json(verified(bench / "manifest.json"))
    identity = {"format_version": 1, "dataset": name, "release": manifest["gt_release_id"],
                "source_lock_sha256": hash_file(lockfile), "views": views}
    dest = out / name
    if dest.exists():
        saved = read_json(dest / "export.json")
        if saved["source"] != identity:
            raise ValueError(f"Different release or view selection already at {dest}; choose another --out")
        required = {f"{v}/{item}" for v in views for item in
                    ("input/dataset.yaml", "input/cells.tsv", "input/genes.tsv", "input/expression.npz", "evaluation_fixture/gt.tsv")}
        if set(saved["files"]) != required:
            raise ValueError(f"Incomplete export inventory at {dest}")
        for rel, digest in saved["files"].items():
            if hash_file(dest / rel) != digest:
                raise ValueError(f"Changed exported file: {dest / rel}")
        return dest

    species, tissue, layer = DATASETS[name]
    spec = read_json(verified(data / "input.json"))
    if spec["is_raw_count"] != (layer == "counts"):
        raise ValueError("Assay disagrees with registered release type")
    genes = list(rows(verified(data / "genes.tsv")))
    gene_ids = [g["feature_id"] for g in genes]
    if not genes or not all(gene_ids) or len(set(gene_ids)) != len(genes):
        raise ValueError("Missing or duplicate feature IDs")
    base_matrix, base_cells = None, None
    if name != "cellsius":
        base_matrix = sparse.load_npz(verified(data / spec["matrix"])).tocsr()
        base_cells = np.array(ids(verified(data / "cells.tsv"), "cell_id"))
        if base_matrix.shape != (len(base_cells), len(genes)):
            raise ValueError("Full expression axes disagree with release")
    scoped = name == "pbmc_tcr" or (name == "cellsius_spheroid" and
            read_json(verified(bench / "source_lock.json")).get("evaluation_profile") == "spheroid_cp_target_f1_v1")
    ref = verified(bench / ("membership.tsv" if scoped else "labels_source.tsv"))
    fields = ["cell_id", "target_id", "reference_kind", "state"] if scoped else ["cell_id", "label"]
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out, prefix="." + name + "-") as tmp:
        stage = Path(tmp) / name
        for v in views:
            folder = data / "views" / v
            cells = ids(verified(folder / "cells.tsv"), "cell_id")
            if name == "cellsius":
                source = read_json(verified(folder / "input.json"))
                matrix = sparse.load_npz(verified(folder / source["matrix"])).tocsr()
            else:
                idx = np.load(verified(folder / "row_indices.npy"), allow_pickle=False)
                if (idx.ndim != 1 or idx.dtype.kind not in "iu" or not len(idx) or idx.min() < 0
                        or idx.max() >= len(base_cells) or base_cells[idx].tolist() != cells):
                    raise ValueError(f"RNA row mapping differs from cell IDs in {v}")
                matrix = base_matrix[idx]
            if (matrix.shape != (len(cells), len(genes)) or len(cells) != sizes[v]
                    or not np.isfinite(matrix.data).all() or (matrix.data < 0).any()):
                raise ValueError(f"Invalid expression in {v}")
            if cells != ids(verified(bench / "views" / v / "cells.tsv"), "cell_id"):
                raise ValueError(f"Benchmark cell list differs in {v}")
            inputs, fixture = stage / v / "input", stage / v / "evaluation_fixture"
            inputs.mkdir(parents=True)
            fixture.mkdir()
            sparse.save_npz(inputs / "expression.npz", matrix)
            write_tsv(inputs / "cells.tsv", ["cell_id"], ({"cell_id": c} for c in cells))
            write_tsv(inputs / "genes.tsv", ["gene_id", "symbol"],
                      ({"gene_id": g["feature_id"], "symbol": g["gene_symbol"]} for g in genes))
            cfg = dict(counts="expression.npz", cells="cells.tsv", genes="genes.tsv", orientation="cells_by_genes",
                       layer=layer, species=species, tissue=tissue, allowed_metadata=[])
            (inputs / "dataset.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf8")
            truth = []
            for r in rows(ref):
                if r["view_id"] != v:
                    continue
                if scoped:
                    truth.append({k: r[k] for k in fields})
                elif r["level_id"] == "fine":
                    truth.append({"cell_id": r["cell_id"], "label": r["original_name"]})
            cell_set = set(cells)
            if {r["cell_id"] for r in truth} != cell_set:
                raise ValueError(f"Incomplete ground truth in {v}")
            if scoped:
                targets = {(r["target_id"], r["reference_kind"]) for r in truth}
                if len(truth) != len(cells) * len(targets) or len({tuple(r[k] for k in fields[:3]) for r in truth}) != len(truth):
                    raise ValueError(f"Incomplete target memberships in {v}")
                if any(r["state"] not in {"positive","negative","unknown","disputed"} for r in truth):
                    raise ValueError(f"Unknown reference state in {v}")
                if name == "pbmc_tcr" and v in MAIT_SIZES:
                    if sum(r["target_id"] == "MAIT" and r["state"] == "positive" for r in truth) != MAIT_SIZES[v]:
                        raise ValueError(f"MAIT count differs from frozen recipe in {v}")
            elif len(truth) != len(cells) or any(not r["label"] for r in truth):
                raise ValueError(f"Missing or duplicate fine labels in {v}")
            write_tsv(fixture / "gt.tsv", fields, truth)
        files = {str(path.relative_to(stage)): hash_file(path) for path in sorted(stage.rglob("*")) if path.is_file()}
        write_json(stage / "export.json", {"source": identity, "files": files})
        os.rename(stage, dest)
    return dest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True, help="Single-agent project containing data/ and benchmark/")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--all-views", action="store_true", help="Include optional views already present in the source release")
    parser.add_argument("datasets", nargs="+", choices=sorted(set(DATASETS) | set(ALIASES)))
    args = parser.parse_args()
    failed = False
    for name in args.datasets:
        try:
            dest = export_release(args.project, args.out, name, args.all_views)
            for v in read_json(dest / "export.json")["source"]["views"]:
                print(f"{dest/v/'input/dataset.yaml'}\t{dest/v/'evaluation_fixture/gt.tsv'}", flush=True)
        except (OSError, ValueError, KeyError) as exc:
            print(f"{name}: {exc}", file=sys.stderr, flush=True)
            failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
