"""Import single-agent RNA views from a copied data tree or a frozen benchmark release."""
import argparse
import csv
from contextlib import ExitStack
import os
from functools import cached_property
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy import sparse
import yaml

from .config import hash_file, read_json, write_json
from .data import table as read_table, write_tsv


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
    "hao": ("human", "peripheral blood", "counts"),
    "larry": ("mouse", "cultured hematopoietic progenitors", "linear_expression"),
    "mmochi": ("human", "peripheral blood", "counts"),
}
INPUT_FILES = ("dataset.yaml", "cells.tsv", "genes.tsv", "expression.npz")
COPIED_OUT = Path("work/copied-inputs")
RELEASE_OUT = Path("work/single-agent-datasets")
ALIASES = {"pgc": "human_pgc", "spheroid": "cellsius_spheroid"}
MAIT_SIZES = dict(zip(("view" + str(i) for i in range(5)), (144, 72, 36, 18, 9)))


def rows(path):
    with path.open(encoding="utf8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def ids(path, key):
    return [r[key] for r in read_table(path, key, [], {key})]


def default_views(name, table):
    selected = [r["view_id"] for r in table if r["role"] == "main"]
    if name == "pbmc_tcr":
        selected = list(MAIT_SIZES)
    if name == "cellsius" and len(selected) != 6:
        raise ValueError("Expected six main CellSIUS conditions")
    if not selected or not set(selected) <= {r["view_id"] for r in table}:
        raise ValueError(f"Required {name} views missing. Rebuild with the pinned single-agent preparer in a fresh project.")
    return selected


def view_sizes(data, verified):
    entries = list(rows(verified(data / "views.tsv")))
    sizes = {r["view_id"]: int(r["n_cells"]) for r in entries}
    if (len(sizes) != len(entries) or not sizes or
            any(not v or Path(v).name != v or v in {".", ".."} or n < 1 for v, n in sizes.items())):
        raise ValueError(f"Invalid view IDs or sizes in {data}/views.tsv")
    return entries, sizes


class ExpressionSource:
    """Read only the axes and selected assay defined by the single-agent preparer."""

    def __init__(self, data, name, verified):
        self.data, self.name, self.verified = data, name, verified
        self.species, self.tissue, self.layer = DATASETS[name]
        self.spec = read_json(verified(data / "input.json"))
        if self.spec["is_raw_count"] is not (self.layer == "counts"):
            raise ValueError(f"Assay disagrees with registered expression layer in {data}")
        if self.spec.get("ground_truth_included") is True:
            raise ValueError(f"Expected a label-free input manifest in {data}")
        self.entries, self.sizes = view_sizes(data, verified)

    @cached_property
    def genes(self):
        return read_table(self.verified(self.data / "genes.tsv"), "feature_id",
                          ["gene_symbol"], {"feature_id", "gene_symbol"})

    @cached_property
    def base_cells(self):
        return np.asarray(ids(self.verified(self.data / "cells.tsv"), "cell_id"))

    def matrix_file(self, folder, spec):
        if spec["orientation"] != "cells_by_genes":
            raise ValueError(f"Unsupported matrix orientation in {folder}")
        expected = ("counts.npz" if self.layer == "counts" else
                    "expression_tpm.npz" if self.name == "villani" else "expression_normalized.npz")
        if spec["matrix"] != expected:
            raise ValueError(f"Expected {expected} in {folder}/input.json")
        return folder / expected

    @cached_property
    def base_matrix(self):
        path = self.verified(self.matrix_file(self.data, self.spec))
        matrix = sparse.load_npz(path).tocsr()
        if matrix.shape != (len(self.base_cells), len(self.genes)):
            raise ValueError(f"Full expression axes disagree in {self.data}")
        return matrix

    def source_files(self, view):
        data = self.data
        folder = data / "views" / view
        common = [data / "input.json", data / "views.tsv", data / "genes.tsv", folder / "cells.tsv"]
        if self.name == "cellsius":
            return common + [folder / "input.json", folder / "counts.npz"]
        return common + [self.matrix_file(data, self.spec), data / "cells.tsv", folder / "row_indices.npy"]

    def write(self, view, inputs):
        folder = self.data / "views" / view
        cells = ids(self.verified(folder / "cells.tsv"), "cell_id")
        if self.name == "cellsius":
            spec = read_json(self.verified(folder / "input.json"))
            if spec.get("view_id", view) != view:
                raise ValueError(f"CellSIUS input belongs to a different view in {folder}")
            matrix = sparse.load_npz(self.verified(self.matrix_file(folder, spec))).tocsr()
        else:
            idx = np.load(self.verified(folder / "row_indices.npy"), allow_pickle=False)
            if (idx.ndim != 1 or idx.dtype.kind not in "iu" or not len(idx) or idx.min() < 0
                    or idx.max() >= len(self.base_cells) or self.base_cells[idx].tolist() != cells):
                raise ValueError(f"RNA row mapping differs from cell IDs in {view}")
            full = self.base_matrix
            matrix = full if len(idx) == full.shape[0] and np.array_equal(idx, np.arange(len(idx))) else full[idx]
        if (matrix.shape != (len(cells), len(self.genes)) or len(cells) != self.sizes[view]
                or matrix.dtype.kind not in "iuf" or not np.isfinite(matrix.data).all() or (matrix.data < 0).any()):
            raise ValueError(f"Invalid expression or dimensions in {view}")
        matrix.check_format(full_check=True)
        inputs.mkdir(parents=True)
        sparse.save_npz(inputs / "expression.npz", matrix)
        write_tsv(inputs / "cells.tsv", ["cell_id"], ({"cell_id": c} for c in cells))
        write_tsv(inputs / "genes.tsv", ["gene_id", "symbol"],
                  ({"gene_id": g["feature_id"], "symbol": g["gene_symbol"]} for g in self.genes))
        cfg = dict(counts="expression.npz", cells="cells.tsv", genes="genes.tsv", orientation="cells_by_genes",
                   layer=self.layer, species=self.species, tissue=self.tissue, allowed_metadata=[])
        (inputs / "dataset.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf8")
        return cells


def copied_inputs(data, out, name, views=None, all_views=False):
    """Publish each requested view atomically. Other views can be added later."""
    source = ExpressionSource(data, name, lambda path: path)
    selected = views or (list(source.sizes) if all_views else default_views(name, source.entries))
    if set(selected) - source.sizes.keys():
        raise ValueError(f"Unknown view for {name}. Available: {', '.join(source.sizes)}")
    hashes = {}
    for v in dict.fromkeys(selected):
        dest = out / name / v
        files = source.source_files(v)
        for path in files:
            if path not in hashes:
                hashes[path] = hash_file(path)
        identity = {"format_version": 1, "dataset": name, "view": v, "source_sha256":
                    {str(p.relative_to(data)): hashes[p] for p in files}}
        if dest.exists():
            receipt = read_json(dest / "import.json")
            if receipt["source"] != identity or set(receipt["outputs"]) != set(INPUT_FILES):
                raise ValueError(f"Source changed or inventory incomplete at {dest}; choose a fresh --out")
            for filename in INPUT_FILES:
                if hash_file(dest / "input" / filename) != receipt["outputs"][filename]:
                    raise ValueError(f"Changed exported file: {dest / 'input' / filename}")
        else:
            dest.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(dir=dest.parent, prefix=".import-") as tmp:
                stage = Path(tmp) / "view"
                source.write(v, stage / "input")
                write_json(stage / "import.json", {"source": identity, "ground_truth": "not_imported",
                           "outputs": {f: hash_file(stage / "input" / f) for f in INPUT_FILES}})
                os.rename(stage, dest)
        yield dest / "input/dataset.yaml"


def write_references(stage, axes, path, outputs, *, membership=False, expected_pairs=None, mait=False):
    """Read a large legacy table once and route its rows to all requested views."""
    fields = ["cell_id", "target_id", "reference_kind", "state"] if membership else ["cell_id", "label"]
    coverage = {v: {} for v in axes}
    positive_mait = dict.fromkeys(axes, 0)
    with ExitStack() as stack:
        writers = {}
        for v in axes:
            directory = stage / v / "evaluation_fixture"
            directory.mkdir(exist_ok=True)
            for profile, filename in outputs.items():
                stream = stack.enter_context((directory / filename).open("w", encoding="utf8", newline=""))
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(fields)
                writers[v, profile] = writer

        for r in rows(path):
            v = r["view_id"]
            if v not in axes:
                continue
            profile = "membership" if membership else r["level_id"]
            if profile not in outputs:
                continue
            cell = r["cell_id"]
            axis = axes[v]
            if cell not in axis:
                raise ValueError(f"Reference cell outside RNA view {v}: {cell}")
            if membership:
                target, kind, state = r["target_id"], r["reference_kind"], r["state"]
                key = (target, kind)
                if not target or not kind or state not in {"positive", "negative", "unknown", "disputed"}:
                    raise ValueError(f"Invalid target membership in {v}")
                record = [cell, target, kind, state]
                if mait and target == "MAIT" and state == "positive":
                    positive_mait[v] += 1
            else:
                key = profile
                if not r["original_name"]:
                    raise ValueError(f"Empty label in {v}/{profile}")
                record = [cell, r["original_name"]]
            seen = coverage[v]
            if key not in seen:
                seen[key] = bytearray(len(axis))
            idx = axis[cell]
            if seen[key][idx]:
                raise ValueError(f"Duplicate ground truth for {cell} in {v}/{key}")
            seen[key][idx] = 1
            writers[v, profile].writerow(record)

    for v, seen in coverage.items():
        required = expected_pairs if membership else set(outputs)
        if not seen or (required is not None and seen.keys() != required) or any(0 in mask for mask in seen.values()):
            raise ValueError(f"Incomplete ground truth in {v}: {path.name}")
        if mait and v in MAIT_SIZES and positive_mait[v] != MAIT_SIZES[v]:
            raise ValueError(f"MAIT count differs from frozen recipe in {v}")


def export_release(project, out, name, all_views=False, *, data_root=None, benchmark_root=None):
    name = ALIASES.get(name, name)
    out = Path(out).resolve()
    if name not in DATASETS:
        raise ValueError(f"Unknown dataset {name}")
    data = (Path(data_root) if data_root else Path(project) / "data") / name
    bench = (Path(benchmark_root) if benchmark_root else Path(project) / "benchmark") / name
    lockfile = bench / "release.lock.json"
    lock = read_json(lockfile)["artifacts"]
    checked = set()

    def verified(path):
        if path.is_relative_to(data):
            rel = str(Path("data") / name / path.relative_to(data))
        else:
            rel = str(Path("benchmark") / name / path.relative_to(bench))
        if path not in checked:
            if rel not in lock or hash_file(path) != lock[rel]:
                raise ValueError(f"Single-agent release checksum mismatch: {rel}")
            checked.add(path)
        return path

    source = ExpressionSource(data, name, verified)
    views = default_views(name, source.entries)
    if all_views:
        views = list(source.sizes)
    manifest = read_json(verified(bench / "manifest.json"))
    identity = {"format_version": 1, "dataset": name, "release": manifest["gt_release_id"],
                "source_lock_sha256": hash_file(lockfile), "views": views}
    scoped = name in {"larry", "pbmc_tcr"} or (name == "cellsius_spheroid" and
            read_json(verified(bench / "source_lock.json")).get("evaluation_profile") == "spheroid_cp_target_f1_v1")
    label_files = {} if scoped else {"fine": "gt.tsv"}
    if name == "hao":
        label_files = {"l1": "labels_l1.tsv", "l2": "gt.tsv", "l3": "labels_l3.tsv"}
    membership_file = "gt.tsv" if scoped else "targets.tsv" if name in {"hao", "mmochi"} else None
    labels = verified(bench / "labels_source.tsv") if label_files else None
    memberships = verified(bench / "membership.tsv") if membership_file else None
    pairs = None
    if name in {"hao", "larry", "mmochi"}:
        pairs = {(r["target_id"], r["reference_kind"]) for r in rows(verified(bench / "targets.tsv"))}
        if not pairs:
            raise ValueError(f"Empty target registry for {name}")
    dest = out / name
    if dest.exists():
        saved = read_json(dest / "export.json")
        if saved["source"] != identity:
            raise ValueError(f"Different release or view selection already at {dest}; choose another --out")
        references = [*label_files.values(), *([membership_file] if membership_file else [])]
        required = {f"{v}/evaluation_fixture/{f}" for v in views for f in references}
        required.update(f"{v}/input/{f}" for v in views for f in INPUT_FILES)
        if set(saved["files"]) != required:
            raise ValueError(f"Incomplete export inventory at {dest}")
        for v in views:
            for path in source.source_files(v):
                verified(path)
            verified(bench / "views" / v / "cells.tsv")
        for rel, digest in saved["files"].items():
            if hash_file(dest / rel) != digest:
                raise ValueError(f"Changed exported file: {dest / rel}")
        return dest

    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out, prefix="." + name + "-") as tmp:
        stage = Path(tmp) / name
        axes = {}
        for v in views:
            cells = source.write(v, stage / v / "input")
            if cells != ids(verified(bench / "views" / v / "cells.tsv"), "cell_id"):
                raise ValueError(f"Benchmark cell list differs in {v}")
            axes[v] = {cell: i for i, cell in enumerate(cells)}
        if labels:
            write_references(stage, axes, labels, label_files)
        if memberships:
            write_references(stage, axes, memberships, {"membership": membership_file},
                             membership=True, expected_pairs=pairs, mait=name == "pbmc_tcr")
        files = {str(path.relative_to(stage)): hash_file(path) for path in sorted(stage.rglob("*")) if path.is_file()}
        write_json(stage / "export.json", {"source": identity, "files": files})
        os.rename(stage, dest)
    return dest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    roots = parser.add_mutually_exclusive_group(required=True)
    roots.add_argument("--project", type=Path, help="Frozen project containing data/ and benchmark/")
    roots.add_argument("--data-root", type=Path, help="Copied directory containing airway/, cellsius/, etc.")
    parser.add_argument("--benchmark-root", type=Path, help="With --data-root, verify against this benchmark/ and export GT")
    parser.add_argument("--out", type=Path, help="Defaults to work/copied-inputs or work/single-agent-datasets (verified GT)")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--all-views", action="store_true", help="Include optional prepared views")
    selection.add_argument("--view", action="append", help="Import just this view (repeatable, expression-only mode)")
    parser.add_argument("--list", action="store_true", help="List view IDs, roles and cell counts without reading matrices")
    parser.add_argument("datasets", nargs="+", choices=sorted(set(DATASETS) | set(ALIASES) | {"all"}))
    args = parser.parse_args(argv)
    locked = args.project is not None or args.benchmark_root is not None
    if args.project and args.benchmark_root:
        parser.error("--benchmark-root is used with --data-root")
    if args.view and (locked or args.list):
        parser.error("--view is for expression-only import with --data-root")
    data_root = (args.data_root or args.project / "data").resolve()
    bench_root = (args.benchmark_root or args.project / "benchmark").resolve() if locked else None
    out = (args.out or (RELEASE_OUT if locked else COPIED_OUT)).resolve()
    failed = False
    if "all" in args.datasets and len(args.datasets) != 1:
        parser.error("Use all by itself or specify individual dataset names")
    names = DATASETS if args.datasets == ["all"] else dict.fromkeys(ALIASES.get(n, n) for n in args.datasets)
    for name in names:
        try:
            if args.list:
                entries, sizes = view_sizes(data_root / name, lambda path: path)
                for r in entries:
                    v = r["view_id"]
                    print(f"{name}\t{v}\t{r['role']}\t{sizes[v]}", flush=True)
                continue
            for root in (data_root, bench_root):
                if root and (out.is_relative_to(root) or (root / name).is_relative_to(out / name)):
                    raise ValueError("Choose an output directory separate from the source data and benchmark")
            if locked:
                print(f"Verifying and exporting {name} ...", file=sys.stderr, flush=True)
                dest = export_release(None, out, name, args.all_views, data_root=data_root, benchmark_root=bench_root)
                for v in read_json(dest / "export.json")["source"]["views"]:
                    print(f"{dest/v/'input/dataset.yaml'}\t{dest/v/'evaluation_fixture/gt.tsv'}", flush=True)
            else:
                for path in copied_inputs(data_root / name, out, name, args.view, args.all_views):
                    print(path, flush=True)
        except (OSError, ValueError, KeyError) as exc:
            print(f"{name}: {exc}", file=sys.stderr, flush=True)
            failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
