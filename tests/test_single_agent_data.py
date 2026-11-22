"""Contract tests for copied inputs and for locked releases with split roots."""
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from rarecell.config import hash_file
from rarecell.dashboard import imported_views
from rarecell.data import read_dataset, write_tsv
from rarecell.single_agent_data import main


def legacy_tree(root, name="villani"):
    data = root / "source" / name
    data.mkdir(parents=True)
    shared = name != "cellsius"
    mat = sparse.csr_matrix([[1.25, 0], [3.5, 4], [0, 5.5]]) if shared else sparse.eye(2, format="csr")
    (data / "genes.tsv").write_text("feature_id\tgene_symbol\ng1\tS\ng2\tS\n")
    spec = {"matrix": "expression_tpm.npz", "orientation": "cells_by_genes", "is_raw_count": False}
    if shared:
        sparse.save_npz(data / spec["matrix"], mat)
        (data / "cells.tsv").write_text("cell_id\nc0\nc1\nc2\n")
    (data / "input.json").write_text(json.dumps(spec if shared else {"is_raw_count": True}))
    entries = []
    for v in (["V000"] if shared else [f"V00{i}" for i in range(6)]):
        folder = data / "views" / v
        folder.mkdir(parents=True)
        (folder / "cells.tsv").write_text("cell_id\tlabel\nc2\tprivate\nc0\tprivate\n")
        if shared:
            np.save(folder / "row_indices.npy", np.array([2, 0]))
        else:
            sparse.save_npz(folder / "counts.npz", mat * (int(v[-1]) + 1))
            (folder / "input.json").write_text(json.dumps({"matrix": "counts.npz", "orientation": "cells_by_genes", "view_id": v}))
        entries.append({"view_id": v, "role": "main", "n_cells": 2})
    write_tsv(data / "views.tsv", ["view_id", "role", "n_cells"], entries)
    return data


def test_incremental_per_view_and_integrity(tmp_path):
    data = legacy_tree(tmp_path, "cellsius")
    out = tmp_path / "work/copied-inputs"
    args = ["--data-root", str(data.parent), "--out", str(out), "cellsius"]
    assert main(args + ["--view", "V001"]) == 0
    item, = imported_views(tmp_path)
    assert item["gt"] == ""
    ds = read_dataset(item["dataset"])
    assert ds.cell_ids == ["c0", "c2"] and ds.raw.sum() == 4
    assert set(ds.cells[0]) == {"cell_id"}
    assert [g["symbol"] for g in ds.genes] == ["S", "S"]
    assert not list(out.rglob("gt.tsv"))
    assert main(args) == 0
    assert len(imported_views(tmp_path)) == 6
    (out / "cellsius/V001/input/cells.tsv").write_text("changed")
    assert main(args + ["--view", "V001"]) == 1


def test_indexed_rna_and_change_detection(tmp_path, capsys):
    data = legacy_tree(tmp_path)
    out = tmp_path / "out"
    args = ["--data-root", str(data.parent), "--out", str(out), "villani"]
    assert main(args + ["--list"]) == 0
    assert capsys.readouterr().out.strip() == "villani\tV000\tmain\t2"
    assert not out.exists()
    assert main(args) == 0
    folder = out / "villani/V000/input"
    loaded = read_dataset(folder)
    assert loaded.raw is None and loaded.info["layer"] == "linear_expression"
    np.testing.assert_array_equal(sparse.load_npz(folder / "expression.npz").toarray(), [[0, 5.5], [1.25, 0]])
    np.save(data / "views/V000/row_indices.npy", [0, 2])
    assert main(args) == 1
    assert "Source changed" in capsys.readouterr().err
    fresh = tmp_path / "fresh"
    assert main(["--data-root", str(data.parent), "--out", str(fresh), "villani"]) == 1
    assert "RNA row mapping" in capsys.readouterr().err
    assert not list(fresh.rglob("dataset.yaml"))


def test_locked_release_with_separate_roots(tmp_path):
    data = legacy_tree(tmp_path)
    bench = tmp_path / "old-project/benchmark/villani"
    view = bench / "views/V000"
    view.mkdir(parents=True)
    (view / "cells.tsv").write_bytes((data / "views/V000/cells.tsv").read_bytes())
    (bench / "labels_source.tsv").write_text("view_id\tcell_id\tlevel_id\toriginal_name\nV000\tc2\tfine\tA\nV000\tc0\tfine\tB\n")
    (bench / "manifest.json").write_text('{"gt_release_id":"test"}')
    hashes = {str(prefix / p.relative_to(src.parent)): hash_file(p)
              for prefix, src in ((Path("data"), data), (Path("benchmark"), bench))
              for p in src.rglob("*") if p.is_file()}
    (bench / "release.lock.json").write_text(json.dumps({"artifacts": hashes}))
    out = tmp_path / "work/single-agent-datasets"
    args = ["--data-root", str(data.parent), "--benchmark-root", str(bench.parent), "--out", str(out), "villani"]
    assert main(args) == 0 and main(args) == 0
    item, = imported_views(tmp_path)
    assert item["gt"] == str(out / "villani/V000/evaluation_fixture/gt.tsv")
    assert "private" not in (out / "villani/V000/input/cells.tsv").read_text()
    (data / "genes.tsv").write_text("other contents")
    assert main(args) == 1
