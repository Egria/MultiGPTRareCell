"""Numerical controls for local-contrast/1. No provider or annotation dependency."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.preprocessing import normalize
from sklearn.random_projection import SparseRandomProjection
from threadpoolctl import threadpool_limits

from rarecell.data import read_dataset
from rarecell.science import Science


def write_table(path, header, values):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow([header])
        writer.writerows(([x] for x in values))


def expression_fixture(folder, rare_size, technical_scale=False, homogeneous=False, shuffled=False):
    """Create real CSR count input; cell identities carry no annotation or target hints."""
    folder.mkdir()
    dataset = folder / "discovery_input"
    dataset.mkdir()
    n = 80
    rng = np.random.RandomState(461)
    cell_ids = ["BC" + str(i).zfill(3) for i in rng.permutation(n)]
    genes = [f"g{i:03}" for i in range(144)]
    counts = np.zeros((n, len(genes)))
    counts[:, 6:14] = 20
    counts[:40, :3] = 100
    counts[40:, 3:6] = 100
    if homogeneous:
        counts[40:] = counts[0]
    rare_indices = [j // 2 if j % 2 == 0 else 40 + j // 2 for j in range(rare_size)]
    counts[np.asarray(rare_indices, dtype=int), 20:25] = 75
    if technical_scale:
        counts = counts * (1 + np.arange(n) % 6)[:, None]
    order = rng.permutation(n) if shuffled else np.arange(n)
    write_table(dataset / "cells.tsv", "cell_id", [cell_ids[i] for i in order])
    write_table(dataset / "genes.tsv", "gene_id", genes)
    sparse.save_npz(dataset / "counts.npz", sparse.csr_matrix(counts[order]))
    spec = {"counts": "counts.npz", "cells": "cells.tsv", "genes": "genes.tsv",
            "layer": "counts", "orientation": "cells_by_genes", "species": "synthetic", "tissue": "synthetic"}
    (dataset / "dataset.yaml").write_text(yaml.safe_dump(spec))
    # Only the evaluation-side fixture has these truth assignments.
    evaluation = folder / "evaluation_fixture"
    evaluation.mkdir()
    truth = {"seed": cell_ids[0], "rare": [cell_ids[i] for i in rare_indices],
             "coarse": {cell_ids[i]: i // 40 for i in range(n)}}
    (evaluation / "truth.json").write_text(json.dumps(truth))

    data = read_dataset(dataset)
    snap = folder / "snapshot"
    snap.mkdir()
    (snap / "views").mkdir()
    rows = np.arange(n) if not shuffled else rng.permutation(n)
    write_table(snap / "cells.tsv", "cell_id", [data.cell_ids[i] for i in rows])
    write_table(snap / "genes.tsv", "gene_id", data.gene_ids)
    x = data.x[rows]
    sparse.save_npz(snap / "data.npz", x)
    sparse.save_npz(snap / "raw.npz", data.raw[rows])
    rp = SparseRandomProjection(n_components=96, random_state=373)
    view = rp.fit_transform(normalize(x)).toarray()
    np.save(snap / "views" / "V_LINEAR.npy", view)
    (snap / "manifest.json").write_text(json.dumps({"views": ["V_LINEAR"]}))
    return snap


class NumericalSearchTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.root = Path(tmp.name)
        self.addCleanup(tmp.cleanup)

    def test_named_contrast_returns_depth_adjusted_high_members(self):
        snap = expression_fixture(self.root / "case", 6, technical_scale=True, shuffled=True)
        truth = json.loads((snap.parent / "evaluation_fixture" / "truth.json").read_text())
        tool = Science(snap, {})
        args = {"members": [truth["seed"]], "genes": ["g020"]}
        result = tool.execute("compare_groups", args)
        self.assertEqual(set(result["nb_high_members"]["g020"]), set(truth["rare"]))
        gene = result["summary"]["nb_mixture"]["genes"]["g020"]
        self.assertEqual((gene["high_foreground"], gene["high_background"]), (1, 5))
        (snap / "raw.npz").unlink()
        fallback = Science(snap, {}).execute("compare_groups", args)
        self.assertEqual(fallback["statistics"], result["statistics"])
        self.assertEqual(fallback["summary"]["nb_mixture"]["status"], "raw_counts_unavailable")

    def run_search(self, rare_size, name="run", search=None, **kwargs):
        snap = expression_fixture(self.root / name, rare_size, **kwargs)
        truth = json.loads((snap.parent / "evaluation_fixture" / "truth.json").read_text())
        tool = Science(snap, {"search": search or {}})
        # The input to the algorithm is a single seed. It cannot access evaluation_fixture.
        result = tool.execute("expand_candidate", {"members": [truth["seed"]]})
        return tool, result, truth

    def test_singleton_recovers_across_coarse_classes(self):
        tool, result, truth = self.run_search(4)
        self.assertEqual(result["status"], "expanded")
        self.assertEqual(set(result["proposed_members"]), set(truth["rare"]))
        self.assertEqual({truth["coarse"][c] for c in result["proposed_members"]}, {0, 1})
        # Check a computed blind coarse partition as well as the generating strata.
        with threadpool_limits(limits=1):
            coarse = KMeans(n_clusters=2, random_state=461, n_init=2).fit_predict(tool.data)
        self.assertEqual(sorted(np.bincount(coarse).tolist()), [40, 40])
        self.assertEqual(len({coarse[tool.cell_index[c]] for c in result["proposed_members"]}), 2)
        history = result["selection_history"]
        self.assertTrue(set(result["added_members"]) - set(history["background"]))
        self.assertEqual(history["dependency"], "selection-dependent")
        self.assertTrue({g["gene_id"] for g in history["signature"]} >= {"g020", "g021"})
        self.assertEqual(len(history["background"]), 24)
        self.assertTrue(all(x["expression"] >= .65 and x["score"] >= .72 for x in result["scores"]))
        self.assertEqual(result["signature_kind"], "exploratory_singleton")

    def test_distinct_population_sizes_and_isolated_seed(self):
        sizes = []
        for n in (1, 2, 5, 8):
            _, result, truth = self.run_search(n, name=f"n{n}")
            self.assertEqual(set(result["proposed_members"]), set(truth["rare"]))
            sizes.append(len(result["proposed_members"]))
            if n == 1:
                self.assertEqual(result["status"], "no_expansion")
        self.assertEqual(sizes, [1, 2, 5, 8])

    def test_zero_signal_and_technical_library_size(self):
        tool, result, truth = self.run_search(0, technical_scale=True, homogeneous=True)
        self.assertEqual(result["status"], "no_signature")
        self.assertEqual(result["proposed_members"], [truth["seed"]])
        self.assertTrue(np.allclose(tool.data.toarray(), tool.data.toarray()[0]))
        self.assertEqual(result["added_members"], [])

    def test_iid_count_noise_has_no_forced_expansion(self):
        snap = expression_fixture(self.root / "null", 0)
        dataset = snap.parent / "discovery_input"
        sizes = []
        for rng_seed in (21, 22, 23, 24, 25):
            rng = np.random.RandomState(rng_seed)
            # All cells share gene-specific rates with independent count noise.
            rate = rng.gamma(shape=2.0, scale=1.0, size=144)
            counts = sparse.csr_matrix(rng.poisson(rate, size=(80, 144)))
            sparse.save_npz(dataset / "counts.npz", counts)
            data = read_dataset(dataset)
            sparse.save_npz(snap / "data.npz", data.x)
            projected = SparseRandomProjection(n_components=96, random_state=373).fit_transform(normalize(data.x))
            np.save(snap / "views" / "V_LINEAR.npy", projected.toarray())
            result = Science(snap, {}).execute("expand_candidate", {"members": [data.cell_ids[0]]})
            sizes.append(len(result["added_members"]))
        self.assertEqual(sizes, [0, 0, 0, 0, 0])

    def test_barcode_ties_after_reordering_and_depth_normalization(self):
        config = {"max_cells": 2, "block_size": 7}
        _, first, _ = self.run_search(8, "first", config)
        _, second, _ = self.run_search(8, "second", config, technical_scale=True, shuffled=True)
        for result in [first, second]:
            self.assertEqual(result["status"], "incomplete")
            self.assertEqual(result["summary"]["n_passing"], 7)
        self.assertEqual(first["added_members"], second["added_members"])
        self.assertEqual(first["added_members"], sorted(first["added_members"]))
        self.assertEqual(first["selection_history"]["background"], second["selection_history"]["background"])

    def test_original_barcode_controls_ties_across_sample_prefixed_ids(self):
        tool, full, truth = self.run_search(8)
        ids = tool.cells.tolist()
        remaining = sorted(full["scores"], key=lambda row: -round(row["score"], 12))
        top_score = round(remaining[0]["score"], 12)
        tied = [row["cell_id"] for row in remaining if round(row["score"], 12) == top_score]
        self.assertGreaterEqual(len(tied), 2)
        # Reverse the first tie by assigning different original barcodes.
        alternate = {cell: f"TAG{100-i:03}" for i, cell in enumerate(sorted(ids))}
        with (tool.snapshot / "cells.tsv").open("w") as f:
            f.write("cell_id\tbarcode\n")
            for cell in ids:
                f.write(f"{cell}\t{alternate[cell]}\n")
        capped = Science(tool.snapshot, {"search": {"max_cells": 1}})
        result = capped.execute("expand_candidate", {"members": [truth["seed"]]})
        self.assertEqual(result["added_members"], [min(tied, key=lambda c: (alternate[c], c))])

    def test_expression_scope_direction_and_logic(self):
        tool, result, truth = self.run_search(4)
        members = [truth["seed"], next(c for c in tool.cells if c not in truth["rare"])]
        background = [c for c in tool.cells if c not in truth["rare"]][:2]
        compare = tool.execute("compare_groups", {"members": [truth["seed"]], "background": background,
                                                  "genes": ["g020", "g006"]})
        records = {r["gene_id"]: r for r in compare["statistics"]}
        self.assertEqual(records["g020"]["fraction_query"], 1)
        self.assertEqual(records["g020"]["fraction_background"], 0)
        self.assertEqual(records["g020"]["direction"], "up")
        self.assertEqual(records["g006"]["direction"], "down")
        args = {"members": members, "genes": ["g020", "g006"], "threshold": 0}
        conjunction = tool.execute("coexpress", dict(args, logic="AND"))
        disjunction = tool.execute("coexpress", dict(args, logic="OR"))
        self.assertEqual(conjunction["matched_members"], [truth["seed"]])
        self.assertEqual(disjunction["matched_members"], sorted(members))
        self.assertEqual(disjunction["gene_counts"], {"g006": 2, "g020": 1})
        higher = tool.execute("coexpress", dict(args, threshold=20))
        self.assertEqual(higher["summary"]["n_matching"], 0)
        self.assertEqual(tool.execute("retrieve_cells", {"members": [truth["seed"]]}), result)

    def test_absence_only_signature_has_no_positive_subscore(self):
        tool, _, truth = self.run_search(4)
        seed = truth["seed"]
        # Use an already-normalized layer to engineer a strictly negative
        # contrast while retaining nonzero RNA in every cell.
        a = sorted(truth["rare"])
        b = [c for c in tool.cells if c not in a]
        values = np.zeros(tool.data.shape)
        values[:, 0] = 2
        values[:, 1] = 5
        for c in a:
            values[tool.cell_index[c], 1] = 0
        layer = sparse.csr_matrix(values)
        sparse.save_npz(tool.snapshot / "data.npz", layer)
        transform = SparseRandomProjection(n_components=96, random_state=373)
        np.save(tool.snapshot / "views" / "V_LINEAR.npy", transform.fit_transform(normalize(layer)).toarray())
        new = Science(tool.snapshot, {})
        out = new.execute("expand_candidate", {"members": [seed], "background": b})
        self.assertTrue(out["added_members"])
        self.assertTrue(all(s["direction"] == "down" for s in out["selection_history"]["signature"]))
        self.assertTrue(all(row["positive_expression"] is None for row in out["scores"]))

    def test_zero_vector_small_background_and_invalid_entities(self):
        snap = expression_fixture(self.root / "edge", 1)
        truth = json.loads((snap.parent / "evaluation_fixture" / "truth.json").read_text())
        seed = truth["seed"]
        tool = Science(snap, {})
        target = next(c for c in tool.cells if c != seed)
        args = {"members": [seed], "background": [target]}
        self.assertEqual(tool.execute("expand_candidate", args)["status"], "insufficient_background")
        for bad in [dict(args, background=[seed]), dict(args, background=[]), {"members": ["almost" + seed]},
                    dict(args, min_score=0), dict(args, layer=None)]:
            with self.assertRaises(ValueError):
                tool.execute("expand_candidate", bad)
        for bad in [None, -1, True, float("inf")]:
            with self.assertRaises(ValueError):
                tool.execute("coexpress", {"members": [seed], "genes": ["g000"], "threshold": bad})
        with self.assertRaises(ValueError):
            tool.execute("coexpress", {"members": [seed], "genes": ["g00o"]})
        csr = sparse.load_npz(snap / "data.npz")
        row = tool.cell_index[seed]
        csr.data[csr.indptr[row]:csr.indptr[row + 1]] = 0
        sparse.save_npz(snap / "data.npz", csr)
        tool = Science(snap, {})
        result = tool.execute("expand_candidate", {"members": [seed]})
        self.assertEqual(result["status"], "no_valid_structure")
        self.assertEqual(result["proposed_members"], [seed])


if __name__ == "__main__":
    unittest.main()
