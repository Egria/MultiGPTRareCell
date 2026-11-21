import csv
import hashlib
import json

import numpy as np
import pytest
import yaml

from rarecell.benchmark_views import generate_abundance_series
from rarecell.data import read_dataset
from rarecell.evaluation import aggregate, evaluate, _partition_scores
from rarecell.synthetic import generate, write_tsv


def dump(path, obj):
    path.write_text(json.dumps(obj), encoding="utf-8")


def freeze(run):
    files = {p.relative_to(run).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in run.rglob("*") if p.is_file() and p.name != "frozen.json"}
    dump(run / "frozen.json", {"files": files, "snapshot_hash": "synthetic-example"})


def scoring_fixture(root):
    # Global assignment forbids reuse: g is a perfect-separation competitor for both A and B;
    # A also has groups h/k, whereas B can only be matched to g.
    run, gt = root / "run", root / "gt.tsv"
    (run / "artifacts").mkdir(parents=True)
    targets = ["A"] * 5 + ["B"] * 3 + ["C"] * 2
    groups = ["g", "g", "g", "h", "k", "g", "g", "g", "x", "x"]
    cells = ["s%02d" % i for i in range(10)]
    owners = {"g": "q", "h": "a", "k": "", "x": "c"}
    part = [{"cell_id": cell, "group_id": g, "candidate_id": owners[g],
             "status": "unresolved" if g == "k" else "provisional"}
            for cell, g in zip(cells, groups)]
    write_tsv(run / "partition.tsv", part, list(part[0]))
    write_tsv(run / "artifacts/baselines.tsv",
              [dict(r, view="V_LINEAR", resolution="0.5") for r in part],
              ["view", "resolution", "cell_id", "group_id"])
    write_tsv(gt, ({"cell_id": c, "label": t} for c, t in zip(cells, targets)), ["cell_id", "label"])
    candidates = {"a": [cells[3]], "b": [cells[6], cells[7]], "c": cells[8:],
                  "q": [c for c, g in zip(cells, groups) if g == "g"]}
    with (run / "candidates.jsonl").open("w") as f:
        for cid in candidates:
            f.write(json.dumps({"id": cid, "initial": cid != "q", "parents": ["b"] if cid == "q" else []})+"\n")
    write_tsv(run / "candidate_members.tsv", ({"candidate_id": cid, "cell_id": cell}
              for cid, mem in candidates.items() for cell in mem), ["candidate_id", "cell_id"])
    dump(run / "usage.json", {"new_queries": 3})
    dump(run / "run_summary.json", {"status": "partial"})
    freeze(run)
    return run, gt


def test_global_matching_diagnostics_and_freeze(tmp_path):
    run, gt = scoring_fixture(tmp_path)
    out = tmp_path / "evaluation"
    result = evaluate(run, gt, out)
    per = {r["target"]: r for r in result["submitted_partition"]["per_target"]}
    assert per["B"]["group_id"] == "g"
    assert per["B"]["recall"] == 1 and per["B"]["precision"] == .5
    assert per["B"]["f1"] == pytest.approx(2/3)
    assert per["A"]["fragments"] == 3 and per["A"]["tp"] == 1
    assert per["A"]["f1"] == pytest.approx(1/3)
    assert per["C"]["f1"] == 1
    assert result["n_candidates"] == 4 and result["n_initial_candidates"] == 3
    assert result["unresolved_fraction"] == .1
    diag = {r["target"]: r for r in result["offline_single_candidate_diagnostics"]}
    assert diag["B"]["best_single_candidate_f1"] == .8
    assert diag["B"]["new_true_members_vs_initial"] == 1
    assert result["baselines"][0]["per_target"] == result["submitted_partition"]["per_target"]
    assert result["usage"]["new_queries"] == 3 and result["run_summary"]["status"] == "partial"
    with (out / "revision_changes.tsv").open() as f:
        assert len(list(csv.DictReader(f,delimiter="\t"))) == 4
    assert evaluate(run, gt, out) == result
    with pytest.raises(ValueError):
        evaluate(run, gt, run / "evaluation")
    (run / "partition.tsv").write_text((run / "partition.tsv").read_text() + "\n")
    with pytest.raises(ValueError, match="digest"):
        evaluate(run, gt, out)


def test_unmatched_tie_and_no_oracle_union():
    gt = {"00": "A", "01": "B", "02": "C"}
    pred = {c: "single" for c in gt}
    first = _partition_scores(gt, pred)
    permuted = _partition_scores(dict(reversed(list(gt.items()))), pred)
    assert first == permuted
    assert len(first["unmatched_targets"]) == 2
    assert first["macro_f1"] == pytest.approx(.5 / 3)
    assert first["mixed_groups"] == 1
    assert _partition_scores(gt, gt)["ari"] == 1
    with pytest.raises(ValueError, match="differ"):
        _partition_scores(gt, {"00": "A"})


def test_freeze_before_reading_gt_and_plan_missing(tmp_path):
    run, gt = scoring_fixture(tmp_path)
    out = tmp_path / "scored"
    plan = tmp_path / "plan.yaml"
    plan.write_text(yaml.safe_dump({"replications":[
        {"arm":"A2", "replication_id":"s1","metrics":"scored/metrics.json"},
        {"arm":"A2", "replication_id":"s2","metrics":"missing/metrics.json"}]}))
    evaluate(run, gt, out)
    dist = aggregate(plan, tmp_path / "distribution.json")
    assert dist["n_planned"] == 2 and dist["n_scored"] == 1
    assert dist["replications"][1]["status"] == "missing"
    assert dist["replications"][0]["run_status"] == "partial"
    (run / "frozen.json").unlink()
    with pytest.raises(FileNotFoundError, match="frozen.json"):
        evaluate(run, tmp_path / "no_gt", out)


def test_synthetic_and_abundance_series(tmp_path):
    dataset = generate(tmp_path / "src", seed=17)
    gt = tmp_path / "src" / "evaluation_fixture" / "gt.tsv"
    ds = read_dataset(dataset)
    assert ds.raw.shape == (160,96) and np.all(ds.raw.data >= 1)
    assert set(p.name for p in dataset.parent.iterdir()) == {"counts.npz","cells.tsv","genes.tsv","dataset.yaml"}
    series = generate_abundance_series(dataset, gt, tmp_path / "nested", "synthetic_rare", (8,4,2), seed=31)
    prev = None
    for i,n in enumerate((8,4,2)):
        view = series / "discovery" / ("view%d" % i)
        src = read_dataset(view / "dataset.yaml")
        assert len(src.cell_ids) == 152+n
        assert np.array_equal(src.raw.toarray(), ds.raw[[ds.cell_ids.index(c) for c in src.cell_ids]].toarray())
        assert all(set(r) == {"cell_id"} for r in src.cells)
        audit = series / "evaluation_fixture" / ("view%d" % i)
        with (audit / "retained_target.tsv").open() as f:
            current = {r["cell_id"] for r in csv.DictReader(f, delimiter="\t")}
        assert len(current) == n and (prev is None or current < prev)
        prev = current
    with pytest.raises(ValueError, match="fewer"):
        generate_abundance_series(dataset, gt, tmp_path/"invalid", "synthetic_rare", (144,72))


def test_frozen_escape_and_evaluation_destination(tmp_path):
    run, gt = scoring_fixture(tmp_path)
    dest = tmp_path / "outputs"
    dest.mkdir()
    (dest / "metrics.json").write_text('{"existing":"unrelated"}')
    with pytest.raises(ValueError, match="fresh output"):
        evaluate(run, gt, dest)
    assert json.loads((dest / "metrics.json").read_text()) == {"existing":"unrelated"}
    artifact = run / "usage.json"
    outside = tmp_path / "outside.json"
    artifact.rename(outside)
    artifact.symlink_to(outside)
    with pytest.raises(ValueError, match="outside discovery"):
        evaluate(run, gt, tmp_path / "scores")


def test_sparse_matching_matches_dense_reference():
    from collections import Counter
    from scipy.optimize import linear_sum_assignment

    rng = np.random.default_rng(2211)
    for _ in range(30):
        n = int(rng.integers(5, 55))
        truth = rng.integers(0, rng.integers(1, 7), n)
        pred = rng.integers(0, rng.integers(1, 9), n)
        score = _partition_scores(dict(enumerate(map(str, truth))), dict(enumerate(map(str, pred))))
        target_sizes, pred_sizes = Counter(truth), Counter(pred)
        exact = np.array([[2 * np.count_nonzero((truth == a) & (pred == b)) / (na + nb)
                           for b, nb in pred_sizes.items()] for a, na in target_sizes.items()])
        i, j = linear_sum_assignment(exact, maximize=True)
        assert score["macro_f1"] == pytest.approx(exact[i, j].sum() / len(target_sizes))
