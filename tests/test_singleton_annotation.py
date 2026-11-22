import json

import numpy as np
import pytest
from scipy import sparse

from rarecell import annotation, job_queue
from rarecell.base_evidence import annotation_parents
from rarecell.cli import parser
from rarecell.config import config, fingerprint, write_json
from rarecell.singleton_markers import panels, SOURCE


def snapshot(path):
    path.mkdir()
    (path / "seurat").mkdir()
    cells = ["a", "b", "c", "d", "e"]
    (path / "cells.tsv").write_text("cell_id\n" + "\n".join(cells) + "\n")
    (path / "genes.tsv").write_text("gene_id\tsymbol\ng0\tLow\ng1\tHigh\ng2\tHigh\ng3\tEqual\ng4\t\n")
    sparse.save_npz(path / "data.npz", sparse.csr_matrix(np.array([
        [0, 8, 7, 2, 2], [40, 0, 1, 2, 1], [20, 4, 3, 2, 1], [0, 0, 0, 0, 0], [0, 1, 0, 0, 0]], dtype=float)))
    local = {"type": "local_neighbours", "members": ["b", "c"]}
    rest = {"type": "complement"}
    seeds = [{"id": c, "candidate_id": "H" + c, "cell_id": c} for c in ["a", "d", "e"]]
    rows = [dict(id="H" + c, members=[c], sources=[{"kind": "outlier", "view": "V_GINI"}],
                 marker_background=bg, markers=[]) for c, bg in (("a", local), ("d", rest), ("e", rest))]
    rows[0]["sources"].append({"kind": "outlier", "view": "V_LINEAR"})
    rows.append(dict(id="Hp", members=["a", "b", "c"], sources=[
        dict(kind="partition", view="V_HVG", partition_id="P", resolution=0.4)]))
    for name, records in (("initial", rows), ("seeds", seeds)):
        (path / (name + ".jsonl")).write_text("".join(json.dumps(r) + "\n" for r in records))
    (path / "seurat/markers.tsv").write_text(
        "resolution\tgroup_id\tgene_id\tsymbol\tavg_log2FC\tp_val_adj\n0.4\tHp\tg1\tHigh\t2\t0.01\n")
    return dict(snapshot_hash="snap", data_version="data", views=["V_HVG"], stages={"seurat_markers": {"status": "ok"}})


def test_frozen_direction_background_symbols_and_ineligible(tmp_path):
    snap = tmp_path / "snapshot"
    snapshot(snap)
    a, d, e = list(panels(snap, 10))
    assert a["genes"] == ["High", "g4"]
    assert a["n_background"] == 2
    assert a["markers"][0] == dict(gene_id="g1", symbol="High", effect=6, pct_in=1, pct_out=.5)
    assert d["reason"] == "no_expression" and e["reason"] == "no_positive_markers"
    assert len(list(panels(snap, 1))[0]["genes"]) == 1


def cached_reply(out, job):
    # Test fixture at the persisted upstream receipt boundary; no external call.
    batch = out / "batches" / fingerprint(job)
    batch.mkdir(parents=True)
    write_json(batch / "response.json", dict(commit=annotation.PIN, valid=True, raw_response="unknown",
               rows=[dict(candidate_id=g["candidate_id"], label="unknown", ontology_id=None) for g in job["groups"]]))


def test_old_receipt_reuse_and_singleton_provenance(tmp_path):
    snap, old, new = [tmp_path / s for s in ("snap", "old", "new")]
    manifest = snapshot(snap)
    old_cfg = {"model": "example", "rscript": "/unavailable"}
    cfg = dict(old_cfg, include_singletons=True)
    partition, = annotation._jobs(snap, manifest, annotation._resolve(old_cfg))[0]
    jobs, skipped = annotation._jobs(snap, manifest, annotation._resolve(cfg))
    assert jobs[0] == partition and len(jobs) == 2 and len(skipped) == 2
    cached_reply(old, partition)
    assert annotation._run(snap, old, old_cfg, manifest)["status"] == "success"
    cached_reply(new, jobs[1])
    status = annotation._run(snap, new, cfg, manifest, reuse=old)
    assert status["status"] == "success" and status["cached_batches"] == 2
    assert status["annotated_singletons"] == 1 and status["newly_completed_batches"] == 0
    ann = annotation.verify_annotation(snap, new, _manifest=manifest)[1]
    candidate = dict(members=["a"], marker_background={"type":"local_neighbours","members":["b","c"]},
                     evidence_ids=["observed"])
    assert ann["marker_source"] == SOURCE and ann["selection_history"]["genes"] == ["High", "g4"]
    assert annotation_parents(ann, candidate) == ["observed"]
    with pytest.raises(ValueError, match="contrast"):
        annotation_parents(ann, dict(candidate, marker_background={"type":"complement"}))
    with pytest.raises(ValueError, match="settings must match"):
        annotation._reuse_receipts(snap, old, new, jobs, annotation._resolve(dict(cfg, model="another")), manifest)


def test_singleton_result_when_seurat_blocked(tmp_path):
    snap, out = tmp_path / "snap", tmp_path / "ann"
    m = snapshot(snap)
    m["views"] = ["V_GINI"]
    cfg = dict(model="example", rscript="/unavailable", include_singletons=True)
    jobs, _ = annotation._jobs(snap, m, annotation._resolve(cfg))
    assert len(jobs) == 1
    cached_reply(out, jobs[0])
    state = annotation._run(snap, out, cfg, m)
    assert state["status"] == "partial" and state["partition_issue"] and not state["errors"]
    assert len(annotation.verify_annotation(snap, out, allow_partial=True, _manifest=m)) == 1


def test_task_settings_and_cache_source(tmp_path):
    snap, old = tmp_path / "snap", tmp_path / "old"
    snap.mkdir()
    old.mkdir()
    root = tmp_path / "tasks"
    job_queue.initialize(root)
    form = dict(kind="annotate", snapshot=str(snap), reuse_annotations=str(old), annotation_model="example",
                include_singletons="yes")
    tid = job_queue.submit(root, tmp_path, form)
    cfg = config(root / tid / "settings.yaml")
    assert cfg["annotation"]["model"] == "example" and cfg["annotation"]["include_singletons"]
    task = job_queue.get_task(root, tid)
    args = parser().parse_args(next(job_queue.steps(root, task)))
    assert args.reuse_annotations == old and args.include_singletons is None
    with pytest.raises(ValueError, match="true or false"):
        config(override={"annotation":{"include_singletons":"yes"}})
