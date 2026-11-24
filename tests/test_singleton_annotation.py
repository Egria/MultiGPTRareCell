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
    old_cfg = {"model": "example", "dual_backgrounds": False, "rscript": "/unavailable"}
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
    cfg = dict(model="example", rscript="/unavailable", include_singletons=True, dual_backgrounds=False)
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


def test_dual_panels_reach_discovery_and_recover_partial(tmp_path):
    from rarecell import contrast_panels, exhaustive
    from rarecell.discovery import Controller

    snap, out = tmp_path / "snap", tmp_path / "ann"
    manifest = snapshot(snap)
    path = snap / "initial.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    rows[-1]["sources"].append(dict(kind="partition", view="V_HVG", partition_id="P2", resolution=1))
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    manifest["prepare"] = {"k": 1}
    manifest["stages"]["seurat_markers"] = {"status": "blocked"}
    write_json(snap / "manifest.json", manifest)
    views = snap / "views"
    views.mkdir()
    (views / "V_HVG_graph_rows.tsv").write_text("cell_id\na\nb\nc\ne\n")
    sparse.save_npz(views / "V_HVG_snn.npz", sparse.csr_matrix(np.array([
        [0, 1, 0, 1], [1, 0, 1, 1], [0, 1, 0, 0], [1, 1, 0, 0]])))
    np.save(views / "V_HVG.npy", np.array([[2, 0], [2, 1], [1, 1], [0, 0], [0, 1.]]))
    options = dict(model="fixture", rscript="/absent", batch_size=3, panel_mode="dual")
    settings = annotation._resolve(options)
    assert settings["include_singletons"] and settings["dual_backgrounds"]
    jobs, skipped = annotation._jobs(snap, manifest, settings)
    pairs = {(p["candidate_id"], p["background_type"]): p for p in
             skipped + [p for job in jobs for p in job["groups"]]}
    assert len(pairs) == 8  # three seeds and one repeated group
    seed_local, seed_global = [pairs["Ha", scope] for scope in ("neighbor", "global")]
    assert [m["symbol"] for m in seed_local["markers"]] == ["High", "g4"]
    assert seed_global["n_background"] == 4 and seed_local["n_background"] == 2
    assert [m["symbol"] for m in seed_global["markers"]] == ["High", "g4", "Equal"]
    assert pairs["Hp", "neighbor"]["n_background"] == 1 and len(pairs["Hp", "global"]["partitions"]) == 2
    assert pairs["Hp", "neighbor"]["background_hash"] == fingerprint(["e"])
    assert pairs["Hd", "neighbor"]["reason"] == pairs["Hd", "global"]["reason"] == "no_expression"
    assert all(len({p["background_type"] for p in job["groups"]}) == 1 for job in jobs)

    def reply(job):
        batch = out / "batches" / fingerprint(job)
        batch.mkdir(parents=True, exist_ok=True)
        labels = [dict(candidate_id=p["candidate_id"], label=p["background_type"] + " label", ontology_id=None)
                  for p in job["groups"]]
        write_json(batch / "response.json", dict(valid=True, commit=annotation.PIN, raw_response="fixture reply", rows=labels))

    for job in jobs[1:]:
        reply(job)
    assert annotation._run(snap, out, options, manifest)["status"] == "partial"
    partial = annotation.verify_annotation(snap, out, _manifest=manifest, allow_partial=True)
    assert len(partial) == 8 and any(p["status"] == "unavailable" for p in partial)
    reply(jobs[0])
    assert annotation._run(snap, out, options, manifest)["status"] == "success"
    annotations = annotation.verify_annotation(snap, out, _manifest=manifest)
    assert len(annotations) == 8
    cfg = config(override={"exhaustive_seed_review": True})
    ctrl = Controller(snap, tmp_path / "run", tmp_path / "cache", cfg, "A3", annotations)
    try:
        ctrl.initialize()
        card = ctrl.card(ctrl.candidate("Ha"))
        assert set(card["tentative_alternatives"]) == {"neighbor label", "global label"}
        assert {r["background_type"]: set(r["genes"]) for r in card["annotation_sources"]} == {
            "neighbor": set(seed_local["genes"]), "global": set(seed_global["genes"])}
        for obj in card["annotation_sources"]:
            entry = ctrl.state.get("Evidence", obj["evidence_id"])
            parent = ctrl.state.get("Evidence", entry["parent_ids"][0])
            assert parent["kind"] == "observation" and parent["method_version"] == contrast_panels.SOURCE
        ctrl.state.put("SeedMeasurement", "Ha", {"candidate_id": "Ha", "status": "no_expansion", "evidence_id": "obs"})
        exhaustive.prepare_batch(ctrl, "T1", [{"candidate_id": "Ha", "previous_batches": 0}])
        batch = ctrl.state.get("SeedBatch", "T1")
        shown = ctrl.state.read_artifact(batch["context_ref"])["seeds"][0]
        assert shown["annotation_sources"] == card["annotation_sources"]
        label_eid = shown["annotation_sources"][0]["evidence_id"]
        result = exhaustive.record(ctrl, "T1", dict(candidate_id="Ha", evidence_ids=[label_eid],
                                   disposition="follow_up", label_reliability="uncertain", reason="Resolve the conflicting labels"))
        assert result["evidence_id"] == "obs" and result["annotation_evidence_ids"] == [label_eid]
    finally:
        ctrl.state.close()
    sparse.save_npz(views / "V_HVG_snn.npz", sparse.csr_matrix((4, 4)))
    fallback = next(p for p in contrast_panels.panels(snap, manifest, settings)
                    if p["candidate_id"] == "Hp" and p["background_type"] == "neighbor")
    assert fallback["background_rule"] == "HVG_nearest_external_to_centroid" and fallback["n_background"] == 1


def test_global_only_reuses_dual_global_receipts_and_builds_actionable_seed_pool(tmp_path):
    from rarecell.discovery import Controller
    from rarecell.protocol import parse_message

    snap, old, new = (tmp_path / x for x in ("snap", "dual", "global"))
    manifest = snapshot(snap)
    manifest["prepare"] = {"k": 1}
    (snap / "views").mkdir()
    np.save(snap / "views/V_HVG.npy", np.ones((5, 2)))
    write_json(snap / "manifest.json", manifest)
    with (snap / "initial.jsonl").open("a") as out:
        out.write(json.dumps(dict(id="Hb", members=["b"], sources=[dict(kind="outlier", view="V_GINI")],
                                  markers=[])) + "\n")
    with (snap / "seeds.jsonl").open("a") as out:
        out.write(json.dumps(dict(id="b", cell_id="b", candidate_id="Hb")) + "\n")
    dual = dict(model="fixture", rscript="/unavailable", dual_backgrounds=True)
    opts = dict(model="fixture", rscript="/unavailable")
    old_jobs, _ = annotation._jobs(snap, manifest, annotation._resolve(dual))
    jobs, skipped = annotation._jobs(snap, manifest, annotation._resolve(opts))
    assert jobs == [job for job in old_jobs if job["groups"][0]["background_type"] == "global"]
    assert all(r["background_type"] == "global" for r in skipped)
    for job in old_jobs:
        group, = job["groups"]
        global_panel = group["background_type"] == "global"
        cid = group["candidate_id"]
        batch = old / "batches" / fingerprint(job)
        batch.mkdir(parents=True)
        label = "Schwann cell" if cid == "Ha" else "Schwann Cell"
        row = dict(candidate_id=cid, label=label, ontology_id="CL:0002573" if global_panel else "CL:0000003")
        write_json(batch / "response.json", dict(commit=annotation.PIN, valid=True, raw_response=label, rows=[row]))
    assert annotation._run(snap, old, dual, manifest)["status"] == "success"
    status = annotation._run(snap, new, opts, manifest, reuse=old)
    assert status["cached_batches"] == len(jobs) and status["newly_completed_batches"] == 0
    annotations = annotation.verify_annotation(snap, new, _manifest=manifest)
    expected = {g: sum(g in r["selection_history"]["genes"] for r in annotations if r["candidate_id"] in {"Ha", "Hb"})
                for r in annotations if r["candidate_id"] in {"Ha", "Hb"} for g in r["selection_history"]["genes"]}
    cfg = config()
    run, cache = tmp_path / "discovery", tmp_path / "cache"
    for i in range(2):
        ctrl = Controller(snap, run, cache, cfg, "A3", annotations)
        try:
            ctrl.initialize()
            ctrl.state.set_meta("policy", "live")
            pool, = ctrl.context(1)["ontology_seed_pools"]
            cid = pool["candidate_id"]
            assert pool["n_seeds"] == 2 and pool["ontology_id"] == "CL:0002573"
            cmd, = parse_message('{"action":"list seed pools"}')
            assert ctrl.execute(cmd, "test")["ontology_seed_pools"] == [pool]
            page = ctrl.read_tool("get_evidence", dict(candidate_id=cid))
            assert page["seed_candidate_ids"] == ["Ha", "Hb"] and page["members"] == ["a", "b"]
            assert page["marker_seed_counts"] == expected
            assert ctrl.state.get_meta("selected", []) == []
            query = ctrl.scientific("coexpress", dict(candidate_id=cid, genes=["g1"]))
            assert query["summary"]["n_query"] == 2
            if i == 0:
                first = cid
            else:
                assert cid == first
                assert ctrl.decide("submit", {"candidate_ids":[cid]}, "test")["candidate_ids"] == [cid]
        finally:
            ctrl.state.close()


def test_ontology_pool_counts_distinct_global_seed_panels_only():
    from rarecell.contrast_panels import ontology_seed_pools

    def row(cid, term, genes, scope="global", status="annotated", label="example"):
        return dict(candidate_id=cid, ontology_id=term, label=label, background_type=scope,
                    status=status, evidence_id=cid + scope, selection_history={"genes": genes})

    term = "CL:0002573"
    rows = [row("A", term, ["G", "G", "X"]), row("B", term, ["G", "Y"]),
            row("A", term, ["Z"], scope="neighbor"), row("C", term, ["Z"], scope="neighbor"),
            row("D", None, ["Z"]), row("E", term, ["Z"], status="unavailable"),
            row("F", term, ["Z"], label="unknown"), row("group", term, ["Z"]),
            row("H", "CL:0005019", ["Y"])]
    result, = ontology_seed_pools(set("ABCDEFH"), rows + rows[:1])
    assert result[0] == term and result[2] == ["A", "B"]
    assert result[4] == {"G": 2, "X": 1, "Y": 1}
