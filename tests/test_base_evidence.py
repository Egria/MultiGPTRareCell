"""Frozen scientific tables in passive/active contexts, with real R numerical output."""
import csv
import os
import shutil
from pathlib import Path

import pytest

from rarecell.base_evidence import DE, HVG
from rarecell.config import ADAPTIVE_MODE, config, fingerprint
from rarecell.discovery import Controller
from rarecell.prepare import prepare


@pytest.mark.skipif(not shutil.which(os.environ.get("RARECELL_RSCRIPT", "Rscript")),
                    reason="requires pinned R/Seurat for real precomputed DE and HVG tables")
def test_shared_seurat_tables_and_annotation_dependencies(tmp_path):
    path = Path(__file__).resolve().parents[1]/"examples"/"minimal_sparse"/"input"/"dataset.yaml"
    cfg = config(override={"prepare": {"views": ["V_HVG"], "n_features": 25, "n_pcs": 8, "k": 7,
                                      "rscript": os.environ.get("RARECELL_RSCRIPT", "Rscript")}})
    snapshot = prepare(path, tmp_path/"snapshot", cfg)
    with (snapshot/"seurat"/"markers.tsv").open() as f:
        markers = list(csv.DictReader(f, delimiter="\t"))
    with (snapshot/"seurat"/"hvg.tsv").open() as f:
        hvgs = list(csv.DictReader(f, delimiter="\t"))
    results = []
    for mode in ("A0","A1","A2","A3", ADAPTIVE_MODE):
        ctrl = Controller(snapshot, tmp_path/mode, tmp_path/"cache", cfg, mode)
        try:
            ctrl.initialize()
            hvg_table = ctrl.context(1)["global_feature_selection"][0]
            def page(eid, offset, limit=2):
                return ctrl.execute({"action": "get_evidence", "args": {
                    "evidence_ids": [eid], "offset": offset, "limit": limit}}, "read")
            hp = page(hvg_table["evidence_id"], 2)
            assert hp["sections"][HVG]["rows"] == hvgs[2:4]
            cid = next(r["group_id"] for r in markers if sum(k["group_id"]==r["group_id"] and
                       k["resolution"]==r["resolution"] for k in markers) >= 4)
            candidate = ctrl.candidate(cid)
            table = ctrl.card(candidate)["seurat_wilcoxon"][0]
            stats = candidate["structure"]["V_HVG"]
            assert ctrl.card(candidate)["structure"]["V_HVG"]["within_fraction"] == pytest.approx(stats["within_fraction"])
            assert 0 <= stats["within_fraction"] <= 1 and stats["graph_coverage"] == 1
            full = [r for r in markers if r["group_id"]==cid and float(r["resolution"])==table["resolution"]]
            info = next(t for t in candidate["seurat_markers"] if t["resolution"]==table["resolution"])
            assert info["marker_hash"] == fingerprint(full)
            dp = page(table["evidence_id"], 1, 3)
            assert dp["sections"][DE]["total"] == len(full)
            assert dp["purpose"] == "selection_dependent" and len(dp["parent_ids"])==1
            observation = ctrl.state.get("Evidence", dp["parent_ids"][0])
            assert hvg_table["evidence_id"] in observation["parent_ids"]
            read = dp["sections"][DE]["rows"]
            assert [r["gene_id"] for r in read] == [r["gene_id"] for r in full[1:4]]
            assert read[0]["avg_log2FC"] == pytest.approx(float(full[1]["avg_log2FC"]))
            assert dp["summary"]["n_foreground"] + dp["summary"]["n_background"] == 160
            assert all(not c.get("seurat_markers") for c in ctrl.ordered() if len(c["members"]) < 3)
            assert ctrl.state.usage()["new_queries"] == ctrl.state.usage()["compute_executions"] == 0
            assert ctrl.science is None
            results.append((hvg_table["evidence_id"], cid, table["evidence_id"], info, full))
        finally:
            ctrl.state.close()
    assert all(x==results[0] for x in results[1:])
    _, cid, de_id, info, _ = results[0]
    # Test the provenance link with an explicitly supplied interpretation; no
    # label service is called or inferred during this controller-level check.
    ann = {"candidate_id": cid, "label":"test interpretation", "source":"OntoAnno::gptcelltype",
           **{k:info[k] for k in ("marker_hash","partition_id","resolution")}}
    for valid in [True, False]:
        row = ann if valid else dict(ann, marker_hash="different foreground markers")
        ctrl = Controller(snapshot, tmp_path/("valid" if valid else "bad"), tmp_path/"cache", cfg, "A1",[row])
        try:
            if valid:
                ctrl.initialize()
                label = ctrl.state.all("Annotation")[0]["evidence_id"]
                assert ctrl.state.get("Evidence", label)["parent_ids"]==[de_id]
            else:
                with pytest.raises(ValueError, match="frozen Seurat marker"):
                    ctrl.initialize()
        finally:
            ctrl.state.close()
