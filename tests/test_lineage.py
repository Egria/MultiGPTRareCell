"""Controller-level lineage and provenance tests with real sparse numerical tools."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy import sparse
from sklearn.preprocessing import normalize
from sklearn.random_projection import SparseRandomProjection

from rarecell.config import config
from rarecell.discovery import Controller


class LineageTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        s = self.root / "snapshot"
        v = s / "views"
        v.mkdir(parents=True)
        self.s = s
        self.cells = ["b1", "b2", "b3", "b4", "b5", "b6"]
        (s / "cells.tsv").write_text("cell_id\n" + "\n".join(self.cells) + "\n")
        (s / "genes.tsv").write_text("gene_id\ng1\ng2\ng3\n")
        x = sparse.csr_matrix([[4.,0,1],[4.,0,1],[0,4,1],[0,4,1],[0,4,1],[0,4,1]])
        sparse.save_npz(s / "data.npz", x)
        proj = SparseRandomProjection(n_components=3, density=1.0, random_state=220, dense_output=True)
        np.save(v / "V_LINEAR.npy", proj.fit_transform(normalize(x)))
        (s / "manifest.json").write_text(json.dumps({"views":["V_LINEAR"], "snapshot_hash":"tiny",
                                                      "data_version":"lineage-test"}))
        rows = [{"id": cid, "members": members, "sources": [], "markers": []}
                for cid, members in [("a",["b1"]),("b",["b3"]),("all",self.cells)]]
        (s / "initial.jsonl").write_text("\n".join(map(json.dumps, rows)) + "\n")
        (s / "seeds.jsonl").write_text("")
        self.cfg = config()

    def controller(self, mode="A2", annotations=(), run="run"):
        c = Controller(self.s, self.root / run, self.root / "cache", self.cfg, mode, annotations)
        self.addCleanup(c.state.close)
        c.initialize()
        return c

    def test_expansion_depth_after_identical_revision_and_split(self):
        c = self.controller()
        base = c.candidate("a")
        first = c.scientific("expand_candidate", {"candidate_id": "a"})
        self.assertEqual(first["added_count"], 1)
        with c.state.transaction():
            alias = c.make_revision(base["members"], ["b"], "unnamed", "provisional",
                                    base["evidence_ids"], "revision", "Reconsider membership")
        duplicate = c.scientific("expand_candidate", {"candidate_id": alias["id"]})
        self.assertEqual(duplicate["evidence_id"], first["evidence_id"])
        self.assertEqual(c.state.usage()["new_queries"], 1)
        self.assertEqual(c.state.usage()["compute_executions"], 1)
        self.assertEqual(len(c.state.all("Investigation")), 2)
        self.assertEqual(c.state.get_meta("depths"), {"a": 1, "b": 1})
        for cid in (alias["id"], "b", first["candidate_id"]):
            with self.assertRaisesRegex(ValueError, "depth"):
                c.scientific("expand_candidate", {"candidate_id": cid, "background":["b4","b5","b6"]})

    def test_result_transaction_failure_does_not_leak_in_memory_revision(self):
        c = self.controller()
        real_put = c.state.put
        def abort_after_candidate(kind, identifier, obj):
            if kind == "Investigation":
                raise OSError("injected DB registration failure")
            return real_put(kind, identifier, obj)
        with patch.object(c.state, "put", side_effect=abort_after_candidate):
            with self.assertRaises(OSError):
                c.scientific("expand_candidate", {"candidate_id":"a"})
        self.assertEqual(set(c.candidates), {"a","b","all"})
        self.assertEqual(c.state.all("QueryObservation"), [])
        self.assertEqual(c.state.get_meta("depths",{}), {})
        result = c.scientific("expand_candidate", {"candidate_id":"a"})
        self.assertIn(result["candidate_id"], c.candidates)
        self.assertEqual(c.state.usage()["compute_executions"], 1)
        self.assertEqual(c.state.usage()["run_reuses"], 1)

    def test_repeated_numeric_query_tracks_selection_without_new_observation(self):
        c = self.controller()
        hit = c.scientific("coexpress", {"candidate_id":"a","genes":["g1"]})
        edit = {"parents":["a"],"groups":[{"members":["b1"],"proposed_label":"tentative"}],
                "evidence_ids":[hit["evidence_id"]],"reason":"re-evaluate this singleton"}
        cmd = {"action":"propose_revision", "args":edit}
        one = c.execute(cmd, "m1")
        two = c.execute(cmd, "m2")
        self.assertEqual(one["candidate_ids"], two["candidate_ids"])
        selected = c.candidate(one["candidate_ids"][0])
        interps = [eid for eid in selected["evidence_ids"] if c.state.get("Evidence",eid)["kind"]=="interpretation"]
        self.assertEqual(len(interps), 1)
        repeat = c.scientific("coexpress", {"members":["b1"],"genes":["g1"]})
        self.assertEqual(repeat["evidence_id"], hit["evidence_id"])
        self.assertNotEqual(repeat["query_use_id"], hit["query_use_id"])
        self.assertIn(interps[0], c.state.get("QueryUse", repeat["query_use_id"])["selection_evidence_ids"])
        self.assertEqual(len(c.state.all("QueryObservation")), 1)
        c.execute({"action":"submit","args":{"status":"provisional"}}, "m3")
        self.assertEqual(c.state.get_meta("submission_status"), "provisional")
        self.assertEqual(c.state.all("Decision")[-1]["status"], "provisional")

    def test_scientific_evidence_pages_actual_hits_and_scores(self):
        c = self.controller()
        hit = c.scientific("coexpress", {"scope":"all","genes":["g1"]})
        tool = "get_evidence"
        args = {"evidence_ids":[hit["evidence_id"]],"limit":1,"offset":1}
        page = c.read_tool(tool, args)
        self.assertEqual(page["sections"]["matched_members"],{"total":2,"rows":["b2"]})
        self.assertEqual(page["sections"]["gene_counts"]["total"],1)
        expand = c.scientific("expand_candidate", {"candidate_id":"a"})
        page = c.read_tool(tool, {"evidence_ids":[expand["evidence_id"]]})
        self.assertEqual(page["sections"]["scores"]["rows"][0]["cell_id"], "b2")
        self.assertGreater(page["sections"]["scores"]["rows"][0]["score"], .72)
        compare = c.scientific("compare_groups", {"candidate_id":"a","genes":["g1"]})
        page = c.read_tool(tool, {"evidence_ids":[compare["evidence_id"]]})
        self.assertEqual(page["sections"]["statistics"]["rows"][0]["mean_query"], 4.0)
        result = c.scientific("compare_groups", {"candidate_id":"all"})
        self.assertEqual(result["status"], "insufficient_background")

    def test_annotation_consistency_and_no_label_mode(self):
        # Test controller consumption of possible annotation outputs only.
        first = {"candidate_id":"a","label":"tentative-1","partition_id":"p1","resolution":0.5}
        second = dict(first, label="tentative-2", partition_id="p2", resolution=1.0)
        annotations = [first, first, second]
        labeled = self.controller("A1", annotations)
        label = labeled.current_label("a")
        self.assertEqual(label["label"], "unknown")
        self.assertEqual(label["labels"], ["tentative-1","tentative-2"])
        self.assertEqual(len(label["evidence_ids"]),2)
        self.assertEqual(len(labeled.state.all("Annotation")),2)
        for eid in label["evidence_ids"]:
            self.assertEqual(labeled.state.get("Evidence",eid)["parent_ids"],labeled.candidate("a")["evidence_ids"])
        fresh = self.controller("A1")
        self.assertEqual(fresh.current_label("a"),label)
        unlabelled = self.controller("A0", annotations, run="passive")
        self.assertEqual(unlabelled.state.all("Annotation"),[])
        self.assertNotIn("tentative-1", json.dumps(unlabelled.context(1)))
        self.assertEqual(set(labeled.candidates),set(unlabelled.candidates))
