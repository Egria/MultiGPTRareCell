"""Exercise the model-facing boundary through actual protocol, controller, tools and SQLite."""
import json
from pathlib import Path
import tempfile
import unittest

from rarecell.config import config, fingerprint
from rarecell.discovery import Controller
from rarecell.protocol import parse_message
from test_science import expression_fixture


class ScientificBoundaryTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.root = Path(tmp.name)
        self.addCleanup(tmp.cleanup)
        self.snap = expression_fixture(self.root / "data", 4)
        # Every input has the same opportunity to act as a seed. Choose the first
        # input barcode for this deterministic test, independently of evaluation data.
        items = (self.snap.parent / "discovery_input" / "cells.tsv").read_text().splitlines()[1:]
        self.seed = items[0]
        self.bg = sorted(items[1:25])
        self.catalog = [
            {"id": "Hseed", "members": [self.seed], "markers": [], "sources": [{"kind": "singleton"}]},
            {"id": "Hback", "members": self.bg, "markers": [], "sources": [{"kind": "group"}]},
            {"id": "Hall", "members": sorted(items), "markers": [], "sources": [{"kind": "group"}]},
        ]
        (self.snap / "seeds.jsonl").write_text(json.dumps({"id": "S0", "cell_id": self.seed, "candidate_id": "Hseed"}) + "\n")
        (self.snap / "initial.jsonl").write_text("".join(json.dumps(x) + "\n" for x in self.catalog))
        genes = (self.snap / "genes.tsv").read_text().splitlines()[1:]
        symbols = {"g020": "Gene-A", "g021": "Twin", "g022": "Twin", "g025": "Mouse-only"}
        gene_rows = [f"{g}\t{symbols.get(g,g)}\t{'mouse' if g == 'g025' else 'synthetic'}\n" for g in genes]
        (self.snap / "genes.tsv").write_text("gene_id\tsymbol\tspecies\n" + "".join(gene_rows))
        (self.snap / "manifest.json").write_text(json.dumps({
            "data_version": "fixture-numeric", "snapshot_hash": fingerprint(self.catalog),
            "species": "synthetic", "views": ["V_LINEAR"],
        }))
        (self.snap.parent / "evaluation_fixture" / "canary.txt").write_text("EVALUATION_CANARY_DO_NOT_EXPOSE")
        self.msg = 0

    def controller(self, name, mode="A2", annotations=()):
        ctl = Controller(self.snap, self.root / name, self.root / "shared_cache", config(), mode, annotations)
        ctl.initialize()
        self.addCleanup(ctl.state.close)
        return ctl

    def issue(self, ctl, action):
        self.msg += 1
        raw = json.dumps(action) if isinstance(action, dict) else action
        parsed, = parse_message(raw)
        return ctl.execute(parsed, f"M{self.msg}")

    def test_equivalent_requests_share_numeric_result_and_observation(self):
        ctl = self.controller("queries")
        q = {"action": "coexpress", "candidate_id": "Hseed", "genes": ["Gene-A", "g006"]}
        expected = self.issue(ctl, q)
        equivalent = [
            {"operation": "检查共表达", "arguments": {"gene_list": "g006; g020", "candidate": "Hseed", "threshold": "0.0"}},
            {"action": "coexpress", "members": [self.seed], "logic": "AND", "genes": ["g020", "g006", "g020"]},
            'Answer: ```json\n{"action":"coexpress","candidate_id":"Hseed","genes":["g020","g006"],}\n```',
        ]
        for raw in equivalent:
            answer = self.issue(ctl, raw)
            self.assertEqual(answer["query_key"], expected["query_key"])
            self.assertEqual(answer["evidence_id"], expected["evidence_id"])
        self.assertEqual(ctl.state.usage()["run_reuses"], 3)
        self.assertEqual(ctl.state.usage()["compute_executions"], 1)
        disjunction = self.issue(ctl, dict(q, logic="OR"))
        self.assertNotEqual(disjunction["query_key"], expected["query_key"])
        q = {"action": "compare_groups", "genes": ["g020"], "candidate_id": "Hseed", "background_id": "Hback"}
        by_id = self.issue(ctl, q)
        actual_cells = self.issue(ctl, {"action": "compare_groups", "members": [self.seed],
                                       "background": self.bg, "genes": ["Gene-A"]})
        self.assertEqual(by_id["evidence_id"], actual_cells["evidence_id"])
        new_bg = self.issue(ctl, {"action": "compare_groups", "members": [self.seed], "genes": ["g020"],
                                  "background": self.bg[1:]})
        self.assertNotEqual(by_id["query_key"], new_bg["query_key"])
        self.assertEqual(ctl.state.usage()["new_queries"], 4)

    def test_unlicensed_queries_and_entities_never_touch_cache(self):
        ctl = self.controller("boundary")
        q = {"action": "coexpress", "candidate": "Hseed", "genes": ["g020"]}
        for invalid in [
            dict(q, genes=["g02O"]), dict(q, genes=["Twin"]), dict(q, genes=["Mouse-only"]), dict(q, genes=[" Gene-A "]),
            dict(q, target_gt="rare"), dict(q, weights=[1]), dict(q, direction="down"),
            dict(q, metadata={"scope": "all"}), dict(q, scope="all"), dict(q, data_version="other"),
            {"action": "coexpress", "members": [" " + self.seed], "genes": ["g020"]},
            {"action": "coexpress", "scope": "those cells", "genes": ["g020"]},
            {"action": "get_evidence", "evidence_ids": ["../evaluation_fixture/canary.txt"]},
            {"action": "get_evidence", "candidate_id": "../evaluation_fixture/truth.json"},
            {"action": "readfile", "path": "../evaluation_fixture/truth.json"},
            {"action": "coexpress", "candidate": "Hseed", "genes": ["g020"], "layer": "raw"},
            {"action": "compare_groups", "candidate": "Hseed", "background": [self.seed]},
        ]:
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                self.issue(ctl, invalid)
        self.assertEqual(ctl.state.usage()["new_queries"], 0)
        self.assertEqual(self.issue(ctl, q)["summary"]["n_matching"], 1)
        for origin, complete in (("agent", False), ("tool", True), ("literature", True)):
            parsed, = parse_message(json.dumps(q), origin=origin, complete=complete)
            with self.assertRaises(ValueError):
                ctl.execute(parsed, "M")
        self.assertEqual(ctl.state.usage()["new_queries"], 1)

    def test_ablation_annotation_visibility_and_warm_cache_permissions(self):
        ann = {"candidate_id": "Hseed", "label": "UNIQUE_ANNOTATION_CANARY"}
        q = {"action": "coexpress", "candidate_id": "Hseed", "genes": ["g020"]}
        warm = self.controller("warm", "A3", [ann])
        self.issue(warm, q)
        shuffled = [dict(ann, candidate_id="Hback")]
        original = {c["id"]: c["set_hash"] for c in warm.ordered()}
        for mode in ("A0", "A1", "A2", "A3"):
            ctl = self.controller("mode" + mode, mode, shuffled)
            self.assertEqual({c["id"]: c["set_hash"] for c in ctl.ordered()}, original)
            content = json.dumps([ctl.context(1), self.issue(ctl, {"action": "get_evidence"}),
                                  self.issue(ctl, {"action": "list_candidates"})])
            self.assertNotIn("EVALUATION_CANARY", content)
            self.assertEqual("UNIQUE_ANNOTATION_CANARY" in content, mode in ("A1", "A3"))
            if mode in ("A0", "A1"):
                for query in (q, {"action": "expand_candidate", "candidate": "Hseed"},
                              {"action": "compare_groups", "candidate": "Hseed"}):
                    with self.assertRaises(ValueError):
                        self.issue(ctl, query)
                self.assertIsNone(ctl.science)
                self.assertEqual(ctl.state.usage()["new_queries"], 0)
            else:
                self.assertEqual(ctl.state.usage()["new_queries"], 0)
                self.assertEqual(self.issue(ctl, q)["summary"]["n_matching"], 1)
                self.assertEqual(ctl.state.usage()["new_queries"], 1)
                self.assertEqual(ctl.state.usage()["global_cache_hits"], 1)
                self.assertEqual(ctl.state.usage()["compute_executions"], 0)

    def test_repeated_annotation_input_does_not_enlarge_evidence_or_context(self):
        ann = {"candidate_id": "Hseed", "label": "provisional interpretation"}
        contexts = []
        for n in (1, 3, 10):
            ctl = self.controller("repeat" + str(n), "A1", [ann] * n)
            contexts.append(ctl.context(1))
            cards = ctl.state.all("Evidence")
            self.assertEqual(len(cards), len(self.catalog) + 1)
            derived, = (c for c in cards if c["kind"] == "interpretation")
            self.assertEqual(derived["parent_ids"], ctl.candidate("Hseed")["evidence_ids"])
            self.assertEqual(len(ctl.state.all("Annotation")), 1)
        self.assertEqual(contexts[0], contexts[1])
        self.assertEqual(contexts[0], contexts[2])

    def test_marker_recheck_tracks_derived_members_and_forbids_extension(self):
        ctl = self.controller("lineage")
        expanded = self.issue(ctl, {"action": "retrieve_cells", "candidate": "Hseed"})
        new_id = expanded["candidate_id"]
        self.assertEqual(expanded["added_count"], 3)
        artifact = ctl.state.read_artifact(expanded["artifact_ref"])
        gene = artifact["selection_history"]["signature"][0]["gene_id"]
        # Recheck one discovery-selected marker by both candidate and exact barcode set.
        rechecks = [{"candidate_id": new_id, "threshold": 0},
                    {"members": ctl.candidate(new_id)["members"], "threshold": 0.1}]
        for scope in rechecks:
            out = self.issue(ctl, dict(scope, action="coexpress", genes=[gene]))
            card = ctl.state.get("Evidence", out["evidence_id"])
            self.assertIn(expanded["evidence_id"], card["parent_ids"])
            self.assertTrue(card["selection_history"][0]["selection_dependent"])
            self.assertNotEqual(card["purpose"], "independent_validation")
            self.assertEqual(card["kind"], "observation")
        with self.assertRaises(ValueError):
            self.issue(ctl, {"action": "expand_candidate", "candidate_id": new_id})
        self.assertEqual(ctl.state.usage()["new_queries"], 3)

    def test_rename_reuses_existing_expansion_without_resetting_depth(self):
        ctl = self.controller("rename")
        expansion = {"action": "expand_candidate", "candidate": "Hseed"}
        first = self.issue(ctl, expansion)
        base = ctl.candidate("Hseed")
        rename = self.issue(ctl, {
            "action": "propose_revision", "parents": ["Hseed"], "evidence_ids": base["evidence_ids"],
            "reason": "Retain singleton pending independent membership support",
            "groups": [{"candidate_ids": ["Hseed"], "proposed_label": "pending"}],
        })
        renamed, = rename["candidate_ids"]
        replay = self.issue(ctl, dict(expansion, candidate=renamed))
        self.assertEqual(replay["evidence_id"], first["evidence_id"])
        self.assertEqual(replay["query_key"], first["query_key"])
        self.assertEqual(ctl.state.usage()["compute_executions"], 1)
        with self.assertRaises(ValueError):
            self.issue(ctl, dict(expansion, candidate=renamed, background_id="Hback"))

    def test_absent_implicit_background_is_scientifically_unavailable(self):
        ctl = self.controller("allcells")
        result = self.issue(ctl, {"action": "compare_groups", "candidate_id": "Hall"})
        self.assertEqual(result["status"], "insufficient_background")


if __name__ == "__main__":
    unittest.main()
