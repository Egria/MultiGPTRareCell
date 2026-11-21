"""Durability/integrity tests, including real process termination and concurrent runs."""

import json
import multiprocessing as mp
import os
from pathlib import Path
import tempfile
import unittest

from rarecell.state import BudgetExhausted, IntegrityError, QueryBlocked, State, digest, key_for


def request(**overrides):
    value = {"action": "compare_groups", "data_version": "fixture-v1", "method_version": "v1",
             "members": ["c1", "c2"], "background": ["c3"], "input_layer": "log1p", "seed": 13}
    value.update(overrides)
    return value


def card(**overrides):
    value = {"object": "candidate-c1-c2", "kind": "observation", "parent_ids": [], "source_family": "RNA",
             "cell_set_hash": digest(["c1","c2"]), "data_version": "fixture-v1", "method_version": "v1",
             "config_version": "config-1", "prompt_version": "none", "selection_history": [],
             "artifact_refs": [], "purpose": "exploration", "summary": {"mean": 2.5}}
    value.update(overrides)
    return value


def compute_counted(directory):
    with (Path(directory) / "calls").open("a") as handle:
        handle.write("compute\n")
    return {"mean": sum([1.0, 4.0]) / 2}


def crash_worker(directory):
    class CrashBeforeStateCommit(State):
        def _finish(self, *args, **kwargs):
            # _install_json has returned. OS death releases the flock and closes SQLite.
            os._exit(83)
    root = Path(directory)
    state = CrashBeforeStateCommit(root / "run", root / "cache")
    state.query(request(), 1, lambda: compute_counted(root))


def concurrent_worker(directory, index, start):
    root = Path(directory)
    state = State(root / ("worker" + str(index)), root / "cache")
    start.wait(15)
    state.query(request(), 1, lambda: compute_counted(root))
    state.close()


class StateTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.state = State(self.root / "run", self.root / "cache")
        self.addCleanup(self.state.close)

    def test_immutable_records_and_atomic_revisions(self):
        state = self.state
        state.put("CellSet", "s1", {"id": "s1", "members": ["c1","c2"]})
        with self.assertRaises(IntegrityError):
            with state.transaction():
                state.put("Revision", "split-left", {"members": ["c1"]})
                state.set_meta("latest", "split-left")
                state.event("revision_created", {"id": "split-left"})
                state.put("CellSet", "s1", {"id": "s1", "members": ["c2"]})
        self.assertEqual(state.all("Revision"), [])
        self.assertEqual(state.get_meta("latest"), None)
        self.assertEqual((state.run / "events.jsonl").read_text(), "")
        self.assertEqual(state.get("CellSet", "s1")["members"], ["c1","c2"])

    def test_evidence_deduplication_and_dag(self):
        s = self.state
        e = s.evidence(card(evidence_id="base"))
        for repeats in (1, 3, 10):
            for _ in range(repeats):
                self.assertEqual(s.evidence(card()), e)
            self.assertEqual(len(s.all("Evidence")), 1)
        label = s.evidence(card(kind="interpretation", parent_ids=[e], summary="OntoAnno tentative name"))
        explanation = s.evidence(card(kind="hypothesis", parent_ids=[label,e,e], summary="possible subgroup"))
        self.assertEqual(s.get("Evidence", explanation)["parent_ids"], sorted([label,e]))
        other = s.evidence(card(method_version="v2", summary={"median": 2.0}))
        self.assertNotEqual(other, e)
        with self.assertRaises(IntegrityError):
            s.evidence(card(evidence_id="base", parent_ids=[label]))  # would close a two-node cycle
        with self.assertRaises(IntegrityError):
            s.evidence(card(evidence_id="self", parent_ids=["self"]))
        with self.assertRaises(IntegrityError):
            s.evidence(card(parent_ids=["missing"]))
        with self.assertRaises(IntegrityError):
            s.evidence(card(evidence_id="base", summary="silently changed"))
        with self.assertRaises(IntegrityError):
            s.evidence(card(evidence_id="new-vote"))
        self.assertEqual(len(s.all("Evidence")), 4)

    def test_selection_dependent_claim_cannot_promote(self):
        s = self.state
        selected = s.evidence(card(selection_history=[{"selected_genes": ["G1"]}]))
        checked = card(parent_ids=[selected], selection_history=[{"selected_genes": ["G1"]}],
                       purpose="selection_dependent", summary={"test_genes": ["G1"], "mean": 2.5})
        checked_id = s.evidence(checked)
        self.assertEqual(s.get("Evidence", checked_id)["purpose"], "selection_dependent")
        with self.assertRaises(IntegrityError):
            s.evidence(dict(checked, purpose="independent_validation"))
        with self.assertRaises(IntegrityError):
            s.evidence(dict(checked, kind="validation"))
        for purpose in (None, "Independent Validation", "arbitrary"):
            with self.assertRaises(IntegrityError):
                s.evidence(dict(checked, purpose=purpose))
        for ref in ("../evaluation/gt.tsv", "artifacts/" + "0" * 64 + ".json"):
            with self.assertRaises(IntegrityError):
                s.evidence(dict(checked, artifact_refs=[ref]))
        with self.assertRaises(ValueError):
            s.put("Evidence", "bypass", dict(checked, purpose="independent_validation"))

    def test_scope_method_and_replication_identity(self):
        q = request()
        self.assertEqual(key_for(q), key_for(dict(reversed(list(q.items())))))
        for key, value in (("background", ["c2"]), ("members", ["c1"]), ("logic", "OR"),
                           ("data_version", "new"), ("method_version", "v2"), ("seed", 37)):
            self.assertNotEqual(key_for(q), key_for(dict(q, **{key: value})))

    def test_run_budget_applies_before_cache(self):
        state = self.state
        calc = lambda: compute_counted(self.root)
        result,key = state.query(request(), 1, calc)
        self.assertEqual(state.read_artifact(state.result_ref(key)), result)
        again,_ = state.query(request(), 0, calc)
        self.assertEqual(again, result)
        b = State(self.root / "b", self.root / "cache")
        self.addCleanup(b.close)
        with self.assertRaises(BudgetExhausted):
            b.query(request(), 0, calc)
        self.assertEqual(b.all("CanonicalQuery"), [])
        b.query(request(), 1, calc)
        with self.assertRaises(BudgetExhausted):
            b.query(request(seed=42), 1, calc)
        self.assertEqual(b.usage()["global_cache_hits"], 1)
        self.assertEqual(b.usage()["new_queries"], 1)
        self.assertEqual(b.usage()["compute_executions"], 0)
        self.assertEqual(state.usage()["run_reuses"], 1)
        self.assertEqual((self.root / "calls").read_text().splitlines(), ["compute"])

    def test_mutated_artifact_and_cache_raise(self):
        s = self.state
        _,key = s.query(request(), 1, lambda: compute_counted(self.root))
        with (s.run / s.result_ref(key)).open("ab") as handle:
            handle.write(b" ")
        with self.assertRaises(IntegrityError):
            s.query(request(), 1, lambda: compute_counted(self.root))
        other = State(self.root / "other", s.cache)
        self.addCleanup(other.close)
        path = s._entry_path(key)
        changed = json.loads(path.read_bytes())
        changed["result"]["mean"] = -9.0
        path.write_text(json.dumps(changed))
        with self.assertRaises(IntegrityError):
            other.query(request(), 1, lambda: compute_counted(self.root))

    def test_process_crash_between_artifact_write_and_db_commit(self):
        p = mp.get_context("spawn").Process(target=crash_worker, args=(str(self.root),))
        p.start()
        p.join(30)
        self.assertEqual(p.exitcode, 83)
        restored = State(self.state.run, self.state.cache)
        self.addCleanup(restored.close)
        self.assertEqual(restored.usage()["recovered_results"], 1)
        self.assertEqual(restored.all("ExecutionAttempt")[0]["status"], "success")
        result, _ = restored.query(request(), 1, lambda: compute_counted(self.root))
        self.assertEqual(result, {"mean": 2.5})
        self.assertEqual((self.root / "calls").read_text().splitlines(), ["compute"])
        self.assertEqual(restored.usage()["compute_executions"], 1)

    def test_parallel_runs_single_physical_execution(self):
        ctx = mp.get_context("spawn")
        start = ctx.Barrier(2)
        workers = [ctx.Process(target=concurrent_worker, args=(str(self.root), i, start)) for i in range(2)]
        for p in workers:
            p.start()
        for p in workers:
            p.join(30)
            self.assertEqual(p.exitcode, 0)
        usage = []
        for i in range(2):
            s = State(self.root / ("worker" + str(i)), self.state.cache)
            usage.append(s.usage())
            s.close()
        self.assertEqual(sum(x["global_cache_hits"] for x in usage), 1)
        self.assertEqual(sum(x["new_queries"] for x in usage), 2)
        self.assertEqual((self.root / "calls").read_text().splitlines(), ["compute"])

    def test_transient_retry_and_local_permanent_failure(self):
        state = self.state
        def filesystem_compute():
            path = self.root / "values"
            if not path.exists():
                path.write_text("2\n3\n")
                raise TimeoutError("staged first read")
            return {"sum": sum(map(float,path.read_text().splitlines()))}

        result,key = state.query(request(), 3, filesystem_compute)
        self.assertEqual(result["sum"], 5.0)
        self.assertEqual(state.usage()["technical_retries"], 1)
        self.assertEqual([a["status"] for a in state.all("ExecutionAttempt")], ["retryable","success"])
        def invalid_input():
            return {"mean": float("nan")}
        with self.assertRaises(QueryBlocked):
            state.query(request(seed=1), 3, invalid_input)
        with self.assertRaises(QueryBlocked):
            state.query(request(seed=1), 3, invalid_input)
        self.assertEqual(state.usage()["compute_executions"], 3)
        self.assertEqual(state.get("CanonicalQuery", key)["status"], "success")
        state.query(request(seed=2), 3, filesystem_compute)

    def test_retry_exhaustion_is_bounded_and_recoverable(self):
        def missing_file():
            return {"sum": sum(map(float, (self.root / "input").read_text().splitlines()))}
        s = self.state
        for _ in range(3):
            with self.assertRaises(QueryBlocked):
                s.query(request(), 1, missing_file)
        self.assertEqual(s.usage()["compute_executions"], 2)
        (self.root / "input").write_text("5\n3\n")
        result,_ = s.query(request(), 1, missing_file, max_attempts=3)
        self.assertEqual(result["sum"], 8)
        self.assertEqual(s.usage()["new_queries"], 1)

    def test_repair_torn_event_export(self):
        self.state.event("first", {"e": 1})
        self.state.event("second", {"e": 2})
        path = self.state.run / "events.jsonl"
        original = path.read_bytes()
        path.write_bytes(original[:-8])
        s = State(self.state.run, self.state.cache)
        self.addCleanup(s.close)
        self.assertEqual(path.read_bytes(), original)
        s.event("third", {"e": 3})
        self.assertEqual([json.loads(line)["seq"] for line in path.read_text().splitlines()], [1,2,3])
