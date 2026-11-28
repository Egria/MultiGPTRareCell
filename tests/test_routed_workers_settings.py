from pathlib import Path

import pytest

from rarecell import job_queue
from rarecell.cli import parser
from rarecell.config import ROUTED_MODE, config, read_json, write_json
from rarecell.routed import RoutedDiscovery
from rarecell.run import discover, resume
from test_routed import setup_input


def test_resume_changes_parallelism_without_reopening_budgets(tmp_path, monkeypatch):
    cfg, snapshot, annotation, ann = setup_input(tmp_path)
    run = tmp_path / "run"
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *args, **kwargs: ann)
    invocations = []

    def routed_loop(self, obtain, pause_after, retry_unknown):
        invocations.append((self.ctrl.cfg["routed_workers"], self.ctrl.budget["max_new_queries"]))
        if len(invocations) == 1:
            return "partial", "routed_query_budget"
        if len(invocations) < 4:
            return "partial", "paused"
        self.ctrl.decide("submit", {"candidate_ids": [], "reason": "Finish the test"}, None, source="agent")
        return "complete", "submitted"

    monkeypatch.setattr(RoutedDiscovery, "run", routed_loop)
    stopped = discover(snapshot, run, tmp_path / "cache", cfg, annotations=annotation)
    original_config = (run / "config_resolved.yaml").read_bytes()
    original_hash = read_json(run / "checkpoint.json")["cfg_hash"]
    assert resume(run, routed_workers=4) == stopped
    assert "routed_workers" not in read_json(run / "checkpoint.json")
    resume(run, max_new_queries=900, routed_workers=4)
    resume(run)
    final = resume(run, routed_workers=1)
    assert final["submitted"]
    assert invocations == [(1, 800), (4, 900), (4, 900), (1, 900)]
    assert resume(run, routed_workers=8) == final
    check = read_json(run / "checkpoint.json")
    assert (check["cfg_hash"], check["routed_workers"]) == (original_hash, 1)
    assert len(check["budget_extensions"]) == 1
    assert (run / "config_resolved.yaml").read_bytes() == original_config


@pytest.mark.parametrize("value", [0, 17, True, "4", 2.5])
def test_invalid_workers(value):
    with pytest.raises(ValueError, match="routed_workers"):
        config(override={"mode": ROUTED_MODE, "routed_workers": value})


def test_worker_options_in_cli_and_queue(tmp_path):
    root = tmp_path / "tasks"
    job_queue.initialize(root)
    dataset = Path(__file__).resolve().parents[1] / "examples/minimal_sparse/input/dataset.yaml"
    task = job_queue.submit(root, tmp_path, dict(kind="pipeline", mode=ROUTED_MODE, dataset=str(dataset),
                                                 routed_workers="4", model="test", annotation_model="test"))
    steps = list(job_queue.steps(root, job_queue.get_task(root, task)))
    args = parser().parse_args(steps[-1])
    assert config(args.config)["routed_workers"] == 4
    assert parser().parse_args([*steps[-1], "--routed-workers", "2"]).routed_workers == 2

    run = tmp_path / "old"
    run.mkdir()
    write_json(run / "checkpoint.json", {"mode": ROUTED_MODE})
    task = job_queue.submit(root, tmp_path, {"kind": "resume", "run": str(run), "routed_workers": "1"})
    args = parser().parse_args(next(job_queue.steps(root, job_queue.get_task(root, task))))
    assert args.cmd == "resume" and args.routed_workers == 1 and args.max_new_queries is None
    with pytest.raises(ValueError, match="only to A3_routed"):
        config(override={"mode": "A2", "routed_workers": 4})


@pytest.mark.parametrize("policy", ["single", "multiple"])
def test_membership_policy_from_dashboard_and_replay_cli(tmp_path, policy):
    root = tmp_path / "tasks"
    job_queue.initialize(root)
    dataset = Path(__file__).resolve().parents[1] / "examples/minimal_sparse/input/dataset.yaml"
    form = dict(kind="pipeline", mode=ROUTED_MODE, dataset=str(dataset),
                routed_conflict_policy=policy, model="test", annotation_model="test")
    task = job_queue.submit(root, tmp_path, form)
    steps = list(job_queue.steps(root, job_queue.get_task(root, task)))
    assert config(parser().parse_args(steps[-1]).config)["routed_conflict_policy"] == policy
    replay = ["reconcile-routed", "--run", "old", "--out", "new", "--snapshot", "snapshot"]
    assert parser().parse_args(replay + ["--routed-conflict-policy", policy]).routed_conflict_policy == policy
    with pytest.raises(ValueError):
        job_queue.submit(root, tmp_path, dict(form, mode="A3"))


def test_dashboard_reruns_frozen_votes_in_new_directory(tmp_path):
    root = tmp_path / "tasks"
    job_queue.initialize(root)
    run, snapshot = tmp_path / "original", tmp_path / "snapshot"
    run.mkdir()
    snapshot.mkdir()
    write_json(run / "run_summary.json", dict(mode=ROUTED_MODE, submitted=True))
    write_json(run / "forced_review_book.json", {})
    gt = tmp_path / "gt.tsv"
    gt.touch()
    form = dict(kind="reconcile", run=str(run), snapshot=str(snapshot), gt=str(gt), routed_conflict_policy="multiple")
    tid = job_queue.submit(root, tmp_path, form)
    task = job_queue.get_task(root, tid)
    reconcile, evaluate = [parser().parse_args(argv) for argv in job_queue.steps(root, task)]
    assert reconcile.run == run and reconcile.snapshot == snapshot and reconcile.direct_votes
    assert reconcile.routed_conflict_policy == "multiple"
    assert evaluate.run == reconcile.out == root / tid / "run" and evaluate.gt == gt
