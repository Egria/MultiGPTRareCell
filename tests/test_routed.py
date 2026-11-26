import csv
import json
from contextlib import closing
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from rarecell import job_queue, providers
from rarecell.cli import parser
from rarecell.config import ROUTED_MODE, config, fingerprint, read_json, read_lines, write_json
from rarecell.discovery import Controller
from rarecell.finalize import assignments
from rarecell.routed import RoutedDiscovery
from rarecell.rule_search import ExpressionSearch
from rarecell.run import discover, resume


def setup_input(tmp_path):
    snap = tmp_path / "snapshot"
    snap.mkdir()
    cells = list("abcdef")
    (snap / "cells.tsv").write_text("cell_id\n" + "\n".join(cells) + "\n")
    (snap / "genes.tsv").write_text("gene_id\tsymbol\ng\tG\nh\tH\nu\tU\n")
    sparse.save_npz(
        snap / "data.npz",
        sparse.csr_matrix(
            [
                [4.0, 1.0, 0.0],
                [3.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 1.0, 1.0],
                [0.0, 1.0, 0.0],
                [0.0, 1.0, 1.0],
            ]
        ),
    )
    (snap / "baselines.tsv").write_text("view\tresolution\tcell_id\tgroup_id\n")
    seeds = [dict(id="S" + c, candidate_id=c, cell_id=c) for c in ("a", "c", "d", "f")]
    (snap / "seeds.jsonl").write_text("".join(json.dumps(s) + "\n" for s in seeds))
    (snap / "initial.jsonl").write_text(
        "".join(
            json.dumps(dict(id=s["candidate_id"], members=[s["cell_id"]], sources=[], markers=[])) + "\n"
            for s in seeds
        )
    )
    manifest = dict(
        schema="rarecell.snapshot.v1",
        data_version="data",
        views=[],
        status="complete",
        stages={},
        files={},
        tissue="pancreas",
        species="human",
    )
    manifest["snapshot_hash"] = fingerprint(manifest)
    write_json(snap / "manifest.json", manifest)
    annotation = tmp_path / "annotations"
    annotation.mkdir()
    write_json(annotation / "annotation_status.json", {"status": "success"})
    ann = [
        dict(candidate_id=cid, label=name, ontology_id=term, status="annotated", background_type=scope)
        for cid, name, term, scope in [
            ("a", "Glia", "CL:1", "global"),
            ("c", "Glial cell", "CL:1", "global"),
            ("d", "Other", "CL:2", "neighbor"),
        ]
    ]
    cfg = config(override={"mode": ROUTED_MODE, "routed_conflict_policy":"legacy", "provider": {"model": "test", "api_key_env": None}})
    return cfg, snap, annotation, ann


def test_public_route_resume_and_blinding(tmp_path, monkeypatch):
    cfg, snap, ann_dir, ann = setup_input(tmp_path)
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **kw: ann)
    blind = []

    def reply(context, cfg):
        phase = context["phase"]
        if phase == "seed_triage":
            req = dict(
                action="route_seeds",
                reviews=[
                    dict(
                        candidate_id=c["candidate_id"],
                        guide_label="Glia" if c["candidate_id"] == "a" else "",
                        reliability="supported" if c["candidate_id"] == "a" else "uncertain",
                        reason="seed assessment",
                    )
                    for c in context["seeds"]
                ],
            )
        elif phase == "blind_name_search":
            assert not set(context) & {"scopes", "seed_review", "annotations", "seed_markers", "drafts"}
            name = context["target"]
            if not context["history"]:
                blind.append(name)
                req = dict(action="survey", genes=["G"])
            else:
                assert len(context["history"]) == 1
                req = dict(
                    action="freeze",
                    rule=["G", 3.0] if name == "Glia" else None,
                    support=73,
                    reason="intensity",
                )
        elif phase == "seed_expansion":
            seed = context["seed_review"]["candidate_id"]
            req = dict(
                action="freeze",
                rule={"scope": "seed"},
                support=55,
                reason="retain core",
                proposed_label={"c": "Glia", "d": "New state", "f": "unknown"}[seed],
                kind="state" if seed == "d" else "lineage",
            )
        else:
            assert context["phase"] == "boundary_review"
            # One same-name pair with no shared cells.
            assert context["scopes"]["shared"] == 0
            req = dict(
                action="resolve",
                choice="merge",
                support=82,
                proposed_label="Glia",
                reason="matching program across both groups",
            )
        return json.dumps(req), {"api_calls": 1, "total_tokens": 12}

    monkeypatch.setattr(providers, "complete", reply)
    run = tmp_path / "run"
    paused = discover(snap, run, tmp_path / "cache", cfg, annotations=ann_dir, pause_after=2)
    assert paused["resumable"] and not paused["submitted"]
    result = resume(run)
    assert result["submitted"] and result["routing"] == dict(
        reviewed_seeds=4, name_searches=3, seed_searches=3, frozen=6, reviews=1, deferred_pairs=0
    )
    assert sorted(blind) == ["Glia", "New state", "Other"]
    assert not any(
        c["sources"][0]["kind"] == "ontology_seed_pool"
        for c in read_lines(run / "candidates.jsonl")
        if c["sources"]
    )
    with (run / "partition.tsv").open() as f:
        assignments_by_cell = {r["cell_id"]: r["proposed_label"] for r in csv.DictReader(f, delimiter="\t")}
    assert {c for c, name in assignments_by_cell.items() if name == "Glia"} == set("abc")
    claims = read_lines(run / "route_claims.jsonl")
    assert any(r["kind"] == "state" and r["n"] == 1 and r["support"] == 55 for r in claims)
    assert "route_resolutions.jsonl" in read_json(run / "frozen.json")["files"]
    assert result["n_pending_seeds"] == 0


def test_reconciliation_flanks_and_unknown_exact_conflict(tmp_path):
    cfg, snap, _, ann = setup_input(tmp_path)
    ctrl = Controller(snap, tmp_path / "run", tmp_path / "cache", cfg, ROUTED_MODE, ann)
    with closing(ctrl.state):
        ctrl.initialize()
        flow = RoutedDiscovery(ctrl)

        def draft(cells, name):
            return flow.draft(
                cells, ["a"], name, "lineage", {"reason": "test"}, ctrl.candidate("a")["evidence_ids"]
            )

        def split(a, b, left, right):
            job = dict(id="pair/" + a + "/" + b, kind="pair", parents=[a, b])
            scopes = flow.scope_members(job)
            return flow.resolve(
                job,
                {
                    "choice": "split_overlap",
                    "left_support": 76,
                    "right_support": 58,
                    "reason": "measure exclusive and contested cells",
                },
                scopes,
                {"left": set(left), "right": set(right)},
                "turn",
            )["retained"]

        a, b = draft("abc", "A"), draft("bcde", "B")
        resolved = split(a, b, "b", "")
        assert [set(ctrl.candidate(i)["members"]) for i in resolved] == [set("ab"), set("de")]
        c, d = draft("ef", "C"), draft("ef", "D")
        (unknown,) = split(c, d, "", "")
        assert ctrl.candidate(unknown)["proposed_label"] == "unknown"
        assert set(ctrl.candidate(unknown)["members"]) == set("ef")
        expr = ExpressionSearch(ctrl)
        with pytest.raises(ValueError):
            expr.mask({"scope": "seed"}, {})
        assert np.flatnonzero(expr.mask({"k": [2, ["g", "h", "u"]]}, {})).tolist() == [0, 2, 3, 5]


def test_routed_gui_budget_override(tmp_path):
    root = tmp_path / "tasks"
    job_queue.initialize(root)
    data = Path(__file__).resolve().parents[1] / "examples/minimal_sparse/input/dataset.yaml"
    jid = job_queue.submit(
        root,
        tmp_path,
        dict(kind="pipeline", mode=ROUTED_MODE, dataset=str(data), model="test", annotation_model="test"),
    )
    steps = list(job_queue.steps(root, job_queue.get_task(root, jid)))
    parsed = parser().parse_args(steps[-1])
    cfg = config(parsed.config)
    assert [s[0] for s in steps] == ["prepare", "annotate", "discover"]
    assert cfg["mode"] == ROUTED_MODE and cfg["budget"]["max_agent_turns"] == 800
    assert (
        config(override={"mode": ROUTED_MODE, "budget": {"max_agent_turns": 3}})["budget"]["max_agent_turns"]
        == 3
    )


def test_resume_legacy_depth_crash_without_repeating_model_turns(tmp_path, monkeypatch):
    from rarecell import run as runner
    from rarecell.config import RESOURCES
    from rarecell.evaluation import evaluate
    from rarecell.state import State

    cfg, snap, ann_dir, ann = setup_input(tmp_path)
    cfg["potential_names"] = ["Shared"]
    monkeypatch.setattr(runner, "verify_annotation", lambda *args, **kwargs: ann)
    calls = []

    def reply(ctx, cfg):
        phase = ctx["phase"]
        if phase == "seed_triage":
            req = dict(action="route_seeds", reviews=[
                dict(candidate_id=seed["candidate_id"], guide_label="Shared", reliability="uncertain",
                     reason="Tentative singleton identity") for seed in ctx["seeds"]])
            calls.append((phase, None))
        elif phase == "blind_name_search":
            # Exactly seed a matches, while four seeds nominated this name.
            req = dict(action="freeze", rule=["G", 3.5] if ctx["target"] == "Shared" else None,
                       support=60, reason="Expression rule")
            calls.append((phase, ctx["target"]))
        else:
            assert phase == "seed_expansion" and "retrieved" in ctx["scopes"]
            req = dict(action="freeze", rule=None, support=0, reason="Abstain")
            calls.append((phase, ctx["seed_review"]["candidate_id"]))
        return json.dumps(req), {"api_calls": 1, "total_tokens": 12}

    monkeypatch.setattr(providers, "complete", reply)
    code_pair = read_json(RESOURCES / "routed_checkpoint_upgrade.json")
    run, cache = tmp_path / "run", tmp_path / "cache"
    # Reproduce the two PR24 depth checks and its unhandled plain ValueError.
    with monkeypatch.context() as legacy:
        legacy.setattr(runner, "code_fingerprint", lambda: code_pair["from"])
        legacy.setattr(Controller, "expansion_depth",
                       lambda self, c, roots, depths: max(depths.get(r, 0) for r in roots))
        legacy.setattr("rarecell.discovery.ExpansionDepthExhausted", ValueError)
        with pytest.raises(ValueError, match="Expansion depth exhausted"):
            discover(snap, run, cache, cfg, annotations=ann_dir)
    assert ("seed_expansion", "a") in calls and ("seed_expansion", "c") not in calls
    assert not (run / "frozen.json").exists()
    with pytest.raises(FileNotFoundError, match="resume discovery"):
        evaluate(run, tmp_path / "gt.tsv", tmp_path / "scores")
    with closing(State(run, cache)) as state:
        messages = state.all("Message")
        before = {row["message_id"]: row["context_hash"] for row in messages}
        old_turn = state.get_meta("finished_turn")
        assert state.get_meta("depths")["c"] == 1
        assert not any(r["candidate_id"] == "c" for r in state.all("Investigation"))
    with monkeypatch.context() as unrecognized:
        unrecognized.setattr(runner, "code_fingerprint", lambda: "unrecognized")
        with pytest.raises(ValueError, match="method version changed"):
            resume(run)

    # Exercise this approved version pair even when a future, separately pinned
    # implementation runs the regression suite.
    with monkeypatch.context() as release:
        release.setattr(runner, "code_fingerprint", lambda: code_pair["to"])
        summary = resume(run)
    assert summary["submitted"] and summary["status"] == "complete"
    assert summary["code_upgrade"] == {"from": code_pair["from"], "to": code_pair["to"], "after_turn": old_turn}
    assert read_json(run / "checkpoint.json")["code_sha256"] == code_pair["to"]
    assert sorted(label for phase, label in calls if phase == "seed_expansion") == ["a", "c", "d", "f"]
    assert len(calls) == len(set(calls))
    with closing(State(run, cache)) as state:
        assert before == {m["message_id"]: m["context_hash"] for m in state.all("Message") if m["message_id"] in before}
        assert len(state.all("ProviderDispatch")) == len(calls)
        assert len(state.all("Investigation")) == 4
    assert (run / "frozen.json").is_file()
    (tmp_path / "gt.tsv").write_text("cell_id\tlabel\n" + "".join(f"{c}\tstate\n" for c in "abcdef"))
    assert evaluate(run, tmp_path / "gt.tsv", tmp_path / "scores")["submitted_partition"]["n_cells"] == 6


def test_routed_zero_depth_keeps_seed_available_for_expression(tmp_path, monkeypatch):
    cfg, snap, ann_dir, _ = setup_input(tmp_path)
    cfg["budget"]["max_expansion_depth"] = 0
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *args, **kwargs: [])
    reviewed = []

    def reply(ctx, cfg):
        if ctx["phase"] == "seed_triage":
            req = dict(action="route_seeds", reviews=[
                dict(candidate_id=s["candidate_id"], reliability="unknown", reason="Unlabeled")
                for s in ctx["seeds"]])
        else:
            assert ctx["phase"] == "seed_expansion" and ctx["scopes"] == {"seed": 1}
            assert ctx["retrieval"]["status"] == "depth_exhausted"
            reviewed.append(ctx["seed_review"]["candidate_id"])
            req = dict(action="freeze", rule={"scope": "seed"}, support=20, reason="Measured core only")
        return json.dumps(req), {"api_calls": 1, "total_tokens": 10}

    monkeypatch.setattr(providers, "complete", reply)
    run, cache = tmp_path / "run", tmp_path / "cache"
    summary = discover(snap, run, cache, cfg, annotations=ann_dir)
    assert summary["submitted"] and sorted(reviewed) == ["a", "c", "d", "f"]
    assert summary["n_expanded_seeds"] == 0 and (run / "frozen.json").exists()


def overlap_reply(ctx):
    if ctx["phase"] == "seed_triage":
        return dict(action="route_seeds", reviews=[
            dict(candidate_id=s["candidate_id"], guide_label="L", reliability="supported", reason="Measured seed")
            for s in ctx["seeds"]])
    if ctx["phase"] == "blind_name_search":
        return dict(action="freeze", rule={"L": ["G", 3], "M": "H"}[ctx["target"]],
                    support=70, reason="Freeze exact expression boundary")
    raise AssertionError(ctx["phase"])


def test_extend_query_then_turn_budget_replays_paid_response(tmp_path, monkeypatch):
    from rarecell import run as runner
    from rarecell.config import RESOURCES, hash_file
    from rarecell.state import State

    cfg, snapshot, annotations, _ = setup_input(tmp_path)
    cfg["potential_names"] = ["L", "M"]
    cfg["budget"].update(max_new_queries=1, max_agent_turns=3)
    monkeypatch.setattr(runner, "verify_annotation", lambda *a, **kw: [])
    dispatches = []

    def reply(ctx, cfg):
        dispatches.append((ctx["phase"], ctx.get("target")))
        if ctx["phase"] == "boundary_review":
            action = dict(action="resolve", choice="left", reason="Retain measured boundary")
        elif ctx.get("target") == "M" and not ctx["history"]:
            action = dict(action="survey", genes=["H"])
        else:
            action = overlap_reply(ctx)
        return json.dumps(action), {"api_calls": 1, "total_tokens": 10}

    monkeypatch.setattr(providers, "complete", reply)
    run, cache = tmp_path / "run", tmp_path / "cache"
    compatible = read_json(RESOURCES / "routed_checkpoint_upgrade.json")
    with monkeypatch.context() as old:
        old.setattr(runner, "code_fingerprint", lambda: compatible["also_from"][-1])
        before = discover(snapshot, run, cache, cfg, annotations=annotations)
    assert before["detail"] == "routed_query_budget" and not before["submitted"]
    assert resume(run) == before
    with pytest.raises(ValueError, match="exhausted cumulative ceiling"):
        resume(run, max_agent_turns=8)
    original_freeze = hash_file(run / "frozen.json")
    with closing(State(run, cache)) as state:
        messages = state.all("Message")
        claim, = state.all("RouteFreeze")
        assert state.get_meta("finished_turn") == 2 and len(messages) == 3
    queue = tmp_path / "tasks"
    job_queue.initialize(queue)
    task = job_queue.submit(queue, tmp_path, dict(kind="resume", run=str(run), max_new_queries="10"))
    args = parser().parse_args(next(job_queue.steps(queue, job_queue.get_task(queue, task))))
    assert args.cmd == "resume" and args.max_new_queries == 10 and args.max_agent_turns is None
    # Simulate a killed process just after the increased budget has been saved.
    with monkeypatch.context() as crash:
        def interrupt(*args):
            raise SystemExit
        crash.setattr(RoutedDiscovery, "apply", interrupt)
        with pytest.raises(SystemExit):
            resume(run, max_new_queries=10)
    middle = resume(run)
    assert middle["detail"] == "routed_turn_budget" and len(dispatches) == 3
    assert hash_file(run / "history" / original_freeze / "frozen.json") == original_freeze
    grant, = middle["budget_extensions"]
    assert grant == dict(source_freeze_sha256=original_freeze, previous_limits=dict(max_new_queries=1, max_agent_turns=3),
                         limits=dict(max_new_queries=10, max_agent_turns=3), after_turn=2, used_queries=1)
    with pytest.raises(ValueError, match="at least as high"):
        resume(run, max_new_queries=9, max_agent_turns=8)
    after = resume(run, max_agent_turns=8)
    assert after["submitted"] and len(dispatches) == 5
    assert len(after["budget_extensions"]) == 2
    assert after["budget_extensions"][-1]["used_queries"] == 2
    assert after["budget_extensions"][-1]["limits"] == dict(max_new_queries=10, max_agent_turns=8)
    with closing(State(run, cache)) as state:
        assert state.all("Message")[:len(messages)] == messages
        assert claim in state.all("RouteFreeze")
    assert resume(run, max_new_queries=10, max_agent_turns=8) == after


@pytest.mark.parametrize("stop", ["bad_decisions", "turn_budget"])
def test_deferred_pair_uses_normal_submission_policy(tmp_path, monkeypatch, stop):
    cfg, snap, ann_dir, _ = setup_input(tmp_path)
    cfg["potential_names"] = ["L", "M"]
    if stop == "turn_budget":
        cfg["budget"]["max_agent_turns"] = 3
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **k: [])
    pair_calls = []

    def reply(ctx, cfg):
        if ctx["phase"] == "boundary_review":
            pair_calls.append(ctx)
            req = [dict(action="resolve", choice=side, reason="Conflicting decisions") for side in ["left", "right"]]
        else:
            req = overlap_reply(ctx)
        return json.dumps(req), {"api_calls": 1, "total_tokens": 10}

    monkeypatch.setattr(providers, "complete", reply)
    run = tmp_path / "run"
    out = discover(snap, run, tmp_path / "cache", cfg, annotations=ann_dir)
    assert out["submitted"] and out["submission_status"] == "provisional"
    assert out["submission_source"] == "controller_routed" and out["detail"] == "submitted_with_deferred_pairs"
    assert out["routing"]["reviews"] == 0 and out["routing"]["deferred_pairs"] == 1
    assert len(pair_calls) == (0 if stop == "turn_budget" else cfg["budget"]["stall_turns"])
    resolution, = read_lines(run / "route_resolutions.jsonl")
    assert resolution["source"] == "controller" and len(resolution["retained"]) == 2
    with (run / "partition.tsv").open() as stream:
        cells = {r["cell_id"]: r for r in csv.DictReader(stream, delimiter="\t")}
    assert cells["b"]["proposed_label"] == "L" and cells["a"]["status"] == "unresolved"
    assert all(cells[c]["proposed_label"] == "M" for c in "cdef")


def test_resume_blocked_pair_preserves_freeze_and_old_decisions(tmp_path, monkeypatch):
    from rarecell import run as runner
    from rarecell.config import RESOURCES, hash_file
    from rarecell.evaluation import evaluate
    from rarecell.state import State

    cfg, snap, ann_dir, _ = setup_input(tmp_path)
    cfg["potential_names"] = ["L", "M"]
    monkeypatch.setattr(runner, "verify_annotation", lambda *a, **k: [])
    calls = []
    def reply(ctx, cfg):
        calls.append(ctx["phase"])
        if ctx["phase"] != "boundary_review":
            req = overlap_reply(ctx)
        elif not ctx["history"]:
            req = {"actions": [dict(operation="survey", arguments={"genes":["G"]}), dict(action="views", rules=[{"scope":"shared"}], genes=["H"]),
                               dict(action="resolve", choice="left", reason="Premature decision")]}
        else:
            assert len(ctx["history"]) == 2 and ctx["history"][-1]["not_executed"] == ["resolve"]
            req = dict(action="resolve", choice="right", reason="Use measured flank")
        return json.dumps(req), {"api_calls": 1, "total_tokens": 10}

    monkeypatch.setattr(providers, "complete", reply)
    code = read_json(RESOURCES / "routed_checkpoint_upgrade.json")
    old = code["also_from"][0]
    run, cache = tmp_path / "run", tmp_path / "cache"
    # PR25 persisted this failure as a terminal partial export. Use the same
    # serialization path, then exercise the public resume command against it.
    def old_apply(self, job, ctx, scopes, req, mid):
        if job and job["kind"] == "pair":
            raise ValueError("Return one JSON action object")
        return normal_apply(self, job, ctx, scopes, req, mid)
    normal_apply = RoutedDiscovery.apply
    with monkeypatch.context() as legacy:
        legacy.setattr(runner, "code_fingerprint", lambda: old)
        legacy.setattr(RoutedDiscovery, "apply", old_apply)
        old_loop = RoutedDiscovery.run
        def legacy_loop(self, obtain, pause_after, retry_unknown):
            # Pause just after a rejected response, then export using the old terminal detail.
            status, _ = old_loop(self, obtain, 4, retry_unknown)
            failed = self.s.get_meta("route_last_error")
            self.s.set_meta("current_issue", "Routed task " + failed["task"] + ": " + failed["error"])
            return status, "routed_action_blocked"
        legacy.setattr(RoutedDiscovery, "run", legacy_loop)
        before = discover(snap, run, cache, cfg, annotations=ann_dir)
    assert not before["submitted"] and not before["resumable"]
    digest = hash_file(run / "frozen.json")
    # Record the same earlier depth upgrade present on DCC.
    previous = {"from":code["from"], "to":old, "after_turn":1}
    checkpoint = read_json(run / "checkpoint.json")
    write_json(run / "checkpoint.json", dict(checkpoint, code_upgrade=previous))
    with closing(State(run, cache)) as state:
        state.set_meta("code_upgrade", previous)
        old_messages = state.all("Message")
    old_count = len(calls)
    with monkeypatch.context() as release:
        release.setattr(runner, "code_fingerprint", lambda: code["to"])
        after = resume(run)
    assert calls[old_count:] == ["boundary_review", "boundary_review"]
    assert after["submitted"] and after["routing"]["reviews"] == 1 and after["routing"]["deferred_pairs"] == 0
    assert after["code_upgrades"] == [previous, {"from":old, "to":code["to"], "after_turn":4}]
    assert hash_file(run / "history" / digest / "frozen.json") == digest
    with closing(State(run, cache)) as state:
        assert state.all("Message")[:len(old_messages)] == old_messages
        assert len(state.all("RouteExchange")) == 2 and len(state.all("RouteResolution")) == 1
    gt = tmp_path / "gt.tsv"
    gt.write_text("cell_id\tlabel\n" + "".join(f"{c}\tX\n" for c in "abcdef"))
    assert evaluate(run, gt, tmp_path / "scores")["submitted_partition"]["n_cells"] == 6
