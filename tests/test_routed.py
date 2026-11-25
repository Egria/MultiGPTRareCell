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
    cfg = config(override={"mode": ROUTED_MODE, "provider": {"model": "test", "api_key_env": None}})
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
        reviewed_seeds=4, name_searches=3, seed_searches=3, frozen=6, reviews=1
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
    assert summary["code_upgrade"] == dict(code_pair, after_turn=old_turn)
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
