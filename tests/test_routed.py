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
