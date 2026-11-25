import csv
import json

from scipy import sparse

from rarecell import providers
from rarecell.config import config, read_json
from rarecell.run import discover, resume
from rarecell.single_choice import build_groups, choose, measure
from test_routed import setup_input


def test_anonymous_three_way_preserves_union_and_ties():
    cells = list("abcd")
    claims = [dict(candidate_id=k, label="unknown", members=list("abc")) for k in "ZYX"]
    ids = {cell: i for i, cell in enumerate(cells)}
    groups = build_groups(claims, {}, {}, ids)
    pairs, cards = measure(sparse.eye(4, format="csr"), cells, {}, groups, {})
    owners, votes, _ = choose(cells, groups, pairs, {})
    assert not cards and all(p["route"] == "blind" for p in votes)
    assert owners == {c: min(g["id"] for g in groups) for c in "abc"}
    assert choose(cells, groups[::-1], pairs[::-1], {})[0] == owners
    assert all(p["counts"] == {"stable_id_tie": 3} for p in votes)


def test_live_submission_reliable_and_uncertain_names(tmp_path, monkeypatch):
    cfg, snap, annotations, _ = setup_input(tmp_path)
    assert config()["routed_conflict_policy"] == "single"
    cfg["routed_conflict_policy"] = "single"
    cfg["potential_names"] = ["A", "B", "C"]
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **kw: [])
    phases = []

    def reply(ctx, cfg):
        phase = ctx["phase"]
        phases.append(phase)
        if phase == "seed_triage":
            result = dict(action="route_seeds", reviews=[
                dict(candidate_id=row["candidate_id"], reliability="uncertain",
                     reason="identity unresolved", guide_label="")
                for row in ctx["seeds"]])
        elif phase == "blind_name_search":
            name = ctx["target"]
            result = dict(action="freeze", rule={"A": "G", "B": "H", "C": "U"}[name],
                          support=70, reason="expressing cells", label_reliability="supported")
            if name == "C":
                del result["label_reliability"]
        elif phase == "seed_expansion":
            result = dict(action="freeze", rule=None, reason="abstain", support=0)
        elif phase == "conflict_label_audit":
            result = {"reviews": {key: "uncertain" for key in ctx["candidates"]}}
        elif phase == "conflict_panel_nomination":
            assert set(ctx) == {"phase", "labels", "vocabulary", "species", "tissue", "instruction"}
            assert set(ctx["labels"]) == {"A", "B"}
            result = {"labels": {k: dict(canonical=k, panels=[[{"A":"G","B":"H"}[k]]])
                                 for k in ctx["labels"]}}
            result["labels"]["A"]["versus"] = {"B": {"panels": [["U"]]}}
        else:
            assert phase == "conflict_measured_review"
            assert all(p["pair_specific_panels"]["A"]["panels"] == [["U"]] for p in ctx["pairs"].values())
            assert all({p["left"], p["right"]} == {"A", "B"} for p in ctx["pairs"].values())
            result = {"rules": {k: dict(default="geometry") for k in ctx["pairs"]}}
        return json.dumps(result), {"api_calls": 1, "total_tokens": 10}

    monkeypatch.setattr(providers, "complete", reply)
    run = tmp_path / "run"
    # 1 batch of triage, 7 searches, followed by the label audit and panel nomination.
    result = discover(snap, run, tmp_path / "cache", cfg, annotations=annotations, pause_after=10)
    assert result["resumable"] and not result["submitted"]
    result = resume(run)
    assert result["submitted"] and result["forced_choice"]["n_claimed"] == 6
    assert phases.count("conflict_panel_nomination") == phases.count("conflict_measured_review") == 1
    audit = read_json(run / result["forced_choice"]["artifact_ref"])
    assert {p["route"] for p in audit["pair_votes"]} == {"blind", "named"}
    rows = list(csv.DictReader((run / "partition.tsv").open(), delimiter="\t"))
    assert len({r["cell_id"] for r in rows}) == len(rows) == 6
    assert all(r["group_id"].startswith("claim/") for r in rows)
