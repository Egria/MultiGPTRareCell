import csv
import json

import pytest
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


@pytest.mark.parametrize("extra_field", [None, "action", "metadata"], ids=["complete", "cached-action", "cached-metadata"])
def test_review_unclaimed_across_routes_resumes_and_replays(tmp_path, monkeypatch, extra_field):
    from rarecell.forced_replay import replay
    from rarecell.config import RESOURCES, fingerprint
    from rarecell import run as runner
    cfg, snap, ann, _ = setup_input(tmp_path)
    cfg.update(routed_conflict_policy="single", potential_names=["A", "B"])
    cfg["budget"]["max_agent_turns"] = 40
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *a, **k: [])
    phases = []

    def reply(ctx, settings):
        phase = ctx["phase"]
        phases.append(phase)
        if phase == "seed_triage":
            result = dict(action="route_seeds", reviews=[dict(candidate_id=row["candidate_id"],
                reliability="unknown", reason="unknown identity", guide_label="") for row in ctx["seeds"]])
        elif phase == "blind_name_search":
            a = ctx["target"] == "A"
            if not ctx["history"]:
                result = dict(action="views", rules=["G" if a else "H"], genes=[])
            else:
                result = dict(action="freeze", rule=["G",4] if a else "U", support=80,
                              reason="freeze the measured core", label_reliability="supported")
        elif phase == "seed_expansion":
            result = dict(action="freeze", rule=None, reason="abstain", support=0)
        elif phase == "conflict_panel_nomination":
            result = {"labels":{n:dict(canonical=n, panels=[["G" if n=="A" else "U"]],
                                       lineage_family=n, broader_labels=[]) for n in ctx["labels"]}}
        elif phase == "unclaimed_evidence_review":
            result = {"reviews":{k:dict(add="all" if card["label"]=="A" else ["c"],
                                        reason="Reviewed prior rule and outside hits") for k,card in ctx["options"].items()}}
            if extra_field and len(ctx["options"]) == 2:
                result["reviews"] = {k: v for k, v in result["reviews"].items() if ctx["options"][k]["label"] == "A"}
            elif extra_field:
                # Both forms sort before "reviews" when SQLite encodes the accepted response.
                key, = ctx["options"]
                result[extra_field] = ("review" if extra_field == "action" else
                                       {key: dict(add=[], reason="A different, unvalidated decision")})
        else:
            assert phase == "unclaimed_conflict_measured_review"
            assert all(v["shared"] == 1 for v in ctx["pairs"].values())
            result = {"rules": {k:dict(default="geometry") for k in ctx["pairs"]}}
        return json.dumps(result), {"api_calls":1,"total_tokens":20}

    monkeypatch.setattr(providers,"complete",reply)
    run=tmp_path/"run"
    # Pause after the accepted addition review and before arbitrating its shared cell.
    cutoff = 13 if extra_field else 11
    with monkeypatch.context() as old:
        if extra_field:
            upgrade = read_json(RESOURCES / "routed_checkpoint_upgrade.json")
            old.setattr(runner, "code_fingerprint", lambda: upgrade["also_from"][-1])
        a = discover(snap, run, tmp_path / "cache", cfg, annotations=ann, pause_after=cutoff)
    assert a["resumable"] and phases[-1] == "unclaimed_evidence_review"
    earlier_phases = len(phases)
    b = resume(run)
    assert b["submitted"] and phases[earlier_phases:] == ["unclaimed_conflict_measured_review"]
    assert phases.count("unclaimed_evidence_review") == (3 if extra_field else 1)
    if extra_field:
        assert b["code_upgrade"]["after_turn"] == cutoff
        assert b["forced_choice"]["review_fallback_items"] == 0
    forced=b["forced_choice"]
    assert (forced["unclaimed_options"],forced["unclaimed_accepted"],forced["unclaimed_competing"])==(2,2,1)
    art=read_json(run/forced["artifact_ref"])
    assert set(c for g in art["groups"] for c in g["members"])==set("abcdf")
    assert sum(len(g["members"]) for g in art["groups"])==5
    assert sum("c" in g["members"] for g in art["groups"])==1
    assert art["unclaimed_conflicts"]["pair_votes"][0]["outcomes"][0][0]=="c"
    out=tmp_path/"replayed"
    replay(run,snap,None,None,out,cfg)
    new=read_json(out/"artifacts/forced_choice.json")
    assert new["groups"]==art["groups"]
    assert new["unclaimed_reviews"]==art["unclaimed_reviews"]
    assert fingerprint(new["unclaimed_conflicts"])==fingerprint(art["unclaimed_conflicts"])


def test_parent_steps_aside_and_remaining_competitors_revote():
    import numpy as np
    from rarecell.lineage_review import relationships
    import pytest
    groups=[dict(id=k,label=k,reliable=(k!="unrelated"),members=np.array([0])) for k in ["broad","child","unrelated"]]
    def pair(a,b,geo):
        return dict(key=a+"/"+b,groups=(a,b),rows=[0],geometry=[geo],named=False)
    pairs=[pair("broad","child",.9),pair("broad","unrelated",.9),pair("child","unrelated",.2)]
    audit=[]
    parents=relationships({"broad":{"canonical":"broad"},"child":{"canonical":"child","broader_labels":["broad"]}})
    assert choose(["cell"],groups,pairs,{})[0]["cell"]=="broad"
    assert choose(["cell"],groups,pairs,{},parents,audit)[0]["cell"]=="unrelated"
    assert audit==[["cell","broad","unrelated",["broad"]]]
    with pytest.raises(ValueError,match="acyclic"):
        relationships({k:{"canonical":k,"broader_labels":[o]} for k,o in [("a","b"),("b","a")]})


def test_cross_fitted_marker_votes_require_reliable_separate_families():
    import numpy as np
    cells=[f"c{i}" for i in range(25)]
    x=sparse.csr_matrix([[5,4,0,0]]*13+[[0,0,4,5]]*12)
    groups=[dict(id="A",label="A",reliable=True,members=np.arange(13)),
            dict(id="B",label="B",reliable=True,members=np.arange(12,25))]
    specs={k:dict(canonical=k,lineage_family=k,panels=[[k+"1"],[k+"2"]]) for k in ["A","B"]}
    idx={n:i for i,n in enumerate(["A1","A2","B1","B2"])}
    pairs,cards=measure(x,cells,idx,groups,specs)
    rules={key:dict(default="right",profiles={p:"right" for p in card["profiles"]}) for key,card in cards.items()}
    owners,audit,_=choose(cells,groups,pairs,rules)
    assert owners["c12"]=="A" and audit[0]["counts"]=={"cross_fitted_modules":1}
    specs["B"]["lineage_family"]="A"
    assert "module_contrast" not in measure(x,cells,idx,groups,specs)[0][0]
    groups[0]["reliable"]=False
    assert not measure(x,cells,idx,groups,specs)[1]


def test_adaptive_can_recover_unclaimed_measured_hits_and_keep_other_owners(tmp_path):
    from rarecell.discovery import Controller
    from rarecell.protocol import parse_message
    from contextlib import closing
    cfg,snapshot,_,_=setup_input(tmp_path)
    ctrl=Controller(snapshot,tmp_path/"run",tmp_path/"cache",cfg,"A3_adaptive")
    def act(action,**args):
        cmd,=parse_message(json.dumps(dict(action=action,**args)))
        return ctrl.execute(cmd,"test")
    with closing(ctrl.state):
        ctrl.initialize()
        ctrl.state.set_meta("selected",["a","c"])
        hit=act("coexpress",scope="all",genes=["G"])
        offer=next(r for r in ctrl.context(1)["unclaimed_evidence"] if r["candidate_id"]=="a")
        assert offer["core_hits"]==offer["n_unclaimed"]==offer["other_claimed_hits"]==1
        assert offer["examples"]==["b"]
        act("propose_revision", parents=["a"],evidence_ids=[hit["evidence_id"]],groups=[dict(candidate_ids=["a"],matched_evidence_id=hit["evidence_id"],
            unclaimed_only=True,proposed_label="provisional")],reason="Include the unclaimed measured hit")
        chosen=ctrl.state.get_meta("selected")
        assert "c" in chosen
        new,=set(chosen)-{"c"}
        assert set(ctrl.candidate(new)["members"])=={"a","b"}


def test_nb_high_needs_counts_and_reports_an_unsupported_split(tmp_path):
    import numpy as np
    import pytest
    from types import SimpleNamespace
    from rarecell.discovery import rows
    from rarecell.rule_search import ExpressionSearch
    _,snap,_,_=setup_input(tmp_path)
    data=SimpleNamespace(snapshot=snap,cells=list("abcdef"),genes=rows(snap/"genes.tsv"))
    e=ExpressionSearch(data)
    with pytest.raises(ValueError,match="raw counts"):
        e.mask({"nb_high":"G"},{})
    sparse.save_npz(snap/"raw.npz",sparse.csr_matrix(np.array([[0,0,0]]+[[0,1,0]]*5)))
    assert e.mask({"nb_high":"G"},{}).tolist()==[False]*6
    assert e.nb_fits[e.index["G"]][0]["status"]=="insufficient_signal"
