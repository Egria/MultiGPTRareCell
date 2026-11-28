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


def test_direct_vote_challenges_respect_cycles_parents_and_stages():
    from copy import deepcopy
    from rarecell.single_choice import apply_direct_votes

    def ballot(left, right, winner, mode, cells):
        return dict(pair=f"{left}/{right}", outcomes=[[c, winner, mode, 0.75] for c in cells])

    groups = [dict(id=k, label=k, reliable=True, members=list("xyzuv") if k == "B" else [])
              for k in "ABCP"]
    result = dict(groups=groups, panels={k: dict(canonical=k) for k in "ABCP"},
                  parent_fallbacks=[["x", "P", "B", ["P"]]],
                  pair_votes=[ballot("A", "B", "A", "cross_fitted_modules", "xyz"),
                              ballot("B", "C", "B", "biology", "y"),
                              ballot("A", "C", "C", "biology", "y"),
                              ballot("B", "C", "B", "biology", "v"),
                              ballot("A", "C", "C", "biology", "z"),
                              ballot("P", "A", "P", "biology", "xz"),
                              ballot("P", "B", "P", "biology", "xz")],
                  unclaimed_conflicts=dict(pair_votes=[ballot("A", "B", "A", "biology", "u")]))
    result["panels"]["A"]["broader_labels"] = ["P"]
    shuffled = deepcopy(result)
    shuffled["pair_votes"].reverse()
    shuffled["groups"].reverse()
    for artifact in (result, shuffled):
        apply_direct_votes(artifact)
        members = {g["id"]: g["members"] for g in artifact["groups"]}
        assert members == {"A": ["u", "x"], "B": ["v", "y"], "C": ["z"], "P": []}
        assert artifact["parent_fallbacks"] == [["x", "P", "A", ["P"]]]
        assert [(c["cell_id"], c["stage"]) for c in artifact["ownership_changes"]] == [
            ("x", "nominated"), ("z", "nominated"), ("u", "unclaimed")]


@pytest.mark.parametrize("historical,policy", [(False, "single"), (True, "single"), (False, "multiple")])
def test_direct_vote_submission_and_offline_replay(tmp_path, monkeypatch, historical, policy):
    from rarecell import forced_routing
    from rarecell.forced_replay import replay

    cfg, snap, ann, _ = setup_input(tmp_path)
    cfg.update(routed_conflict_policy=policy, potential_names=["A", "B", "C"])
    monkeypatch.setattr("rarecell.run.verify_annotation", lambda *args, **kw: [])
    measure = forced_routing.measure

    def measured(x, cells, index, groups, specs, focus):
        pairs, cards = measure(x, cells, index, groups, specs, focus)
        # A defeats B on the panel. Geometry favors B over C and C over A.
        for pair in pairs:
            card = cards[pair["key"]]
            sides = (card["left"], card["right"])
            geometry = {("A", "B"): .1, ("A", "C"): .45, ("B", "C"): .95}
            pair["geometry"][:] = geometry[sides] if sides in geometry else 1 - geometry[sides[::-1]]
        return pairs, cards

    def reply(ctx, cfg):
        phase = ctx["phase"]
        if phase == "seed_triage":
            result = dict(action="route_seeds", reviews=[dict(candidate_id=r["candidate_id"],
                          reliability="unknown", reason="unknown", guide_label="") for r in ctx["seeds"]])
        elif phase == "seed_expansion":
            result = dict(action="freeze", rule=None, reason="abstain", support=0)
        elif phase == "blind_name_search":
            result = dict(action="freeze", rule={"any": ["G", "H"]}, reason="observed", support=80,
                          label_reliability="supported")
        elif phase == "conflict_panel_nomination":
            result = {"labels": {n: dict(canonical=n, panels=[["G"]]) for n in ctx["labels"]}}
        else:
            assert phase == "conflict_measured_review"
            rules = {}
            for key, card in ctx["pairs"].items():
                direction = "geometry"
                if {card["left"], card["right"]} == {"A", "B"}:
                    direction = "left" if card["left"] == "A" else "right"
                rules[key] = dict(default=direction, profiles={p: direction for p in card["profiles"]})
            result = {"rules": rules}
        return json.dumps(result), {"api_calls": 1, "total_tokens": 10}

    run = tmp_path / "run"
    monkeypatch.setattr(forced_routing, "measure", measured)
    with monkeypatch.context() as live:
        live.setattr(providers, "complete", reply)
        if historical:
            live.setattr(forced_routing, "apply_direct_votes", lambda result: None)
        summary = discover(snap, run, tmp_path / "cache", cfg, annotations=ann)
    assert summary["submitted"]
    source = read_json(run / summary["forced_choice"]["artifact_ref"])
    assert [g["label"] for g in source["groups"] if g["members"]] == ["B" if historical else "A"]
    if policy == "multiple":
        from rarecell.evaluation import evaluate

        assert summary["n_ambiguous_cells"] == 6 and summary["n_memberships"] == 12
        memberships = list(csv.DictReader((run / "memberships.tsv").read_text().splitlines(), delimiter="\t"))
        assert {r["proposed_label"] for r in memberships if r["is_primary"] == "1"} == {"A"}
        assert {r["proposed_label"] for r in memberships if r["is_primary"] == "0"} == {"C"}
        assert all(r["group_id"].startswith("vote/") for r in memberships if r["proposed_label"] == "C")
        ambiguous = list(csv.DictReader((run / "ambiguous_cells.tsv").read_text().splitlines(), delimiter="\t"))
        assert {r["cell_id"] for r in ambiguous} == set("abcdef")
        for row in ambiguous:
            assert json.loads(row["labels"]) == ["A", "C"]
            assert json.loads(row["group_ids"]) == [r["group_id"] for r in memberships if r["cell_id"] == row["cell_id"]]
        gt = tmp_path / "reference.tsv"
        gt.write_text("cell_id\ttarget_id\treference_kind\tstate\n" + "".join(
            c + "\tT\ttest\t" + state + "\n" for c,state in zip("abcdef", ["positive"]*3+["negative","unknown","disputed"])))
        score = evaluate(run, gt, tmp_path / "eval")
        assert score["submitted_memberships"]["per_target"][0]["tp"] == 3
        assert score["submitted_memberships"]["per_target"][0]["fp"] == 1
        assert (tmp_path / "eval/membership_targets.tsv").is_file()
        gt.write_text("cell_id\tlabel\n" + "".join(c+"\t"+("T" if c < "d" else "U")+"\n" for c in "abcdef"))
        scores = evaluate(run, gt, tmp_path / "fine")
        assert scores["submitted_memberships"]["ari"] is None
        assert scores["submitted_memberships"]["macro_f1"] > scores["submitted_partition"]["macro_f1"]
    else:
        assert not (run / "memberships.tsv").exists()
    for change in (False, True):
        out = tmp_path / str(change)
        replay(run, snap, None, None, out, cfg, direct_votes=change)
        if policy == "multiple":
            assert (out / "ambiguous_cells.tsv").is_file()
            with (out / "memberships.tsv").open() as stream:
                replayed = list(csv.DictReader(stream, delimiter="\t"))
            assert [(r["cell_id"], r["vote_group_id"], r["is_primary"]) for r in replayed] == [
                (r["cell_id"], r["vote_group_id"], r["is_primary"]) for r in memberships]
        new = read_json(out / "artifacts/forced_choice.json")
        label = "B" if historical and not change else "A"
        assert [(g["label"], g["members"]) for g in new["groups"] if g["members"]] == [(label, list("abcdef"))]
        assert new["pair_votes"] == source["pair_votes"]
        assert len(new.get("ownership_changes", [])) == (0 if label == "B" else 6)


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


def test_multiple_membership_evidence_precedence_and_cycles():
    from rarecell.overlap_membership import retained_groups
    labels = {g: g for g in "ABCP"}
    # Biology removes B. Its geometry win cannot then remove C.
    ballots = [("A", "B", "A", "biology", 0), ("B", "C", "B", "geometry", 1),
               ("A", "C", "A", "geometry", .99), ("P", "C", "P", "geometry", 1)]
    assert retained_groups("A", ballots, labels, {"C": {"P"}}) == {"A", "C"}
    ballots[2] = ("A", "C", "A", "geometry", 1)
    assert retained_groups("A", ballots, labels, {"C": {"P"}}) == {"A"}
    cycle = [("A", "B", "A", "biology", .5), ("B", "C", "B", "biology", .5),
             ("A", "C", "C", "cross_fitted_modules", .5)]
    assert retained_groups("A", cycle, labels, {}) == {"A", "B", "C"}
    cycle.append(("P", "C", "P", "biology", .5))
    assert retained_groups("P", cycle, labels, {}) == {"P"}


def test_multiple_membership_unclaimed_ballots_and_empty_ambiguity(tmp_path):
    from rarecell.overlap_membership import write_memberships
    groups = [dict(id="A", label="alpha", reliable=False, members=["c"]),
              dict(id="B", label="beta", reliable=False, members=[])]
    result = dict(groups=groups, panels={}, unclaimed_conflicts={"pair_votes": [
        {"pair": "A/B", "outcomes": [["c", "A", "stable_id_tie", .5]]}]})
    part = [dict(cell_id="c", group_id="claim/a", proposed_label="alpha", status="provisional"),
            dict(cell_id="u", group_id="residual/u", proposed_label="unknown", status="unresolved")]
    assert write_memberships(tmp_path, result, part)["n_memberships"] == 3
    result.pop("unclaimed_conflicts")
    assert write_memberships(tmp_path, result, part)["n_ambiguous_cells"] == 0
    assert len((tmp_path / "ambiguous_cells.tsv").read_text().splitlines()) == 1
