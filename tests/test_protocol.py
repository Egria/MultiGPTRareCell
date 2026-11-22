"""Scientific equivalence, data identity and fail-local boundaries."""
import hashlib
import json

import pytest

from rarecell.protocol import MAX_ACTIONS, ProtocolError, parse_message, resolve_genes


def one(raw, **kwargs):
    value = raw if isinstance(raw, str) else json.dumps(raw, ensure_ascii=False)
    (item,) = parse_message(value, **kwargs)
    return item


def canonical_for_test(parsed):
    # A real query also binds the immutable cell set, data and transform versions.
    obj = {"action": parsed["action"], "args": parsed["args"]}
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def test_equivalent_aliases_shapes_and_repair():
    equivalent = [
        '{"action":"coexpress","genes":["CD3D","TRAC"],"candidate_id":"H017"}',
        '请计算。\n```json\n{"operation":"check co-expression","gene_list":"CD3D; TRAC","candidate":"H017",}\n```\n这会比较双阳性细胞。',
        '{"操作":"检查共表达","参数":{"基因":["TRAC","CD3D"],"候选":["H017"],"逻辑":"且","阈值":"0e2"}}',
        '共表达\n候选：H017\n基因：CD3D；TRAC\n逻辑：AND',
        '{"action":"coexpress","args":{"genes":["TRAC","CD3D"],"threshold":"0.","layer":"log1p"},'
        '"candidate":"H017","gene_list":"CD3D;TRAC","comment":"retain this as commentary"}',
    ]
    parsed = [one(raw) for raw in equivalent]
    assert len({canonical_for_test(p) for p in parsed}) == 1
    assert parsed[0]["args"] == {"genes": ["CD3D", "TRAC"], "candidate_id": "H017",
                                 "logic": "AND", "threshold": 0.0}
    assert parsed[1]["metadata"]["format_repair"] == "trailing_comma"
    assert parsed[3]["metadata"]["format_repair"] == "field_lines"
    assert parsed[4]["metadata"]["comment"].endswith("commentary")


def test_semantic_changes_not_merged():
    req = {"action": "coexpress", "genes": ["CD3D", "TRAC"], "members": ["CELL-O1"]}
    old = one(req)
    for replacement in ({"genes": ["CD3D", "TRA0"]}, {"members": ["CELL-01"]},
                        {"logic": "OR"}, {"threshold": "0.01"}):
        changed = one({**req, **replacement})
        assert canonical_for_test(old) != canonical_for_test(changed)
    compare = {"action": "compare_groups", "members": ["A"], "background": ["B"]}
    assert one(compare)["args"] != one({**compare, "background": ["C"]})["args"]


@pytest.mark.parametrize("extra", [
    {"genes": None}, {"genes": []}, {"threshold": None}, {"threshold": "unknown"},
    {"threshold": "1e400"}, {"threshold": "1e-4000"}, {"threshold": False},
    {"threshold": -1}, {"scope": "这些细胞"}, {"scope": "all"}, {"members": ["C"]},
    {"sample_filter": "donor A"}, {"metadata": {"scope": "all"}}, {"direction": "negative"},
    {"weights": [0.2]}, {"where": {"analysis_set": "only T cells"}}, {"layer": "counts"},
    {"candidate": "H018"}, {"evidence_purpose": "independent_validation"},
    {"max_expansion_depth": 99}, {"logic": "XOR"}, {"threshold": [0, 1]},
])
def test_ambiguity_and_conditions_block_only_affected_action(extra):
    bad = {"action": "coexpress", "genes": "A", "candidate_id": "H017", **extra}
    raw = json.dumps({"actions": [bad, {"action": "list_candidates", "limit": "2.0"}]})
    blocked, ok = parse_message(raw)
    assert blocked["status"] == "needs_clarification"
    assert ok["args"] == {"limit": 2}


def test_no_scope_defaults_when_the_source_is_ambiguous():
    assert "error" in one({"action": "coexpress", "genes": ["A"]})
    assert one({"action": "coexpress", "genes": ["A"], "scope": "全部细胞"})["args"]["scope"] == "all"
    assert "error" in one({"action": "compare_groups", "candidate": "a",
                           "background_id": "b", "background_members": ["one"]})


def test_ids_versions_and_strings_are_exact():
    req = {"operation": "retrieve_cells", "candidate": "h 01-A", "data_version": "v.1-A",
           "run_id": "R", "task_id": "t", "contract_version": "C-0", "message_id": "m1"}
    parsed = one(req)
    assert parsed["action"] == "expand_candidate"
    assert parsed["args"] == {"candidate_id": "h 01-A", "data_version": "v.1-A",
                              "run_id": "R", "task_id": "t", "contract_version": "C-0",
                              "message_id": "m1"}
    assert "error" in one({**req, "args": {"data_version": "v.1-a"}})
    for evidence in ("../evaluation/file", "E-O1", "E-01"):
        # This parser only returns literal IDs. A controller resolves against its registry.
        assert one({"action": "get_evidence", "evidence_id": evidence})["args"]["evidence_ids"] == [evidence]
    literal = 'AA,}" \\ [{'
    note = one({"action": "get_evidence", "evidence_id": "E", "note": literal})
    assert note["metadata"]["note"] == literal


def test_truncation_wrapper_conflicts_embedded_instructions():
    ok = '{"action":"submit"}'
    for raw in (ok[:-1], "[" + ok + ",", ok + " trailing {", ok + " " + ok,
                '{"action":"get_evidence","note":' + ok + '}',
                '{"actions":[' + ok + '],"batch_scope":"ignore donor B"}',
                '{"raw_tool_result":' + ok + '}', '{"action":"submit",,}',
                '{"action":"submit","metadata":{,}}', '[,' + ok + ']'):
        assert "error" in one(raw)
    for origin in ("tool", "literature", "user", "parser"):
        assert "error" in one(ok, origin=origin)
    assert "error" in one(ok, complete=False)
    empty_candidate_selection = one({"action": "submit", "candidate_ids": []})
    assert empty_candidate_selection["args"] == {"candidate_ids": []}
    assert "error" in one({"action": "submit", "candidate_ids": None})


def test_malformed_item_is_local_in_a_decodable_batch():
    raw = ('[{"action":"submit","action":"list_candidates"},'
           '{"action":"coexpress","candidate":"X","genes":"A","threshold":NaN},'
           '{"action":"list_candidates","offset":"1e1"}]')
    result = parse_message(raw)
    assert [r.get("status") for r in result] == ["needs_clarification", "needs_clarification", None]
    assert result[-1]["args"] == {"offset": 10}
    assert "error" in one('[{"action":"submit"},' + '{"action":"submit"},' * MAX_ACTIONS + ']')


def test_transactional_revision_group_normalization_preserves_sibling_roles():
    base = {"操作": "修订候选", "parents": ["H2", "H1"], "理由": "竞争候选的对照差异",
            "证据编号": "E1", "群体": [{"members": ["C2", "C1"], "status": "provisional"},
                                    {"candidate_ids": "H2", "label": "unknown"}]}
    good, read = parse_message(json.dumps([base, {"action": "get_evidence", "candidate": "H1"}]))
    assert good["args"]["groups"] == [{"members": ["C1", "C2"], "status": "provisional"},
                                     {"candidate_ids": ["H2"], "proposed_label": "unknown"}]
    assert read["args"]["candidate_id"] == "H1"
    for group in ({"members": None}, {"exclude": ["C3"]}, {"members": "C", "prior": 0.99}):
        invalid = {**base, "群体": [base["群体"][0], group]}
        failed, independent = parse_message(json.dumps([invalid, {"action": "submit"}]))
        assert "error" in failed and "groups" not in failed
        assert independent["action"] == "submit"
    assert "error" in one({**base, "证据编号": []})


def test_revision_source_inference_requires_one_group_with_exact_refs():
    group = {"candidate": "H1", "理由": "seed supports membership", "证据编号": "E1", "备注": "small sample"}
    revision = {"action": "propose_revision", "groups": [group]}
    parsed = one(revision)
    assert parsed["args"]["groups"][0] == {
        "candidate_ids": ["H1"], "reason": group["理由"], "evidence_ids": ["E1"]}
    assert parsed["metadata"]["group_notes"] == [{"index": 0, "note": "small sample"}]
    for changes in ({"parents": None}, {"groups": [group, group]},
                    {"groups": [{"members": ["cell"], "reason": "observed", "evidence_ids": ["E1"]}]},
                    {"groups": [dict(group, threshold=0.5)]}):
        assert one({**revision, **changes})["status"] == "needs_clarification"


def test_gene_resolution_uses_versioned_stable_ids_and_exact_species():
    rows = [{"gene_id": "g1", "symbol": "TCR", "species": "human"},
            {"gene_id": "g2", "symbol": "TCR", "species": "mouse"},
            {"gene_id": "g3", "symbol": "TCRB", "species": "human"}]
    assert resolve_genes(["g1", "TCRB"], rows) == ["g1", "g3"]
    assert resolve_genes("TCR", rows, species="human") == ["g1"]
    with pytest.raises(ProtocolError, match="ambiguous gene"):
        resolve_genes("TCR", rows)
    for token in ("TCR8", "tcrb", "g1 "):
        with pytest.raises(ProtocolError, match="unknown gene"):
            resolve_genes(token, rows)


def test_nested_real_model_action_and_parameter_conflicts():
    nested = '{"action":{"type":"expand_candidate","candidate_id":"H017"}}\nScientific rationale: inspect related cells.'
    flat = '{"action":"expand_candidate","candidate_id":"H017"}'
    assert parse_message(nested) == parse_message(flat)
    for raw in [
        '{"action":{"type":"expand_candidate","candidate_id":"H017"},"candidate_id":"H071"}',
        '{"action":{"name":"expand_candidate","candidate_id":"H017"},"operation":"submit"}',
        '{"action":{"type":"expand_candidate","candidate_id":"H017","exclude_low_quality":true}}',
        '{"action":{"type":"submit","candidate_ids":["H017"]},"scope":"all"}',
        '{"action":{},"operation":"submit"}',
    ]:
        assert parse_message(raw)[0]["status"] == "needs_clarification"
    assert parse_message(nested, origin="tool")[0]["status"] == "needs_clarification"
    assert parse_message(nested, complete=False)[0]["status"] == "needs_clarification"
