from rarecell.seed_review import screen


def test_signed_shared_marker_screen_includes_isolated_seeds():
    def seed(cid, genes, views):
        return {"id": cid, "markers": [{"gene_id": g, "effect": v} for g, v in genes],
                "sources": [{"view": v} for v in views]}
    candidates = [seed("C", [("A", -1), ("B", 1)], ["Y"]),
                  seed("B", [("A", 1), ("B", -1)], ["Y", "Z"]),
                  seed("D", [], []), seed("A", [("A", 2), ("B", -3)], ["Z", "Y"])]
    ranked = screen(candidates)
    assert ranked == screen(reversed(candidates))
    assert [(r["candidate_id"], r["peer"]) for r in ranked] == [("A", "B"), ("B", "A"), ("C", None), ("D", None)]
    assert ranked[0]["shared_markers"] == ["+A", "-B"] and ranked[0]["matching_seeds"] == 2
