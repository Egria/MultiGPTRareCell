"""Triage singleton anchors using frozen RNA observations, before any labels are known."""
from collections import defaultdict
from itertools import combinations


def screen(candidates):
    """Compare each seed's top three signed markers with all other seeds in linear space."""
    index = defaultdict(list)
    subsets = {}
    cards = {}
    for c in sorted(candidates, key=lambda c: c["id"]):
        cid = c["id"]
        markers = {(r["gene_id"], "+" if r["effect"] > 0 else "-"): r.get("symbol") or r["gene_id"]
                   for r in c.get("markers", [])[:3] if r["effect"]}
        subsets[cid] = [keys for n in range(1, len(markers) + 1) for keys in combinations(sorted(markers), n)]
        for keys in subsets[cid]:
            index[keys].append(cid)
        cards[cid] = {"candidate_id": cid, "views": sorted({s["view"] for s in c["sources"] if s.get("view")}),
                      "markers": markers}
    for cid, card in cards.items():
        markers = card.pop("markers")
        shared = [keys for keys in subsets[cid] if len(index[keys]) > 1]
        # Prefer longer shared patterns, then rarer patterns and positive expression.
        keys = min(shared, key=lambda k: (-len(k), len(index[k]), -sum(sign == "+" for _, sign in k), k),
                   default=())
        hits = index.get(keys, [])
        peer = next((other for other in hits[:2] if other != cid), None)
        card.update(peer=peer, shared_markers=[sign + markers[(gene, sign)] for gene, sign in keys],
                    matching_seeds=len(hits))
    return sorted(cards.values(), key=lambda r: (-len(r["shared_markers"]), -len(r["views"]),
                                                  r["matching_seeds"], r["candidate_id"]))
