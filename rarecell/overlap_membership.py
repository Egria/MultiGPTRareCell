"""Export a cell-level overlapping submission from both frozen ballot ledgers."""
import json
from collections import defaultdict

from .lineage_review import relationships
from .single_choice import DIRECT_MODES, _reachable


def _undefeated(rivals, edges):
    incoming, outgoing = defaultdict(set), defaultdict(set)
    for winner, loser in edges:
        if winner in rivals and loser in rivals:
            incoming[loser].add(winner)
            outgoing[winner].add(loser)
    # An external win against any member defeats the entire cycle.
    return {g for g in rivals if not incoming[g]
            or _reachable(g, incoming) <= _reachable(g, outgoing)}


def retained_groups(primary, ballots, labels, ancestors):
    """Keep the primary plus plausible, specific rivals from the frozen ballots."""
    rivals = {g for a, b, _, _, _ in ballots for g in (a, b)}
    broader = {name for g in rivals for name in ancestors.get(labels[g], ())}
    eligible = {g for g in rivals if labels[g] not in broader}
    biology, unanimous_geometry = [], []
    for left, right, winner, mode, left_fraction in ballots:
        edge = (winner, right if winner == left else left)
        if mode in DIRECT_MODES:
            biology.append(edge)
        elif mode == "geometry" and left_fraction in (0, 1):
            unanimous_geometry.append(edge)
    viable = _undefeated(eligible, biology)
    return {primary} | _undefeated(viable, unanimous_geometry)


def write_memberships(directory, result, partition):
    from .finalize import tsv

    labels = {g["id"]: g["label"] for g in result["groups"]}
    reliable = {g["id"]: g["label"] if g["reliable"] else None for g in result["groups"]}
    ancestors = relationships(result["panels"])
    owners = {c: g["id"] for g in result["groups"] for c in g["members"]}
    # Retain the partition group IDs for joins. A group with no primary owners
    # still has its immutable vote-group ID and may contribute secondary members.
    group_ids = {g: "vote/" + g for g in labels}
    group_ids.update({owners[r["cell_id"]]: r["group_id"] for r in partition if r["cell_id"] in owners})
    by_cell = defaultdict(list)
    for ledger in (result, result.get("unclaimed_conflicts", {})):
        for pair in ledger.get("pair_votes", []):
            left, right = pair["pair"].split("/")
            for cell, winner, mode, fraction in pair["outcomes"]:
                by_cell[cell].append((left, right, winner, mode, fraction))
    if by_cell.keys() - owners.keys():
        raise ValueError("Reviewed conflict has no primary owner")
    rows, ambiguous = [], []
    for row in partition:
        cell = row["cell_id"]
        primary = owners.get(cell)
        if primary is None:
            selections = [None]
        else:
            kept = retained_groups(primary, by_cell.get(cell, ()), reliable, ancestors)
            selections = [primary, *sorted(kept - {primary})]
        memberships = []
        for g in selections:
            entry = dict(cell_id=cell, group_id=group_ids[g] if g else row["group_id"],
                         proposed_label=labels[g] if g else row["proposed_label"],
                         vote_group_id=g or "", is_primary=int(g == primary), status=row["status"])
            rows.append(entry)
            memberships.append(entry)
        if len(memberships) > 1:
            ambiguous.append(dict(cell_id=cell, primary_group_id=row["group_id"], n_groups=len(memberships),
                                  group_ids=json.dumps([e["group_id"] for e in memberships]),
                                  labels=json.dumps([e["proposed_label"] for e in memberships], ensure_ascii=False)))
    tsv(directory / "memberships.tsv",
        ["cell_id", "group_id", "proposed_label", "vote_group_id", "is_primary", "status"], rows)
    tsv(directory / "ambiguous_cells.tsv",
        ["cell_id", "primary_group_id", "n_groups", "group_ids", "labels"], ambiguous)
    return dict(n_ambiguous_cells=len(ambiguous), n_memberships=len(rows),
                n_membership_groups=len({r["group_id"] for r in rows}))
