"""Print ground-truth cells and their containing Seurat/OntoAnno groups after discovery."""
import argparse
import csv
import sys
from pathlib import Path

from rarecell.annotation import verify_annotation
from rarecell.config import read_json, read_lines


def memberships(snapshot, annotations, targets):
    by_group = {(float(r["resolution"]), r["candidate_id"]): r for r in annotations if r.get("resolution") is not None}
    singles = {r["candidate_id"]: r for r in annotations if r.get("marker_source") == "singleton_positive_contrast"}
    seeds = {r["cell_id"]: r["candidate_id"] for r in read_lines(snapshot/"seeds.jsonl")}
    with (snapshot/"baselines.tsv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["view"] != "V_HVG" or row["cell_id"] not in targets:
                continue
            annotation = by_group.get((float(row["resolution"]), row["group_id"]), {})
            direct = singles.get(seeds.get(row["cell_id"]), {})
            yield dict(row, seed_label=direct.get("label", ""),
                       seed_marker_input=",".join(direct.get("selection_history", {}).get("genes", [])),
                       seed_candidate_id=seeds.get(row["cell_id"], ""),
                       parent_label=annotation.get("label", ""), ontology_id=annotation.get("ontology_id", ""),
                       marker_input=",".join(annotation.get("selection_history", {}).get("genes", [])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", required=True, type=Path, help="dashboard task containing run/checkpoint.json")
    parser.add_argument("--gt", required=True, type=Path, help="evaluation fixture with cell_id and label")
    parser.add_argument("--label", required=True, help="exact target label in the evaluation fixture")
    args = parser.parse_args()
    checkpoint = read_json(args.task/"run/checkpoint.json")
    if not checkpoint.get("annotations"):
        parser.error("this run has no upstream annotations")
    snapshot = Path(checkpoint["snapshot"])
    ann = verify_annotation(snapshot, Path(checkpoint["annotations"]), allow_partial=True)
    with args.gt.open(newline="", encoding="utf-8-sig") as stream:
        targets = {r["cell_id"] for r in csv.DictReader(stream, delimiter="\t") if r["label"] == args.label}
    if not targets:
        parser.error("ground truth contains no cells with this label")
    rows = sorted(memberships(snapshot, ann, targets), key=lambda r: (r["cell_id"], float(r["resolution"])))
    if targets != {r["cell_id"] for r in rows}:
        parser.error("one or more ground-truth cells lack V_HVG membership in this snapshot")
    writer = csv.DictWriter(sys.stdout, ["cell_id", "seed_candidate_id", "resolution", "group_id",
                            "parent_label", "ontology_id", "marker_input", "seed_label", "seed_marker_input"], delimiter="\t",
                            lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
