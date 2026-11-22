"""Register already-frozen Seurat output; no expression calculations or model calls."""
import csv
import math
from collections import defaultdict

from .config import fingerprint
from .singleton_markers import SOURCE as SINGLETON_CONTRAST


HVG = "seurat_vst_features"
DISPERSION = "seurat_dispersion_features"
DE = "seurat_wilcoxon_markers"
NUMERIC = ("avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")


def _numeric(value):
    if value in ("NA", "NaN", ""):
        return None
    n = float(value)
    if not math.isfinite(n):
        raise ValueError("Invalid nonfinite Seurat marker statistic")
    return n


def _table(path):
    with path.open(encoding="utf8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class FrozenSeurat:
    def __init__(self, ctrl):
        self.ctrl = ctrl
        self.hvg = None
        self.ranges = defaultdict(list)
        self.fields = []
        self.universe_size = 0
        self.marker_path = ctrl.snapshot / "seurat" / "markers.tsv"
        manifest = ctrl.manifest
        stages = manifest.get("stages", {})
        if stages.get("V_HVG", {}).get("status") != "ok":
            return
        files = manifest["files"]
        genes = _table(ctrl.snapshot / "seurat" / "hvg.tsv")
        informative = sorted(r["cell_id"] for r in _table(ctrl.snapshot / "cell_qc.tsv")
                             if r["zero_expression"] == "0")
        self.universe_size = len(informative)
        self.universe_hash = fingerprint(informative)
        selection = stages["V_HVG"].get("selection_method", "vst")
        statistic = DISPERSION if selection == "dispersion" else HVG
        descriptor = {"statistic": statistic, "n_features": len(genes), "selection_method": selection,
                      "n_informative": self.universe_size, "source_sha256": files["seurat/hvg.tsv"],
                      "population_qc_sha256": files["cell_qc.tsv"]}
        self.hvg = self.register(statistic, [], self.universe_hash, descriptor, genes,
                                 "Seurat/5.3.0 FindVariableFeatures", "exploration")
        ctrl.state.set_meta("global_hvg", {"evidence_id": self.hvg, "n_features": len(genes),
                                          "statistic": statistic, "first_features": genes[:3]})
        if stages.get("seurat_markers", {}).get("status") == "ok":
            self.index_file()

    def index_file(self):
        # R emits contiguous group runs; index byte spans without retaining every
        # partition's gene table in RAM. Repeated spans for one group are valid.
        with self.marker_path.open("rb") as f:
            self.fields = next(csv.reader([f.readline().decode()], delimiter="\t"))
            if self.fields[:2] != ["resolution", "group_id"] or not set(NUMERIC).issubset(self.fields):
                raise ValueError("Unexpected prepared Seurat marker columns")
            previous = None
            for line in iter(f.readline, b""):
                res, cid, _ = line.split(b"\t", 2)
                key = (float(res), cid.decode())
                end = f.tell()
                if previous != key:
                    self.ranges[key].append([end-len(line), end])
                else:
                    self.ranges[key][-1][1] = end
                previous = key

    def table_for(self, key):
        table = []
        with self.marker_path.open("rb") as f:
            for start, end in self.ranges[key]:
                f.seek(start)
                lines = f.read(end-start).decode().splitlines()
                table.extend(csv.DictReader(lines, delimiter="\t", fieldnames=self.fields))
        return table

    def register(self, obj, parents, sh, summary, table, method, purpose):
        ctrl = self.ctrl
        ref = ctrl.state.artifact({"statistic": summary["statistic"], "rows": table})
        eid = ctrl.evidence(obj, "observation", parents, "input_expression", sh, summary, [ref],
                            method=method, purpose=purpose,
                            history=[{"source_sha256": summary["source_sha256"],
                                      "selection": summary["statistic"]}])
        ctrl.state.put("PreparedTable", eid, {"evidence_id": eid, "ref": ref,
                                              "statistic": summary["statistic"]})
        return eid

    def markers(self, candidate, sources, observed):
        info = []
        count = len(candidate["members"])
        background_count = self.universe_size - count
        if not self.fields or count < 3 or background_count < 3:
            return info
        for source in sources:
            if source.get("kind") != "partition" or source.get("view") != "V_HVG":
                continue
            key = (float(source["resolution"]), candidate["id"])
            raw = self.table_for(key)
            # The OntoAnno cache hashes the exact per-source TSV row dictionaries.
            # Exclude partition identifiers from numerical rows to deduplicate a
            # repeated test on exactly the same cells across multiple resolutions.
            values = [dict(gene_id=r["gene_id"], symbol=r.get("symbol",""),
                           **{k: _numeric(r[k]) for k in NUMERIC}) for r in raw]
            summary = {"statistic": DE, "n_genes": len(values), "n_foreground": count,
                       "background_rule": "all informative input cells outside candidate",
                       "n_background": background_count, "population_hash": self.universe_hash,
                       "source_sha256": self.ctrl.manifest["files"]["seurat/markers.tsv"],
                       "only_pos": True, "min_pct": 0.01, "logfc_threshold": 0}
            eid = self.register(DE + "/" + candidate["id"], [observed], candidate["set_hash"],
                                summary, values, "Seurat/5.3.0 FindMarkers(wilcox)", "selection_dependent")
            info.append({"evidence_id": eid, "marker_hash": fingerprint(raw),
                         "partition_id": source["partition_id"], "resolution": key[0],
                         "n_genes": len(values), "first_markers": values[:2]})
        return info


def annotation_parents(ann, candidate):
    if "marker_hash" not in ann and "source" not in ann:
        # Direct callers can register a plain interpretation in unit-level
        # controller scenarios. The CLI accepts only verified OntoAnno rows.
        return candidate["evidence_ids"]
    if ann.get("marker_source") == SINGLETON_CONTRAST:
        if (len(candidate["members"]) != 1 or fingerprint(candidate["marker_background"]) != ann["background_hash"]
                or fingerprint(ann["markers"]) != ann["marker_hash"]):
            raise ValueError("Annotation has no matching frozen singleton contrast")
        return [candidate["evidence_ids"][0]]
    matches = {item["evidence_id"] for item in candidate.get("seurat_markers", [])
               if item["marker_hash"] == ann.get("marker_hash")
               and item["partition_id"] == ann.get("partition_id")
               and item["resolution"] == ann.get("resolution") and item["n_genes"] > 0}
    if len(matches) != 1:
        raise ValueError("Annotation has no corresponding frozen Seurat marker observation")
    return list(matches)


def prepared_page(state, card, offset, limit):
    try:
        binding = state.get("PreparedTable", card["evidence_id"])
    except KeyError:
        return None
    ref = binding["ref"]
    if card["artifact_refs"] != [ref]:
        raise ValueError("Prepared observation artifact mismatch")
    data = state.read_artifact(ref)
    name = binding["statistic"]
    if data["statistic"] != name or name not in (DE, HVG, DISPERSION):
        raise ValueError("Unknown prepared observation table")
    rows = data["rows"]
    return {"evidence_id": card["evidence_id"], "offset": offset, "limit": limit, "status": "ok",
            "summary": card["summary"], "purpose": card["purpose"], "parent_ids": card["parent_ids"],
            "sections": {name: {"total": len(rows), "rows": rows[offset:offset+limit]}}}
