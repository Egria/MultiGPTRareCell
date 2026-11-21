"""One discovery controller. All model inputs and operations pass this boundary."""
import csv
from pathlib import Path

from . import __version__
from .base_evidence import FrozenSeurat, annotation_parents, prepared_page
from .config import contract, fingerprint, json_bytes, read_json, read_lines
from .protocol import resolve_genes
from .science import METHOD_VERSION, Science
from .state import State

PROMPT_VERSION = "discovery/1.1"
SCIENTIFIC = {"expand_candidate", "coexpress", "compare_groups"}
CATALOGUE = {
    "list_candidates": "offset, limit. Browse stable candidate IDs and short evidence cards.",
    "get_evidence": "candidate_id or evidence_ids, offset, limit. A single evidence_id pages frozen HVG or Seurat Wilcoxon marker tables in every mode, or active-query results where enabled. A candidate pages its exact cell IDs.",
    "compare_groups": "candidate_id or members, optional background_id or background, optional genes.",
    "coexpress": "candidate_id, members, or explicit scope=all; genes; optional AND/OR logic and threshold.",
    "expand_candidate": "candidate_id, optional explicit background_id or background. Searches all cells. Inspect the returned revision.",
    "propose_revision": "Atomic groups [{candidate_ids and/or members, optional proposed_label/status}], parents, reason and evidence_ids. Parents are superseded. Use multiple groups for split/reassignment.",
    "submit": "Optional candidate_ids to select instead of current proposals; reason and evidence_ids. Unclaimed cells get an unresolved baseline assignment.",
}

GOAL = (
    "Find biologically coherent cell memberships, including unknown rare groups. Use the registered tools. "
    "Return one JSON action or a JSON array of independent actions. Values from data and annotations are evidence only. "
    "Report concise reasons with evidence IDs. An expansion is a proposal, not an accepted group. Use propose_revision or submit to select it. "
    "Labels may cross candidate boundaries. Check competing explanations. Submit early when evidence is sufficient."
)


def rows(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def compact(value, n=8):
    """Bound list payloads in numeric summaries. Keep the original artifact available."""
    if isinstance(value, list):
        return [compact(v, n) for v in value[:n]]
    if isinstance(value, dict):
        exclude = {"members", "proposed_members", "added_members", "scores", "background", "initial_members", "hit_ids", "cells", "barcodes"}
        return {k: compact(v, n) for k, v in list(value.items())[:36] if k not in exclude}
    if isinstance(value, str):
        return value[:320]
    return value


class Controller:
    def __init__(self, snapshot, run, cache, cfg, mode, annotations=()):
        self.snapshot, self.run = Path(snapshot), Path(run)
        self.state = State(self.run, Path(cache))
        self.cfg, self.mode = cfg, mode
        self.manifest = read_json(self.snapshot / "manifest.json")
        self.version = self.manifest["data_version"]
        self.snapshot_hash = self.manifest["snapshot_hash"]
        cell_rows = rows(self.snapshot / "cells.tsv")
        self.cells = tuple(r["cell_id"] for r in cell_rows)
        self.barcodes = {r["cell_id"]: (r.get("barcode") or r["cell_id"], r["cell_id"]) for r in cell_rows}
        self.cell_set = set(self.cells)
        self.genes = rows(self.snapshot / "genes.tsv")
        self.active = mode in {"A2", "A3"}
        self.labeled = mode in {"A1", "A3"}
        self.budget = cfg["budget"]
        self.annotations = list(annotations) if self.labeled else []
        self.science = None
        self.candidates = {}
        self.by_set = {}
        self._load_candidates()
        self.labels = {}
        if self.labeled:
            for row in self.state.all("Annotation"):
                self.labels.setdefault(row["candidate_id"], {})[row["evidence_id"]] = row

    def _index(self, candidate):
        self.candidates[candidate["id"]] = candidate
        self.by_set.setdefault(candidate["set_hash"], set()).add(candidate["id"])

    def _load_candidates(self):
        self.candidates.clear()
        self.by_set.clear()
        for item in self.state.all("Candidate"):
            self._index(item)

    def evidence(self, obj, kind, parents, family, members_hash, summary, refs, history=None,
                 purpose=None, method=None, prompt=None):
        return self.state.evidence({
            "object": obj, "kind": kind, "parent_ids": sorted(set(parents)), "source_family": family,
            "cell_set_hash": members_hash, "data_version": self.version, "method_version": method or ("LLMRareCell/" + __version__),
            "config_version": self.snapshot_hash, "prompt_version": prompt or PROMPT_VERSION,
            "selection_history": history or [], "artifact_refs": refs, "summary": summary,
            "purpose": purpose or ("exploration" if kind == "observation" else "related_check"),
        })

    def initialize(self):
        try:
            self._initialize()
        except BaseException:
            self._load_candidates()
            self.labels.clear()
            raise

    def _initialize(self):
        if self.state.get_meta("initialized", False):
            return
        with self.state.transaction():
            self.state.put("ProjectContract", contract()["version"], contract())
            self.state.put("Dataset", self.version, {"data_version": self.version, "snapshot_hash": self.snapshot_hash})
            for view in self.manifest["views"]:
                self.state.put("View", view, {"id": view, "snapshot_hash": self.snapshot_hash})
            prepared = FrozenSeurat(self)
            for raw in read_lines(self.snapshot / "initial.jsonl"):
                cid = raw["id"]
                members = sorted(set(raw["members"]))
                if not members or not set(members) <= self.cell_set or cid in self.candidates:
                    raise ValueError("Invalid initial candidate " + cid)
                sh = fingerprint(members)
                marker_context = {k: raw[k] for k in ("marker_kind", "marker_status", "marker_background", "structure") if k in raw}
                ref = self.state.artifact({"candidate_id": cid, "sources": raw["sources"],
                                           "markers": raw.get("markers", []), **marker_context})
                basis = [prepared.hvg] if prepared.hvg and any(s.get("view") == "V_HVG" for s in raw["sources"]) else []
                observed = self.evidence(cid, "observation", basis, "input_expression", sh,
                                         {"size": len(members), "sources": compact(raw["sources"], 4),
                                          "markers": compact(raw.get("markers", []), 5), **compact(marker_context)}, [ref])
                c = dict(raw, id=cid, members=members, set_hash=sh, evidence_ids=[observed],
                         parents=[], roots=[cid], initial=True, proposed_label="unknown", status="provisional")
                tables = prepared.markers(c, raw["sources"], observed)
                if tables:
                    c["seurat_markers"] = tables
                    c["evidence_ids"].extend(dict.fromkeys(t["evidence_id"] for t in tables))
                self.save_candidate(c)
                self.state.put("Task", "review/" + cid, {"id": "review/" + cid, "candidate_id": cid,
                               "cell_set_hash": sh, "initial_status": "pending", "mode": self.mode})
            for seed in read_lines(self.snapshot / "seeds.jsonl"):
                self.state.put("Seed", seed["id"], seed)
            for ann in self.annotations:
                c = self.candidate(ann["candidate_id"])
                ref = self.state.artifact(ann)
                if ann.get("cell_set_hash", c["set_hash"]) != c["set_hash"]:
                    raise ValueError("Annotation membership mismatch for " + c["id"])
                aid = "A" + fingerprint(ann)
                projection = {k: ann[k] for k in ("label", "ontology_id", "resolution", "partition_id",
                                                  "marker_hash", "model", "cache_key") if k in ann}
                eid = self.evidence(aid, "interpretation", annotation_parents(ann, c), "input_expression",
                                    c["set_hash"], projection, [ref], history=ann.get("selection_history"),
                                    method=ann.get("author_commit"), prompt=ann.get("prompt_version"))
                row = dict(ann, evidence_id=eid, id=aid)
                self.state.put("Annotation", aid, row)
                self.labels.setdefault(c["id"], {})[eid] = row
            self.state.put("Task", "discovery", {"id": "discovery", "goal": GOAL, "mode": self.mode})
            self.state.set_meta("initialized", True)
            self.state.event("initialized", {"mode": self.mode, "candidates": len(self.candidates)})

    def save_candidate(self, c):
        self.state.put("CellSet", c["set_hash"], {"cell_set_hash": c["set_hash"], "members": c["members"]})
        self.state.put("Candidate", c["id"], c)
        self._index(c)

    def candidate(self, cid):
        try:
            return self.candidates[cid]
        except (KeyError, TypeError):
            raise ValueError("Unknown candidate ID %s" % cid) from None

    def current_label(self, cid):
        items = self.labels.get(cid, {})
        if not items:
            return {}
        labels = sorted({item["label"] for item in items.values()})
        sources = []
        for eid, record in sorted(items.items()):
            summary = {k: record[k] for k in ("label", "partition_id", "resolution", "ontology_id") if k in record}
            sources.append(dict(summary, evidence_id=eid))
        return {"label": labels[0] if len(labels) == 1 else "unknown", "labels": labels,
                "evidence_ids": sorted(items), "sources": sources}

    def card(self, c):
        label = self.current_label(c["id"])
        card = {"id": c["id"], "n": len(c["members"]), "sources": compact(c["sources"], 2),
                "evidence_ids": c["evidence_ids"], "markers": compact(c.get("markers", []), 3),
                "parents": c["parents"], "status": c["status"]}
        card["selected"] = c["id"] in self.state.get_meta("selected", [])
        if c.get("structure"):
            card["structure"] = {view: {k: v for k, v in stats.items()
                                        if k in {"status", "within_fraction", "graph_coverage", "n_zero_degree"}}
                                 for view, stats in c["structure"].items()}
        if "marker_kind" in c:
            card["marker_kind"] = c["marker_kind"]
        if c.get("seurat_markers"):
            card["seurat_wilcoxon"] = [{k: row[k] for k in ("evidence_id","resolution","n_genes","first_markers")}
                                       for row in c["seurat_markers"][:2]]
        if label:
            card.update(tentative_label=label["label"], tentative_alternatives=label["labels"],
                        annotation_sources=label["sources"])
        if not c["initial"]:
            card["proposed_label"] = c["proposed_label"]
        return card

    def priority(self, c):
        first_barcode = min(self.barcodes[cell] for cell in c["members"])
        return (len(c["members"]), first_barcode, c["id"])

    def ordered(self):
        return sorted(self.candidates.values(), key=self.priority)

    def context(self, turn):
        s = self.state
        actions = {k: v for k, v in CATALOGUE.items() if self.active or k not in SCIENTIFIC}
        recent = s.get_meta("receipts", [])
        queries = s.all("CanonicalQuery")
        usage = s.usage()
        ctx = {"contract": contract(), "goal": GOAL, "mode": self.mode, "actions": actions,
               "bindings": {"run_id": s.get_meta("run_id"), "task_id": "discovery",
                            "data_version": self.version, "contract_version": contract()["version"]},
               "species": self.manifest.get("species"), "tissue": self.manifest.get("tissue"),
               "views": self.manifest["views"], "candidate_count": len(self.candidates),
               "selected_count": len(s.get_meta("selected", [])),
               "remaining": {"turns_including_this": max(0, self.budget["max_agent_turns"] - turn + 1),
                             "new_queries": self.budget["max_new_queries"] - usage["new_queries"] if self.active else 0,
                             "lineage_max_depth": self.budget["max_expansion_depth"]},
               "current_candidates": [], "recent_receipts": [], "prior_queries": [],
               "selected_ids": s.get_meta("selected", [])[:25]}
        ctx["prompt_version"] = PROMPT_VERSION
        ctx["turn_instruction"] = (
            "Inspect the candidate cards and the most recent results. Choose an action that resolves a specific membership question. "
            "Previously viewed pages are available by ID. Repeating a request yields no new evidence. "
            "Use submit when you can decide or no further query is justified. Write the chosen action and its named arguments."
        )
        if turn == self.budget["max_agent_turns"]:
            ctx["turn_instruction"] = "This is the final turn. Submit your best supported assignment with current evidence."
        if not self.active:
            ctx["mode_constraint"] = "Choose or combine complete precomputed candidate sets. No new scientific queries."
        # Reserve the entire contract, catalogue and budgets in every turn.
        limit = self.budget["context_bytes"]
        if len(json_bytes(ctx)) > limit:
            raise ValueError("context_bytes cannot hold complete discovery contract and required fields")
        def append_if_fits(key, obj):
            ctx[key].append(obj)
            if len(json_bytes(ctx)) > limit:
                ctx[key].pop()
                return False
            return True
        hvg = s.get_meta("global_hvg")
        if hvg:
            ctx["global_feature_selection"] = []
            if not append_if_fits("global_feature_selection", hvg):
                ctx.pop("global_feature_selection")
        seen_receipts = set()
        shown = set()
        for receipt in reversed(recent[-12:]):
            key = fingerprint(receipt)
            if key in seen_receipts:
                continue
            seen_receipts.add(key)
            if not append_if_fits("recent_receipts", receipt):
                continue
            shown.update(c["id"] for c in receipt.get("candidates", []))
            if "candidate" in receipt:
                shown.add(receipt["candidate"]["id"])
        ctx["recent_receipts"].reverse()
        for c in self.ordered():
            if c["id"] in shown:
                continue
            if len(ctx["current_candidates"]) >= 5 or not append_if_fits("current_candidates", self.card(c)):
                break
        for q in queries[-25:]:
            entry = {"key": q.get("canonical_query_key", q.get("key")), "action": q.get("canonical", {}).get("action")}
            if not append_if_fits("prior_queries", entry):
                break
        return ctx

    def resolve_members(self, args, id_key="candidate_id", list_key="members"):
        if id_key in args:
            return self.candidate(args[id_key])["members"]
        members = args.get(list_key)
        if not members or not set(members) <= self.cell_set:
            raise ValueError("Missing scope or unknown exact cell ID")
        return sorted(set(members))

    def scientific(self, action, args):
        if not self.active:
            raise ValueError("New scientific queries are disabled in " + self.mode)
        all_scope = args.get("scope") == "all"
        members = sorted(self.cells) if all_scope else self.resolve_members(args)
        bound = {"members": members}
        if "background_id" in args or "background" in args:
            bound["background"] = self.resolve_members(args, "background_id", "background")
            if set(members).intersection(bound["background"]):
                raise ValueError("Reference background intersects foreground")
        elif action == "compare_groups":
            complement = sorted(self.cell_set.difference(members))
            if complement:
                bound["background"] = complement
        if "genes" in args:
            bound["genes"] = sorted(set(resolve_genes(args["genes"], self.genes, species=self.manifest.get("species"))))
        if action == "coexpress":
            bound.update(logic=args.get("logic", "AND"), threshold=float(args.get("threshold", 0)))
        canonical = {"action": action, "args": bound, "data_version": self.version,
                     "transform_version": self.snapshot_hash, "method_version": [__version__, METHOD_VERSION],
                     "parameters": self.cfg["search"], "layer": "log1p"}
        signature = fingerprint(canonical)
        c = self.candidate(args["candidate_id"]) if "candidate_id" in args else None
        depths = self.state.get_meta("depths", {})
        seen = self.state.get_meta("expansions", {})
        if action == "expand_candidate":
            if c is None:
                raise ValueError("Expansion needs an exact lineage anchor")
            if signature not in seen and max(depths.get(r, 0) for r in self.lineage_roots(c)) >= self.budget["max_expansion_depth"]:
                raise ValueError("Expansion depth exhausted for this lineage")
        def compute():
            if self.science is None:
                self.science = Science(self.snapshot, self.cfg)
            return self.science.execute(action, bound)
        result, key = self.state.query(canonical, self.budget["max_new_queries"], compute)
        ref = self.state.result_ref(key)
        try:
            with self.state.transaction():
                return self._register_result(action, args, bound, c, signature, result, key, ref)
        except BaseException:
            # A failed DB commit or a failed compound edit must also discard the
            # uncommitted in-memory search/revision index for this live controller.
            self._load_candidates()
            raise

    def _selection_parents(self, args, bound):
        parents = set()
        for cid, member_key in (("candidate_id", "members"), ("background_id", "background")):
            if cid in args:
                matches = [args[cid]]
            else:
                matches = self.by_set.get(fingerprint(bound.get(member_key, [])), ())
            for match in matches:
                parents.update(self.candidate(match)["evidence_ids"])
        return parents

    def _register_result(self, action, args, bound, c, signature, result, key, ref):
        parents = self._selection_parents(args, bound)
        sh = fingerprint(bound["members"])
        summary = compact(result.get("summary", {"status": result.get("status")}))
        hist = [{"query_key": key, "selection_dependent": True, "selection_history_ref": ref}]
        try:
            eid = self.state.get("QueryObservation", key)["evidence_id"]
        except KeyError:
            eid = self.evidence(key, "observation", parents, "input_expression", sh, summary, [ref],
                                hist, purpose="selection_dependent", method=METHOD_VERSION)
            self.state.put("QueryObservation", key, {"query_key": key, "evidence_id": eid})
        # One numerical observation may be re-read through a newly selected
        # subset. Record that context without manufacturing a second observation.
        selection_ids = sorted(parents.difference({eid}))
        use = {"canonical_query_key": key, "evidence_id": eid, "selection_evidence_ids": selection_ids}
        use_id = "U" + fingerprint(use)
        self.state.put("QueryUse", use_id, dict(use, id=use_id))
        response = {"action": action, "query_key": key, "evidence_id": eid, "query_use_id": use_id,
                    "status": result.get("status", "ok"), "summary": summary, "artifact_ref": ref}
        if action != "expand_candidate":
            return response

        roots = self.lineage_roots(c)
        depths = self.state.get_meta("depths", {})
        seen = self.state.get_meta("expansions", {})
        depth = max(depths.get(r, 0) for r in roots)
        if signature in seen:
            used_depth = max(depth, seen[signature])
        else:
            if depth >= self.budget["max_expansion_depth"]:
                raise ValueError("Expansion depth exhausted for this lineage")
            used_depth = depth + 1
            seen[signature] = used_depth
        for root in roots:
            depths[root] = used_depth
        self.state.set_meta("depths", depths)
        self.state.set_meta("expansions", seen)
        if result["added_members"]:
            revision = self.make_revision(result["proposed_members"], [c["id"]], "unknown", "provisional",
                                          [eid], "expansion", key)
            response.update(candidate_id=revision["id"], added_count=len(result["added_members"]))
        # Same numerical query can be addressed through distinct revisions.
        investigation_id = "I" + fingerprint([signature, c["id"]])
        self.state.put("Investigation", investigation_id,
                       {"id": investigation_id, "query_key": key, "candidate_id": c["id"],
                        "proposed_revision_id": response.get("candidate_id"), "evidence_id": eid})
        return response

    def lineage_roots(self, c):
        """Follow parents plus equal member sets, including equivalences registered later."""
        roots, visited = set(), set()
        stack = [c["id"]]
        while stack:
            cid = stack.pop()
            if cid in visited:
                continue
            visited.add(cid)
            current = self.candidate(cid)
            roots.update(current["roots"])
            stack.extend(current["parents"])
            stack.extend(self.by_set.get(current["set_hash"], ()))
        return roots

    def make_revision(self, members, parents, label, status, evidence_ids, kind, basis):
        if not members or not set(members) <= self.cell_set:
            raise ValueError("Revision contains empty or unknown membership")
        members = sorted(set(members))
        sh = fingerprint(members)
        parents = sorted(set(parents))
        roots = {r for p in parents for r in self.lineage_roots(self.candidate(p))}
        for same in self.by_set.get(sh, ()):
            roots.update(self.lineage_roots(self.candidate(same)))
        if kind == "revision":
            # Carry forward expansion provenance even if a split cites the
            # retrieval evidence alongside a different coarse parent.
            inputs = set(evidence_ids)
            ancestors = set()
            while inputs:
                eid = inputs.pop()
                if eid not in ancestors:
                    ancestors.add(eid)
                    inputs.update(self.state.get("Evidence", eid)["parent_ids"])
            for investigation in self.state.all("Investigation"):
                if investigation["evidence_id"] in ancestors:
                    roots.update(self.lineage_roots(self.candidate(investigation["candidate_id"])))
        if not roots:
            raise ValueError("Revision requires traceable parent candidates")
        descriptor = [sh, parents, label, status, sorted(set(evidence_ids)), kind, basis]
        cid = "R" + fingerprint(descriptor)[:24]
        if cid in self.candidates:
            return self.candidate(cid)
        evidence_ids = set(evidence_ids)
        if kind == "revision":
            # A model interpretation remains connected to the marker/structural
            # observations and to any tentative labels available on its inputs.
            for parent in parents:
                evidence_ids.update(self.candidate(parent)["evidence_ids"])
                evidence_ids.update(self.current_label(parent).get("evidence_ids", []))
            interpretation = self.evidence(cid, "interpretation", evidence_ids, "discovery_agent",
                                           sh, {"reason": basis, "proposed_label": label,
                                                "status": status, "size": len(members)}, [],
                                           history=[{"parent_candidates": parents}])
            evidence_ids.add(interpretation)
        record = {"id": cid, "members": members, "set_hash": sh, "parents": parents, "roots": sorted(roots),
                  "initial": False, "sources": [{"kind": kind}], "evidence_ids": sorted(set(evidence_ids)),
                  "proposed_label": label, "status": status, "markers": []}
        self.save_candidate(record)
        return record

    def checked_evidence(self, ids):
        for eid in ids:
            try:
                self.state.get("Evidence", eid)
            except KeyError:
                raise ValueError("Unknown evidence ID " + eid) from None
        return ids

    def decide(self, action, args, message_id):
        selected = set(self.state.get_meta("selected", []))
        evidence_ids = self.checked_evidence(args.get("evidence_ids", []))
        reason = args.get("reason", "Submit current provisional assignment")
        if action == "submit":
            selected = set(args.get("candidate_ids", selected))
            for cid in selected:
                evidence_ids = evidence_ids + self.candidate(cid)["evidence_ids"]
                evidence_ids = evidence_ids + self.current_label(cid).get("evidence_ids", [])
        else:
            parents = [self.candidate(cid)["id"] for cid in args["parents"]]
            if not parents or not evidence_ids:
                raise ValueError("Revisions require parent candidates and registered evidence")
            chosen = []
            occupied = set()
            for group in args["groups"]:
                refs = [self.candidate(cid) for cid in group.get("candidate_ids", [])]
                member_set = {m for c in refs for m in c["members"]}
                explicit = set(group.get("members", []))
                if explicit and not self.active and not self.by_set.get(fingerprint(sorted(explicit))):
                    raise ValueError("Passive mode can select or combine only precomputed member sets")
                member_set |= explicit
                if occupied.intersection(member_set):
                    raise ValueError("New groups in a single revision overlap")
                occupied.update(member_set)
                inherited = set(parents) | {c["id"] for c in refs}
                chosen.append(self.make_revision(member_set, inherited, group.get("proposed_label", "unknown"),
                                                 group.get("status", "provisional"), evidence_ids, "revision", reason)["id"])
            selected.difference_update(parents)
            selected.update(chosen)
            evidence_ids = list(evidence_ids) + [eid for cid in chosen for eid in self.candidate(cid)["evidence_ids"]]
        decision = {"action": action, "reason": reason, "evidence_ids": sorted(set(evidence_ids)),
                    "selected": sorted(selected), "message_id": message_id, "args_hash": fingerprint(args)}
        if action == "submit":
            decision["status"] = args.get("status", "final")
            if evidence_ids:
                submitted = sorted(selected)
                interpretation = self.evidence("submitted_groups", "interpretation", evidence_ids, "discovery_agent",
                                               fingerprint(sorted({m for cid in submitted
                                                                   for m in self.candidate(cid)["members"]})),
                                               {"reason": reason, "status": decision["status"],
                                                "group_count": len(selected)}, [],
                                               history=[{"selected_candidates": submitted}])
                decision["interpretation_id"] = interpretation
        did = "D" + fingerprint(decision)[:24]
        self.state.put("Decision", did, dict(decision, id=did))
        self.state.set_meta("selected", sorted(selected))
        if action == "submit":
            self.state.set_meta("submitted", True)
            self.state.set_meta("submission_status", decision["status"])
        return {"decision_id": did, "action": action, "selected_count": len(selected),
                "candidate_ids": sorted(selected)[:15]}

    def read_tool(self, action, args):
        offset = args.get("offset", 0)
        limit = min(args.get("limit", 8), 25)
        if action == "list_candidates":
            ordered = self.ordered()
            return {"total": len(ordered), "offset": offset, "candidates": [self.card(c) for c in ordered[offset:offset+limit]]}
        if "candidate_id" in args:
            c = self.candidate(args["candidate_id"])
            return {"candidate": self.card(c), "offset": offset, "total": len(c["members"]),
                    "members": c["members"][offset:offset+limit], "markers": compact(c.get("markers", [])[offset:offset+limit], limit)}
        ids = args.get("evidence_ids")
        if ids and len(ids) == 1:
            card = self.state.get("Evidence", self.checked_evidence(ids)[0])
            prepared = prepared_page(self.state, card, offset, limit)
            if prepared is not None:
                return prepared
            if self.active and card["kind"] == "observation":
                try:
                    fact = self.state.get("QueryObservation", card["object"])
                except KeyError:
                    fact = None
                if fact and fact["evidence_id"] == ids[0]:
                    key = fact["query_key"]
                    ref = self.state.result_ref(key)
                    if ref not in card["artifact_refs"]:
                        raise ValueError("Scientific observation/result reference mismatch")
                    data = self.state.read_artifact(ref)
                    history = data.get("selection_history", {})
                    sections = {}
                    for field in ("matched_members", "statistics", "scores"):
                        if field in data:
                            sections[field] = data[field]
                    if "gene_counts" in data:
                        sections["gene_counts"] = [{"gene_id": g, "n_detected": n}
                                                  for g, n in sorted(data["gene_counts"].items())]
                    for field in ("signature", "background"):
                        if field in history:
                            sections[field] = history[field]
                    return {"evidence_id": ids[0], "query_key": key, "offset": offset, "limit": limit,
                            "status": data["status"], "summary": compact(data.get("summary", {})),
                            "sections": {field: {"total": len(values), "rows": values[offset:offset+limit]}
                                         for field, values in sections.items()}}
        if ids is None:
            ids = [row["evidence_id"] for row in self.state.all("Evidence")]
        self.checked_evidence(ids)
        items = []
        for eid in ids[offset:offset+limit]:
            card = self.state.get("Evidence", eid)
            items.append({k: card[k] for k in ("evidence_id", "object", "kind", "summary", "parent_ids", "purpose", "source_family", "artifact_refs")})
        return {"total": len(ids), "offset": offset, "evidence": items}

    def execute(self, parsed, message_id):
        if "error" in parsed:
            raise ValueError(parsed["error"])
        action = parsed["action"]
        args = dict(parsed["args"])
        bindings = dict(run_id=self.state.get_meta("run_id"), task_id="discovery", data_version=self.version,
                        contract_version=contract()["version"], message_id=message_id)
        for key, expected in bindings.items():
            if key in args and args.pop(key) != expected:
                raise ValueError("Conflicting explicit " + key)
        if action in SCIENTIFIC:
            return self.scientific(action, args)
        if action in {"list_candidates", "get_evidence"}:
            return self.read_tool(action, args)
        if action in {"propose_revision", "submit"}:
            try:
                with self.state.transaction():
                    return self.decide(action, args, message_id)
            except BaseException:
                self._load_candidates()
                raise
        raise ValueError("Unregistered tool")
