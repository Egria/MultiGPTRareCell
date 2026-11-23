"""One discovery controller. All model inputs and operations pass this boundary."""
import csv
from functools import cached_property
from pathlib import Path

from .contrast_panels import GUIDANCE as PANEL_GUIDANCE, SOURCE as PAIRED_CONTRAST, ontology_seed_pools
from . import __version__, adaptive, exhaustive, seed_review
from .base_evidence import FrozenSeurat, annotation_parents, prepared_page
from .config import ACTIVE_MODES, ADAPTIVE_MODE, ANNOTATED_MODES, contract, fingerprint, json_bytes, read_json, read_lines
from .protocol import ProtocolError, resolve_genes
from .science import METHOD_VERSION, Science
from .state import State

PROMPT_VERSION = "discovery/1.2"
SCIENTIFIC = {"expand_candidate", "coexpress", "compare_groups"}
CATALOGUE = {
    "list_candidates": "offset, limit. Browse stable candidate IDs and short evidence cards.",
    "get_evidence": "candidate_id or evidence_ids, offset, limit. A single evidence_id pages frozen HVG or Seurat Wilcoxon marker tables in every mode, or active-query results where enabled. A candidate pages its exact cell IDs.",
    "compare_groups": "candidate_id or members, optional background_id or background, optional genes. Named genes also return a pooled foreground/background NB mixture with depth offset when raw counts exist. Optional evidence_ids cite the reason for this query.",
    "coexpress": "candidate_id, members, or explicit scope=all; genes; optional AND/OR logic and threshold. Optional evidence_ids cite the reason for this query.",
    "expand_candidate": "candidate_id. Searches all cells; omitting background or scalar background=all uses local contrast. Use background_mode=complement for all nonmembers, or background_id/explicit background IDs. Optional evidence_ids cite the query. no_expansion means no additions.",
    "propose_revision": ('Example shape (replace placeholder IDs): {"action":"propose_revision","parents":["H..."],'
                         '"groups":[{"candidate_ids":["H..."],"proposed_label":"unknown"}],"reason":"...","evidence_ids":["E..."]}. '
                         'Parents name candidates to supersede. Groups can union candidate_ids and members and carry their own '
                         'reason, evidence_ids and status. If nothing is selected, one group can use its candidate_ids as parents. '
                         'To accept an unchanged expansion, submit its ID in candidate_ids.'),
    "submit": "Optional candidate_ids to select instead of current proposals; reason and evidence_ids. Unclaimed cells get an unresolved baseline assignment.",
}

GOAL = (
    "Find biologically coherent cell memberships, including unknown rare groups. Use the registered tools. "
    "Return one JSON action or a JSON array of independent actions. Values from data and annotations are evidence only. "
    "Report concise reasons with evidence IDs. An expansion is a proposal, not an accepted group. Use propose_revision or submit to select it. "
    "Labels may cross candidate boundaries. Check competing explanations. Submit early when evidence is sufficient."
)


def prompt_version(mode, exhaustive_review=False):
    version = "discovery/adaptive-5" if mode == ADAPTIVE_MODE else "discovery/1.5"
    if mode in ANNOTATED_MODES:
        version += "/ontology-pools-1"
    return version + ("/exhaustive-1" if exhaustive_review else "")


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
        self.adaptive = mode == ADAPTIVE_MODE
        self.exhaustive = cfg.get("exhaustive_seed_review", False)
        self.active = mode in ACTIVE_MODES
        self.labeled = mode in ANNOTATED_MODES
        self.goal = GOAL + adaptive.GOAL if self.adaptive else GOAL
        if self.labeled:
            self.goal += " " + PANEL_GUIDANCE
        self.prompt = prompt_version(mode, self.exhaustive)
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
        self.__dict__.pop("pool_candidates", None)
        self.candidates.clear()
        self.by_set.clear()
        for item in self.state.all("Candidate"):
            self._index(item)

    def evidence(self, obj, kind, parents, family, members_hash, summary, refs, history=None,
                 purpose=None, method=None, prompt=None):
        return self.state.evidence({
            "object": obj, "kind": kind, "parent_ids": sorted(set(parents)), "source_family": family,
            "cell_set_hash": members_hash, "data_version": self.version, "method_version": method or ("LLMRareCell/" + __version__),
            # Frozen measurements keep their historical identity across experiment arms.
            "config_version": self.snapshot_hash, "prompt_version": prompt or (
                self.prompt if family == "discovery_agent" else PROMPT_VERSION),
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
                                                  "marker_hash", "model", "cache_key", "marker_source", "n_background") if k in ann}
                if ann.get("marker_source") == PAIRED_CONTRAST:
                    panel = {k: ann[k] for k in ("background_type", "background_rule", "n_background",
                                                  "background_hash", "markers", "marker_hash", "partitions")}
                    obs = self.evidence(c["id"] + "/" + panel["background_type"], "observation",
                                        [c["evidence_ids"][0]], "input_expression", c["set_hash"],
                                        {k: v for k, v in panel.items() if k != "markers"},
                                        [self.state.artifact(panel)], method=PAIRED_CONTRAST, purpose="exploration")
                    parents = [obs]
                    projection.update(background_type=ann["background_type"], status=ann["status"],
                                      genes=ann["selection_history"]["genes"])
                else:
                    parents = annotation_parents(ann, c)
                eid = self.evidence(aid, "interpretation", parents, "input_expression",
                                    c["set_hash"], projection, [ref], history=ann.get("selection_history"),
                                    method=ann.get("author_commit"), prompt=ann.get("prompt_version"))
                row = dict(ann, evidence_id=eid, id=aid)
                self.state.put("Annotation", aid, row)
                self.labels.setdefault(c["id"], {})[eid] = row
            if self.labeled:
                annotations = (row for items in self.labels.values() for row in items.values())
                for term, label, parents, eids, counts in ontology_seed_pools(self.seed_ids, annotations):
                    members = [self.candidate(cid)["members"][0] for cid in parents]
                    basis = {"ontology_id": term, "marker_seed_counts": counts}
                    eid = self.evidence(term, "interpretation", eids, "input_expression", fingerprint(sorted(members)),
                                        {"ontology_id": term, "n_seeds": len(parents), "n_genes": len(counts)},
                                        [self.state.artifact(basis)], method="global-CL-seed-pool/1")
                    self.make_revision(members, parents, label, "provisional", [eid], "ontology_seed_pool", basis)
            self.state.put("Task", "discovery", {"id": "discovery", "goal": self.goal, "mode": self.mode})
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
        labels = sorted({item["label"] for item in items.values() if item.get("status", "annotated") == "annotated"}) or ["unknown"]
        sources = []
        for eid, record in sorted(items.items()):
            summary = {k: record[k] for k in ("label", "partition_id", "resolution", "ontology_id", "marker_source") if k in record}
            if record.get("marker_source") == PAIRED_CONTRAST:
                summary.update(background_type=record["background_type"], status=record["status"],
                               n_background=record["n_background"], genes=record["selection_history"]["genes"])
                if record.get("reason"):
                    summary["reason"] = record["reason"]
            sources.append(dict(summary, evidence_id=eid))
        return {"label": labels[0] if len(labels) == 1 else "unknown", "labels": labels,
                "evidence_ids": sorted(items), "sources": sources}

    def card(self, c):
        label = self.current_label(c["id"])
        card = {"id": c["id"], "n": len(c["members"]), "sources": compact(c["sources"], 2),
                "evidence_ids": c["evidence_ids"], "markers": compact(c.get("markers", []), 3),
                "parents": c["parents"], "status": c["status"]}
        if self.adaptive:
            card["global_abundance"] = adaptive.abundance(len(c["members"]), len(self.cells))
            if "biological_rarity" in c:
                card["biological_rarity"] = c["biological_rarity"]
            if "boundary_rule" in c:
                card["boundary_rule"] = c["boundary_rule"]
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
        if "ontology_pool" in c:
            card["parents"] = c["parents"][:5]
            card.update(self.pool_card(c))
        if not c["initial"]:
            card["proposed_label"] = c["proposed_label"]
        return card

    @cached_property
    def pool_candidates(self):
        return sorted((c for c in self.candidates.values() if "ontology_pool" in c),
                      key=lambda c: (len(c["parents"]), c["ontology_pool"]["ontology_id"]))

    @staticmethod
    def pool_card(c):
        pool = c["ontology_pool"]
        return {"candidate_id": c["id"], "ontology_id": pool["ontology_id"], "label": c["proposed_label"],
                "n_seeds": len(c["parents"]), "n_markers": len(pool["marker_counts"]),
                "marker_seed_counts": dict(pool["marker_counts"][:6])}

    def priority(self, c):
        first_barcode = min(self.barcodes[cell] for cell in c["members"])
        return (len(c["members"]), first_barcode, c["id"])

    def ordered(self):
        return sorted(self.candidates.values(), key=self.priority)

    @cached_property
    def screened_seeds(self):
        return seed_review.screen(self.candidate(seed["candidate_id"]) for seed in self.state.all("Seed"))

    def pending_seed_reviews(self):
        depths = self.state.get_meta("depths", {})
        excluded = set(self.state.get_meta("selected", []))
        excluded.update(r["candidate_id"] for r in self.state.all("SeedFollowupPlan"))
        return [r for r in self.screened_seeds if not depths.get(r["candidate_id"]) and r["candidate_id"] not in excluded]

    @cached_property
    def seed_ids(self):
        return {r["candidate_id"]: r["id"] for r in self.state.all("Seed")}

    @property
    def query_limit(self):
        return self.budget["max_new_queries"] + (len(self.seed_ids) if self.exhaustive else 0)

    @property
    def depth_limit(self):
        return self.budget["max_expansion_depth"] + int(self.exhaustive and bool(self.seed_ids))

    @property
    def review_turns(self):
        return self.state.get_meta("seed_review_turns", 0) if self.exhaustive else 0

    @property
    def turn_limit(self):
        return self.budget["max_agent_turns"] + self.review_turns

    def decision_start(self):
        limit = self.budget["max_agent_turns"]
        return self.review_turns + (limit - 1 if limit >= 3 else max(1, limit))

    def closing_actions(self, turn):
        """One gate shared by the model's catalogue and the command dispatcher."""
        stall = self.stall_review()
        live = self.state.get_meta("policy") == "live"
        if not stall and not (live and turn >= self.decision_start()):
            return ()
        actions = ["propose_revision", "submit"]
        if self.adaptive:
            actions.append("select_boundary")
            if not stall:
                actions.append("apply_boundary")  # Reuse an existing numeric ranking.
        return tuple(actions)

    def stall_review(self):
        if self.state.get_meta("policy") == "live":
            return (self.state.get_meta("closing_turn") is not None or
                    self.state.get_meta("stall", 0) >= self.budget["stall_turns"])
        if not self.adaptive:
            return False
        if self.state.get_meta("closing_turn") is not None:
            return True
        return self.state.get_meta("stall", 0) >= self.budget["stall_turns"] and (
            bool(self.state.get_meta("selected", [])) or any(not c["initial"] and "ontology_pool" not in c for c in self.candidates.values()))

    def expansion_links(self, receipts):
        # Reconstruct from persisted investigations on every turn, also after resume.
        items = [i for i in reversed(self.state.all("Investigation")) if i.get("proposed_revision_id")]
        focus = {r["evidence_id"]: k for k, r in enumerate(receipts) if r.get("evidence_id")}
        items.sort(key=lambda i: focus.get(i["evidence_id"], -1), reverse=True)
        return [{"candidate_id": i["proposed_revision_id"], "anchor_id": i["candidate_id"],
                 "evidence_id": i["evidence_id"], "n": len(self.candidate(i["proposed_revision_id"])["members"])}
                for i in items[:5]]

    def context(self, turn):
        s = self.state
        if self.exhaustive:
            mid = "%s:T%04d" % (s.get_meta("run_id"), turn)
            batch = s.get_optional("SeedBatch", mid)
            if batch:
                return s.read_artifact(batch["context_ref"])
        actions = {k: v for k, v in CATALOGUE.items() if self.active or k not in SCIENTIFIC}
        if self.labeled and self.pool_candidates:
            actions["list_ontology_pools"] = ("offset, limit. Page provisional pools of singleton seeds with the same "
                "global CL ID. marker_seed_counts counts distinct seed panels containing a gene. "
                "get_evidence(candidate_id, offset, limit) pages the full marker union and seed IDs. "
                "Use a pool ID directly to investigate membership, revise it or submit it.")
        if self.exhaustive:
            actions["list_seed_reviews"] = "offset, limit. Read written seed assessments, follow-up and retained seeds first."
        if self.adaptive:
            actions.update(adaptive.TOOLS)
            actions["get_evidence"] += " For a single expansion evidence_id, page the ranking, optionally sorted by metric."
            actions["coexpress"] += " Alternatively specify min_detected=k of the distinct genes, with optional expression threshold."
            actions["propose_revision"] += " Optional group biological_rarity (yes/no/uncertain) and rarity_reason record context-specific interpretation. Use matched_evidence_id in a group for all exact hits of one measured coexpress query. Each group can subtract exclude_candidate_ids from its union (e.g. remove a rare candidate from a broad group)."
        recent = s.get_meta("receipts", [])
        queries = s.all("CanonicalQuery")
        usage = s.usage()
        stalled_turns = s.get_meta("stall", 0)
        stalls_left = max(0, self.budget["stall_turns"] - stalled_turns)
        ctx = {"contract": contract(), "goal": self.goal, "mode": self.mode, "actions": actions,
               "bindings": {"run_id": s.get_meta("run_id"), "task_id": "discovery",
                            "data_version": self.version, "contract_version": contract()["version"]},
               "species": self.manifest.get("species"), "tissue": self.manifest.get("tissue"),
               "views": self.manifest["views"], "candidate_count": len(self.candidates),
               "selected_count": len(s.get_meta("selected", [])),
               "remaining": {"turns_including_this": max(0, self.turn_limit - turn + 1),
                             "consecutive_no_progress_turns_before_stop": stalls_left,
                             "new_queries": self.query_limit - usage["new_queries"] if self.active else 0,
                             "lineage_max_depth": self.depth_limit},
               "current_candidates": [], "recent_receipts": [], "prior_queries": [],
               "selected_ids": s.get_meta("selected", [])[:25]}
        ctx["prompt_version"] = self.prompt
        ctx["turn_instruction"] = (
            "Inspect the candidate cards and the most recent results. Choose an action that resolves a specific membership question. "
            "Previously viewed pages are available by ID. Repeating a request yields no new evidence. "
            "Use submit when you can decide or no further query is justified. Write the chosen action and its named arguments."
        )
        if turn == self.turn_limit:
            ctx["turn_instruction"] = "This is the final turn. Submit your best supported assignment with current evidence."
        elif stalled_turns:
            ctx["turn_instruction"] += (
                " The previous turn made no progress. A repeated page, query, or rejected action cannot reset the stall counter."
                " Use an unread evidence page or a distinct permitted investigation only if it could change membership."
            )
            if stalls_left == 1:
                ctx["turn_instruction"] += (
                    " One more turn without progress ends this run. Submit your best supported candidate_ids now"
                    " if no such investigation is needed. An empty selection is valid when none are supported."
                )
            if self.active and ctx["remaining"]["new_queries"] > 0:
                depths = s.get_meta("depths", {})
                selected = set(s.get_meta("selected", []))
                pending = sorted({seed["candidate_id"] for seed in s.all("Seed")
                                  if not depths.get(seed["candidate_id"]) and seed["candidate_id"] not in selected})
                if pending:
                    ctx["pending_seeds"] = {"count": len(pending), "candidate_ids": pending[:3]}
        if self.adaptive:
            ctx["rarity_search"] = {"population_n": len(self.cells), "minimum_nonempty_pct": 100 / len(self.cells),
                                  "resolvable_global_bands_pct": [p for p in (5, 1, 0.1) if p * len(self.cells) >= 100]}
            pending = sorted({r["candidate_id"] for r in s.all("Seed")} - set(s.get_meta("selected", []))
                             - {r for r, depth in s.get_meta("depths", {}).items() if depth})
            ctx["pending_seeds"] = {"count": len(pending), "candidate_ids": pending[:5]}
            if stalls_left == 1 and stalled_turns and not self.stall_review():
                ctx["turn_instruction"] = (
                    "A repeat will end exploration. Compare the held boundary options and their additional cells "
                    "in boundary_reviews. Submit the better-supported selection now, or take a distinct action."
                )
            ctx["boundary_reviews"] = []
            reviews = s.all("BoundaryReview")
            ctx["boundary_review_count"] = len(reviews)
            ctx["recent_evidence_pages"] = s.get_meta("evidence_pages", [])[-5:]
        closing = self.closing_actions(turn)
        live_active = self.active and s.get_meta("policy") == "live"
        if live_active:
            ctx.pop("pending_seeds", None)  # The ranked screen below replaces the legacy ID sample.
        if closing:
            ctx["actions"] = {k: actions[k] for k in closing}
            ctx["phase"] = "decision"
            ctx["turn_instruction"] = (
                "Exploration is closed. Submit now using exact candidate_ids, or omit the IDs to publish a saved "
                "selection. Use an explicit empty list and a reason if you decline all groups. You may revise "
                "memberships with existing evidence then submit in the same reply. New queries and browsing are closed."
            )
            ctx["reviewed_anchors"] = []
        if s.get_meta("policy") == "live":
            ctx["remaining"]["exploration_turns_including_this"] = 0 if closing else self.decision_start() - turn
            if not closing and turn + 1 == self.decision_start():
                ctx["turn_instruction"] += " This is your last exploration turn. Batch the measurements needed for your decision now."
        if live_active:
            queue = self.screened_seeds if closing else self.pending_seed_reviews()
            ctx["seed_review"] = {"screened": len(self.screened_seeds), "priority": [],
                                  "basis": "Shared signed markers and RNA views guide follow-up. These are not independent "
                                           "validation. Use their evidence IDs to propose a joint group if warranted."}
            if not closing:
                ctx["seed_review"]["pending_followup"] = len(queue)
        if self.adaptive and self.budget.get("max_reported_tokens"):
            consumed = [m["usage"].get("total_tokens") for m in s.all("Message")]
            if all(n is not None for n in consumed):
                left = self.budget["max_reported_tokens"] - sum(consumed)
                ctx["remaining"]["reported_tokens"] = max(0, left)
                if consumed and left <= consumed[-1]:
                    ctx["turn_instruction"] += (
                        " The remaining token budget is at most the cost of your last turn. "
                        "Submit your best supported selection in this response if possible.")
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
        if self.labeled and self.pool_candidates:
            ctx["ontology_pool_count"] = len(self.pool_candidates)
            visible = ctx["ontology_seed_pools"] = []
            # Reserve space for scientific receipts as well as the first-page pool index.
            available = min(limit // 4, max(0, (limit - len(json_bytes(ctx))) // 2))
            for candidate in self.pool_candidates:
                entry = self.pool_card(candidate)
                available -= len(json_bytes(entry)) + 1
                if available < 0 or not append_if_fits("ontology_seed_pools", entry):
                    break
            # Show complete dictionaries for smaller seed pools when their union fits.
            for entry in visible:
                counts = self.candidate(entry["candidate_id"])["ontology_pool"]["marker_counts"]
                if len(counts) <= len(entry["marker_seed_counts"]):
                    continue
                complete = dict(counts)
                extra = len(json_bytes(complete)) - len(json_bytes(entry["marker_seed_counts"]))
                if extra <= available:
                    entry["marker_seed_counts"] = complete
                    available -= extra
        if self.exhaustive:
            ctx["written_seed_reviews"] = []
            for item in exhaustive.page(self, 0, 8)["reviews"]:
                if not append_if_fits("written_seed_reviews", dict(item, reason=item["reason"][:240])):
                    break
        if closing:
            for item in reversed(s.all("Investigation")[-12:]):
                record = {k: item[k] for k in ("candidate_id", "evidence_id", "proposed_revision_id")}
                if not append_if_fits("reviewed_anchors", record):
                    break
        if live_active:
            for row in queue[:4]:
                refs = [self.candidate(cid)["evidence_ids"][0] for cid in (row["candidate_id"], row["peer"]) if cid]
                item = dict(row, evidence_ids=refs)
                label = self.current_label(row["candidate_id"])
                if label:
                    item["annotation_sources"] = label["sources"]
                ctx["seed_review"]["priority"].append(item)
                if len(json_bytes(ctx)) > limit:
                    ctx["seed_review"]["priority"].pop()
                    break
        links = self.expansion_links(recent) if self.active else []
        if links:
            ctx["recent_expansions"] = []
            if not append_if_fits("recent_expansions", links[0]):
                if not append_if_fits("recent_expansions", {"candidate_id": links[0]["candidate_id"]}):
                    del ctx["recent_expansions"]
                    links = []
        if self.adaptive:
            for review in reversed(reviews[-5:]):
                if not append_if_fits("boundary_reviews", adaptive.review_card(self, review)):
                    break
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
        for link in links[1:]:
            if not append_if_fits("recent_expansions", link):
                break
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
        self.checked_evidence(args.get("evidence_ids", []))
        all_scope = args.get("scope") == "all"
        members = sorted(self.cells) if all_scope else self.resolve_members(args)
        bound = {"members": members}
        # Scalar "all" in an expansion is the deployed agent's shorthand for
        # dataset-wide retrieval. An explicit list always retains exact ID semantics.
        background = args.get("background")
        search_all = action == "expand_candidate" and background == "all" and "all" not in self.cell_set
        if action == "expand_candidate" and args.get("background_mode") == "complement":
            complement = sorted(self.cell_set.difference(members))
            if not complement:
                raise ValueError("Complement reference is empty")
            bound["background"] = complement
        elif background == "all" and not search_all:
            bound["background"] = ["all"]
        elif "background_id" in args or ("background" in args and not search_all):
            bound["background"] = self.resolve_members(args, "background_id", "background")
        elif action == "compare_groups":
            complement = sorted(self.cell_set.difference(members))
            if complement:
                bound["background"] = complement
        if "background" in bound and set(members).intersection(bound["background"]):
            raise ValueError("Reference background intersects foreground")
        if "genes" in args:
            bound["genes"] = sorted(set(resolve_genes(args["genes"], self.genes, species=self.manifest.get("species"))))
        if action == "coexpress":
            bound["threshold"] = float(args.get("threshold", 0))
            if "min_detected" in args:
                if not self.adaptive or not 1 <= args["min_detected"] <= len(bound["genes"]):
                    raise ValueError("min_detected requires A3_adaptive and 1..number of distinct resolved genes")
                k = args["min_detected"]
                if k == len(bound["genes"]) or k == 1:
                    bound["logic"] = "AND" if k == len(bound["genes"]) else "OR"
                else:
                    bound["min_detected"] = k
            else:
                bound["logic"] = args.get("logic", "AND")
        canonical = {"action": action, "args": bound, "data_version": self.version,
                     "transform_version": self.snapshot_hash, "method_version": [__version__, METHOD_VERSION],
                     "parameters": self.cfg["search"], "layer": "log1p"}
        if action == "compare_groups" and bound.get("genes"):
            canonical["method_version"].append("nb-mixture/1")
        if self.adaptive and action == "expand_candidate":
            canonical["method_version"].append("ranked-expansion/1")
        signature = fingerprint(canonical)
        c = self.candidate(args["candidate_id"]) if "candidate_id" in args else None
        depths = self.state.get_meta("depths", {})
        seen = self.state.get_meta("expansions", {})
        if action == "expand_candidate":
            if c is None:
                raise ValueError("Expansion needs an exact lineage anchor")
            if signature not in seen and max(depths.get(r, 0) for r in self.lineage_roots(c)) >= self.depth_limit:
                raise ValueError("Expansion depth exhausted for this lineage")
        def compute():
            if self.science is None:
                self.science = Science(self.snapshot, self.cfg, adaptive=self.adaptive)
            return self.science.execute(action, bound)
        result, key = self.state.query(canonical, self.query_limit, compute)
        ref = self.state.result_ref(key)
        try:
            with self.state.transaction():
                response = self._register_result(action, args, bound, c, signature, result, key, ref)
                if action == "expand_candidate":
                    response.update(anchor_id=c["id"], search_scope="all",
                                    background_rule=result["selection_history"]["background_rule"])
                return response
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
        if args.get("evidence_ids"):
            use["cited_evidence_ids"] = sorted(set(args["evidence_ids"]))
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
            if depth >= self.depth_limit:
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

    def make_revision(self, members, parents, label, status, evidence_ids, kind, basis, biological_rarity=None):
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
        if biological_rarity:
            descriptor.append(biological_rarity)
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
            interpretation_summary = {"reason": basis, "proposed_label": label, "status": status, "size": len(members)}
            if biological_rarity:
                interpretation_summary["biological_rarity"] = biological_rarity
            interpretation = self.evidence(cid, "interpretation", evidence_ids, "discovery_agent",
                                           sh, interpretation_summary, [],
                                           history=[{"parent_candidates": parents}])
            evidence_ids.add(interpretation)
        record = {"id": cid, "members": members, "set_hash": sh, "parents": parents, "roots": sorted(roots),
                  "initial": False, "sources": [{"kind": kind}], "evidence_ids": sorted(set(evidence_ids)),
                  "proposed_label": label, "status": status, "markers": []}
        if biological_rarity:
            record["biological_rarity"] = biological_rarity
        if kind == "ontology_seed_pool":
            counts = sorted(basis["marker_seed_counts"].items(), key=lambda item: (-item[1], item[0]))
            record["ontology_pool"] = {"ontology_id": basis["ontology_id"], "marker_counts": counts}
        if kind == "boundary":
            record["boundary_rule"] = basis
        self.save_candidate(record)
        return record

    def checked_evidence(self, ids):
        for eid in ids:
            try:
                self.state.get("Evidence", eid)
            except KeyError:
                raise ValueError("Unknown evidence ID " + eid) from None
        return ids

    def decide(self, action, args, message_id, *, source="agent"):
        if self.exhaustive and action == "submit" and exhaustive.coverage(self)["remaining"]:
            raise ValueError("Submission requires a measurement and written assessment for every seed")
        selected = set(self.state.get_meta("selected", []))
        evidence_ids = self.checked_evidence(args.get("evidence_ids", []))
        reason = args.get("reason", "Submit current provisional assignment")
        if action == "submit":
            selected = set(args.get("candidate_ids", selected))
            if self.adaptive:
                anchors = [self.candidate(cid)["parents"][0] for cid in selected
                           if self.candidate(cid).get("boundary_rule")]
                if len(anchors) != len(set(anchors)):
                    raise ValueError("Choose one boundary per anchor or revise the alternatives into disjoint groups")
            for cid in selected:
                evidence_ids = evidence_ids + self.candidate(cid)["evidence_ids"]
                evidence_ids = evidence_ids + self.current_label(cid).get("evidence_ids", [])
        else:
            groups = args["groups"]
            parent_ids = args.get("parents")
            if not parent_ids:
                if selected or len(groups) != 1 or not groups[0].get("candidate_ids"):
                    raise ProtocolError("Specify parents to supersede when a selection exists, or when sources are ambiguous")
                parent_ids = groups[0]["candidate_ids"]
            parents = [self.candidate(cid)["id"] for cid in parent_ids]
            reason = args.get("reason") or "\n".join(dict.fromkeys(group["reason"] for group in groups))
            chosen = []
            occupied = set()
            measured_evidence = set()
            revision_inputs = set(parents)
            for group in groups:
                group_reason = group.get("reason", reason)
                group_evidence = sorted(set(evidence_ids) | set(self.checked_evidence(group.get("evidence_ids", []))))
                if not group_evidence:
                    raise ValueError("Each revision group requires registered evidence")
                refs = [self.candidate(cid) for cid in group.get("candidate_ids", [])]
                revision_inputs.update(c["id"] for c in refs)
                member_set = {m for c in refs for m in c["members"]}
                explicit = set(group.get("members", []))
                if explicit and not self.active and not self.by_set.get(fingerprint(sorted(explicit))):
                    raise ValueError("Passive mode can select or combine only precomputed member sets")
                rarity = None
                if "biological_rarity" in group or "rarity_reason" in group:
                    if not self.adaptive:
                        raise ValueError("biological rarity metadata requires A3_adaptive")
                    rarity = {"assessment": group.get("biological_rarity", "uncertain"),
                              "reason": group.get("rarity_reason", group_reason)}
                matched = group.get("matched_evidence_id")
                if matched:
                    if not self.adaptive:
                        raise ValueError("matched_evidence_id requires A3_adaptive")
                    measured = adaptive.coexpression_members(self, matched)
                    member_set.update(measured)
                    measured_evidence.add(matched)
                member_set |= explicit
                exclusions = group.get("exclude_candidate_ids", [])
                if exclusions and not self.adaptive:
                    raise ValueError("exclude_candidate_ids requires A3_adaptive")
                for excluded in exclusions:
                    other = self.candidate(excluded)
                    member_set.difference_update(other["members"])
                if occupied.intersection(member_set):
                    raise ValueError("New groups in a single revision overlap")
                occupied.update(member_set)
                inherited = set(parents) | {c["id"] for c in refs} | set(exclusions)
                chosen.append(self.make_revision(member_set, inherited, group.get("proposed_label", "unknown"),
                                                 group.get("status", "provisional"), group_evidence + ([matched] if matched else []),
                                                 "revision", group_reason, biological_rarity=rarity)["id"])
            if self.adaptive:
                # Revising from a held alternative also supersedes the old working
                # boundary of that anchor. An exclusion alone leaves its draft intact.
                anchors = {self.candidate(i)["parents"][0] for i in revision_inputs
                           if self.candidate(i).get("boundary_rule")}
                selected = {i for i in selected if not (self.candidate(i).get("boundary_rule") and
                                                        self.candidate(i)["parents"][0] in anchors)}
            selected.difference_update(parents)
            selected.update(chosen)
            evidence_ids = list(evidence_ids) + sorted(measured_evidence) + [eid for cid in chosen for eid in self.candidate(cid)["evidence_ids"]]
        decision = {"action": action, "reason": reason, "evidence_ids": sorted(set(evidence_ids)),
                    "selected": sorted(selected), "message_id": message_id, "args_hash": fingerprint(args)}
        if action == "propose_revision" and not args.get("parents"):
            decision["parents_from_candidate_ids"] = parents
        if action == "submit":
            decision["status"] = args.get("status", "final")
            if self.adaptive:
                decision["assignment_preview"] = adaptive.preview(self, {"candidate_ids": sorted(selected)})
            if self.state.get_meta("policy") == "live" or self.adaptive:
                decision["source"] = source
                self.state.set_meta("submission_source", source)
            if evidence_ids and source == "agent":
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
        if "metric" in args:
            ids = args.get("evidence_ids", [])
            if not self.adaptive or len(ids) != 1 or not any(r["evidence_id"] == ids[0] for r in self.state.all("Investigation")):
                raise ValueError("metric requires one A3_adaptive expansion evidence_id")
        offset = args.get("offset", 0)
        limit = min(args.get("limit", 8), 25)
        if action == "list_candidates":
            ordered = self.ordered()
            return {"total": len(ordered), "offset": offset, "candidates": [self.card(c) for c in ordered[offset:offset+limit]]}
        if "candidate_id" in args:
            c = self.candidate(args["candidate_id"])
            page = {"candidate": self.card(c), "offset": offset, "total": len(c["members"]),
                    "members": c["members"][offset:offset+limit],
                    "markers": compact(c.get("markers", [])[offset:offset+limit], limit)}
            if "ontology_pool" in c:
                page["seed_candidate_ids"] = c["parents"][offset:offset+limit]
                page["marker_seed_counts"] = dict(c["ontology_pool"]["marker_counts"][offset:offset+limit])
            if self.adaptive:
                reviews = [r for r in self.state.all("BoundaryReview")
                           if c["id"] in (r["anchor_id"], r["requested"]["candidate_id"], r.get("wider", {}).get("candidate_id"))]
                if reviews:
                    page["boundary_review_count"] = len(reviews)
                    page["boundary_reviews"] = [adaptive.review_card(self, r) for r in reversed(reviews[-3:])]
            return page
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
                    if "ranking" in data:
                        ranked = data["ranking"]
                        if "metric" in args:
                            metric = args["metric"]
                            ranked = sorted((r for r in ranked if r[metric] is not None),
                                            key=lambda r: (-round(r[metric], 12), self.barcodes[r["cell_id"]]))
                        sections["ranking"] = ranked
                    elif "metric" in args:
                        raise ValueError("metric requires adaptive expansion ranking evidence")
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
        if self.labeled and action == "list_ontology_pools":
            start = args.get("offset", 0)
            size = min(25, args.get("limit", 8))
            return {"total": len(self.pool_candidates), "offset": start,
                    "ontology_seed_pools": [self.pool_card(c) for c in self.pool_candidates[start:start+size]]}
        if self.exhaustive and action == "review_seed":
            return exhaustive.record(self, message_id, args)
        if self.exhaustive and action == "list_seed_reviews":
            return exhaustive.page(self, args.get("offset", 0), min(25, args.get("limit", 8)))
        if action in SCIENTIFIC:
            return self.scientific(action, args)
        if action == "preview_submission" and self.adaptive:
            return adaptive.preview(self, args)
        if action in {"list_candidates", "get_evidence"}:
            return self.read_tool(action, args)
        if action in {"propose_revision", "submit"} or (self.adaptive and action in {"apply_boundary", "select_boundary"}):
            try:
                with self.state.transaction():
                    if action == "apply_boundary":
                        return adaptive.apply_boundary(self, args)
                    if action == "select_boundary":
                        return adaptive.select_boundary(self, args)
                    return self.decide(action, args, message_id)
            except BaseException:
                self._load_candidates()
                raise
        raise ValueError("Unregistered tool")
