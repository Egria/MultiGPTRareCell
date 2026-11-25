"""Seed triage, isolated expression searches and reconciliation of frozen memberships."""

from itertools import combinations
import math

from . import providers
from .config import fingerprint, json_bytes
from .discovery import ExpansionDepthExhausted
from .protocol import MAX_ACTIONS, MAX_CHARS, _decode, _materialize, _json_region, _trim_trailing_commas
from .rule_search import ExpressionSearch
from .state import BudgetExhausted, IntegrityError, QueryBlocked

RULES = (
    'A rule is "gene" (detected), ["gene",cutoff] (log1p expression >= cutoff), '
    '{"otsu":"gene"} (high side of a positive-expression Otsu split), {"k":[2,["A","B","C"]]}, '
    '{"nb_high":"gene"} (local two-component negative-binomial high state on raw counts, empty if unsupported), '
    '{"all":[rules]}, {"any":[rules]} or {"not":rule}. {"scope":"alias"} is allowed only for listed scopes.'
)
TOOLS = {
    "survey": "genes; returns detection count, positive expression quantiles at 10/50/75/90/95/99/100% and positive Otsu split.",
    "views": "rules (1..12), genes; returns group sizes and [in detected, in mean, outside detected, outside mean] per gene.",
    "markers": "rule, optional exclude genes and limit (1..40); genome-wide mean expression differences and detection fractions.",
}
FREEZE = (
    "rule (null for abstention), support (0..100), reason, kind (lineage or state). "
    "Give a proposed_label on seed searches; unknown is allowed. A score describes RNA support, not a probability."
)
GENERIC = {"unknown", "unresolved", "uncertain", "unassigned", "none", "n/a", "unclassified"}


def label(value):
    if not isinstance(value, str):
        return ""
    value = " ".join(value.split())
    if len(value) > 160:
        raise ValueError("Hypothesis name exceeds 160 characters")
    return "" if value.casefold() in GENERIC else value


def required_text(data, key):
    value = data.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Supply a nonempty {key}")
    return value.strip()


def support_score(value):
    if type(value) not in (float, int) or not math.isfinite(value) or not 0 <= value <= 100:
        raise ValueError("Supply a numerical RNA support score from 0 to 100")
    return value


def response(raw):
    if len(raw) > MAX_CHARS:
        raise ValueError("Model action too large")
    data = _materialize(_decode(_trim_trailing_commas(_json_region(raw))))
    if isinstance(data, dict) and set(data) == {"actions"}:
        data = data["actions"]
    items = data if isinstance(data, list) else [data]
    if not 1 <= len(items) <= MAX_ACTIONS or not all(isinstance(item, dict) for item in items):
        raise ValueError("Supply one action object or a nonempty list of query objects")
    for item in items:
        for wrapper in ("args", "arguments"):
            if wrapper in item:
                args = item.pop(wrapper)
                if not isinstance(args, dict) or item.keys() & args.keys() or {"action", "op", "operation"} & args.keys():
                    raise ValueError("Action arguments must be an object with unambiguous fields")
                item.update(args)
        names = [item[key] for key in ("action", "op", "operation") if key in item]
        if names:
            if not all(isinstance(name, str) and name == names[0] for name in names):
                raise ValueError("Specify a single unambiguous action name")
            item["action"] = names[0]
    return items[0] if len(items) == 1 else items


class RoutedDiscovery:
    def __init__(self, ctrl):
        self.ctrl = ctrl
        self.s = ctrl.state
        self.engine = None
        self.terms = {}
        for ann in ctrl.annotations:
            name = label(ann.get("label"))
            term = ann.get("ontology_id")
            if name and ann.get("status") == "annotated" and isinstance(term, str) and term.startswith("CL:"):
                self.terms.setdefault(name.casefold(), set()).add(term)

    def key(self, name):
        terms = self.terms.get(name.casefold(), set())
        return next(iter(terms)) if len(terms) == 1 else name.casefold()

    def record_route(self, row, expected):
        cid = row.get("candidate_id")
        if cid not in expected:
            raise ValueError("Review only an assigned seed")
        reliability = row.get("reliability", row.get("label_reliability"))
        if reliability == "unassessed":
            reliability = "uncertain"
        if reliability not in ("supported", "uncertain", "contradicted", "unknown"):
            raise ValueError("Invalid seed reliability")
        guide = label(row.get("guide_label"))
        if reliability == "supported" and not guide:
            raise ValueError("A supported seed requires a guide_label")
        rec = dict(
            candidate_id=cid, guide_label=guide, reliability=reliability, reason=required_text(row, "reason")
        )
        self.s.put("SeedRoute", cid, rec)

    def jobs(self):
        """Rebuild the queue from immutable seed assessments and completed unnamed searches."""
        names, anchored = {}, []

        def nominate(name, parent, priority):
            name = label(name)
            if not name:
                return
            k = self.key(name)
            job = names.setdefault(
                k, dict(id="name/" + k, kind="name", name=name, parents=set(), priority=priority)
            )
            if parent:
                job["parents"].add(parent)
            job["priority"] = min(priority, job["priority"])

        for name in self.ctrl.cfg.get("potential_names", []):
            nominate(name, None, 0)
        for row in self.s.all("SeedRoute"):
            cid, guide = row["candidate_id"], row["guide_label"]
            supported = row["reliability"] == "supported"
            nominate(guide, cid, 0 if supported else 1)
            for ann in self.ctrl.labels.get(cid, {}).values():
                if ann.get("status") == "annotated":
                    nominate(ann.get("label"), cid, 2)
            if not supported:
                anchored.append(dict(id="seed/" + cid, kind="seed", name=guide, parents=[cid], priority=1))
        for result in self.s.all("RouteFreeze"):
            if result["job"].startswith("seed/") and result.get("candidate_id"):
                nominate(result["label"], result["candidate_id"], 1)
        for job in names.values():
            job["parents"] = sorted(job["parents"])
        # Give low-confidence names (including neighbor-only nominations) an opportunity
        # even when a large unknown-seed queue consumes the remaining global budget.
        primary = sorted((j for j in names.values() if j["priority"] == 0), key=lambda j: j["id"])
        secondary = sorted(
            (j for j in names.values() if j["priority"] > 0), key=lambda j: (j["priority"], j["id"])
        )
        anchored.sort(key=lambda j: j["id"])
        for i in range(max(len(secondary), len(anchored))):
            primary.extend(secondary[i : i + 1])
            primary.extend(anchored[i : i + 1])
        return primary

    def triage_context(self, pending):
        seeds = []
        context = dict(
            phase="seed_triage",
            seeds=seeds,
            actions={
                "route_seeds": "reviews: one object per seed with candidate_id, guide_label, reliability (supported|uncertain|contradicted|unknown), reason."
            },
            instruction="Audit each individual seed against its RNA markers and proposed names. Use biological knowledge. A concrete alternative can be worth a whole-data search even when its seed is unreliable.",
        )
        for cid in pending[: self.ctrl.budget["seed_review_batch_size"]]:
            c = self.ctrl.candidate(cid)
            seeds.append(
                dict(
                    candidate_id=cid,
                    markers=c.get("markers", [])[:10],
                    annotations=self.ctrl.current_label(cid).get("sources", []),
                )
            )
            if len(json_bytes(context)) > self.ctrl.budget["context_bytes"] - 512 and len(seeds) > 1:
                seeds.pop()
                break
        return context

    def scope_members(self, job):
        if job["kind"] == "name":
            return {}
        if job["kind"] == "pair":
            left, right = [set(self.ctrl.candidate(i)["members"]) for i in job["parents"]]
            return {
                k: sorted(v)
                for k, v in dict(
                    left=left,
                    right=right,
                    shared=left & right,
                    left_only=left - right,
                    right_only=right - left,
                ).items()
            }
        cid = job["parents"][0]
        anchor = self.ctrl.candidate(cid)["members"]
        saved = self.s.get_optional("RouteRetrieval", cid)
        if not saved:
            try:
                result = self.ctrl.scientific("expand_candidate", {"candidate_id": cid})
            except ExpansionDepthExhausted:
                saved = dict(candidate_id=cid, status="depth_exhausted", depth_limit=self.ctrl.depth_limit)
            else:
                raw = self.s.read_artifact(result["artifact_ref"])
                saved = dict(
                    candidate_id=cid,
                    ref=result["artifact_ref"],
                    evidence_id=result["evidence_id"],
                    n_added=len(raw.get("added_members", [])),
                    status=result["status"],
                )
            self.s.put("RouteRetrieval", cid, saved)
        if saved["status"] == "depth_exhausted":
            return {"seed": anchor}
        measured = self.s.read_artifact(saved["ref"])
        return {"seed": anchor, "retrieved": sorted(set(anchor) | set(measured.get("added_members", [])))}

    def search_context(self, job, scopes, history):
        pair = job["kind"] == "pair"
        context = dict(
            phase="boundary_review"
            if pair
            else "blind_name_search"
            if job["kind"] == "name"
            else "seed_expansion",
            actions=dict(TOOLS),
            rules=RULES,
            population_size=len(self.ctrl.cells),
            expression_layer="log1p of normalized expression",
            target=job["name"],
        )
        if pair:
            context.update(
                drafts=[self.claim_card(cid) for cid in job["parents"]],
                instruction="Review shared and exclusive cells, including disjoint same-name fragments. Query separate markers for each flank. Merge coherent groups, retain distinct hypotheses or allocate only the contested cells. Reason from measured evidence.",
            )
            context["actions"]["resolve"] = (
                "choice: merge|left|right|separate|split_overlap; reason. Supply "
                "support (0..100) on a merge and a proposed_label when merging different labels. For split_overlap supply left_support and right_support (0..100 for each revised flank) and optional left_rule and right_rule on shared cells. Cells satisfying both or neither remain unresolved. "
                "Separate retains disjoint groups; for intersecting groups select a boundary explicitly."
            )
        else:
            context["actions"]["freeze"] = FREEZE
            if self.ctrl.cfg.get("routed_conflict_policy") == "single":
                context["actions"]["freeze"] += " Include label_reliability (supported or uncertain) for identity independent of boundary precision."
            context["instruction"] = (
                "Use biological knowledge to choose genes and evaluate competing explanations. "
                "Compare expression intensity and joint detection across the whole sample, inspect neighboring cutoffs, "
                "and check genes outside the selection rule. Measure plausible alternate boundaries via views, including "
                "with and without a competing-marker veto. On raw counts consider a local nb_high state. "
                "Previously measured overlapping alternatives can be reviewed for unclaimed hits. Freeze exact membership or abstain."
            )
        if scopes:
            context["scopes"] = {k: len(v) for k, v in scopes.items()}
        if job["kind"] == "seed":
            seed = job["parents"][0]
            retrieval = self.s.get("RouteRetrieval", seed)
            if retrieval["status"] == "depth_exhausted":
                context["retrieval"] = dict(status="depth_exhausted", depth_limit=retrieval["depth_limit"],
                    reason="No structural retrieval was run. Inspect the original seed and expression rules, or abstain.")
            else:
                measured = self.s.read_artifact(retrieval["ref"])
                context["retrieval"] = {
                    "summary": measured["summary"],
                    "signature": measured.get("selection_history", {}).get("signature", [])[:12],
                }
            context.update(
                seed_review=self.s.get("SeedRoute", seed),
                seed_markers=self.ctrl.candidate(seed).get("markers", [])[:12],
                annotations=self.ctrl.current_label(seed).get("sources", []),
            )
        context["history"] = list(history)
        context["remaining_queries"] = max(
            0, self.ctrl.budget["queries_per_route"] - sum("query_key" in h for h in history)
        )
        if not context["remaining_queries"]:
            context["actions"] = {k: v for k, v in context["actions"].items() if k not in TOOLS}
        return context

    def claim_card(self, cid):
        c = self.ctrl.candidate(cid)
        return dict(
            candidate_id=cid,
            label=c["proposed_label"],
            n=len(c["members"]),
            evidence=c.get("route_evidence", {}),
        )

    def next_task(self):
        pending = sorted(set(self.ctrl.seed_ids) - {r["candidate_id"] for r in self.s.all("SeedRoute")})
        if pending:
            return None, self.triage_context(pending)
        frozen = {r["job"] for r in self.s.all("RouteFreeze")}
        job = next((j for j in self.jobs() if j["id"] not in frozen), None)
        if job is None and self.ctrl.cfg.get("routed_conflict_policy") == "single":
            return None, None
        return job or next(self.pending_pairs(), None), None

    def pending_pairs(self):
        reviewed = {r["job"] for r in self.s.all("RouteResolution")}
        for a, b in combinations(sorted(self.s.get_meta("route_active", [])), 2):
            jid = "pair/" + a + "/" + b
            if jid in reviewed:
                continue
            left, right = self.ctrl.candidate(a), self.ctrl.candidate(b)
            same = label(left["proposed_label"]) and self.key(left["proposed_label"]) == self.key(right["proposed_label"])
            if same or not set(left["members"]).isdisjoint(right["members"]):
                yield dict(id=jid, kind="pair", name="Compare submitted boundaries", parents=[a, b])

    def defer_pair(self, job, reason, mid=None):
        # Keep both original versions active. Final assignment will use the normal
        # overlap policy on the surviving drafts; no new biological decision is inferred.
        self.s.put("RouteResolution", job["id"], dict(
            job=job["id"], choice="defer_to_submission", source="controller", reason=reason,
            message_id=mid, retained=job["parents"],
        ))

    def measured(self, req, scopes):
        op = req["action"]
        fields = {
            "survey": ("genes",),
            "views": ("rules", "genes"),
            "markers": ("rule", "exclude", "limit"),
            "freeze": ("rule",),
        }
        query = {k: req[k] for k in fields[op] if k in req}
        if op == "markers" and type(query.get("limit")) is int and query["limit"] > 40:
            query["limit"] = 40
        query["action"] = op
        canonical = dict(
            method="routed-expression/1",
            snapshot=self.ctrl.snapshot_hash,
            query=query,
            scopes={k: fingerprint(v) for k, v in scopes.items()},
        )

        def compute():
            if self.engine is None:
                self.engine = ExpressionSearch(self.ctrl)
            return self.engine.measure(query, scopes)

        result, key = self.s.query(canonical, self.ctrl.query_limit, compute)
        return result, key

    def active_change(self, remove, add):
        active = set(self.s.get_meta("route_active", []))
        self.s.set_meta("route_active", sorted((active - set(remove)) | set(add)))

    def draft(self, members, parents, name, kind, basis, evidence):
        c = self.ctrl.make_revision(members, parents, name, "provisional", evidence, "routed_claim", basis)
        return c["id"]

    def freeze(self, job, req, members, key, mid):
        name = job["name"] if job["kind"] == "name" else label(req.get("proposed_label")) or "unknown"
        kind = req.get("kind", "lineage")
        support = support_score(req.get("support"))
        if kind not in ("lineage", "state"):
            raise ValueError("Supply kind lineage or state")
        note = dict(
            job=job["id"],
            label=name,
            kind=kind,
            support=support,
            reason=required_text(req, "reason"),
            rule=req["rule"],
            n=len(members),
            message_id=mid,
            query_key=key,
        )
        if self.ctrl.cfg.get("routed_conflict_policy") == "single" and "label_reliability" in req:
            reliability = req["label_reliability"]
            if reliability not in ("supported", "uncertain"):
                raise ValueError("label_reliability must be supported or uncertain")
            note["label_reliability"] = reliability
        cid = None
        if members:
            sh = fingerprint(sorted(members))
            ref = self.s.result_ref(key)
            refs = [ref] + [
                self.s.result_ref(r["exchange"]["query_key"])
                for r in self.s.all("RouteExchange")
                if r["job"] == job["id"] and "query_key" in r["exchange"]
            ]
            eid = self.ctrl.evidence(
                job["id"], "observation", [], "input_expression", sh, note, list(dict.fromkeys(refs))
            )
            parents = job["parents"]
            cid = self.draft(members, parents, name, kind, note, [eid])
            if kind == "lineage":
                self.active_change([], [cid])
        self.s.put("RouteFreeze", job["id"], dict(note, candidate_id=cid))
        return dict(action="freeze", candidate_id=cid, n=len(members), support=support)

    def resolve(self, job, req, scopes, splits, mid, split_keys=()):
        left, right = job["parents"]
        choice = req.get("choice")
        why = required_text(req, "reason")
        a, b = self.ctrl.candidate(left), self.ctrl.candidate(right)
        evidence = sorted(set(a["evidence_ids"] + b["evidence_ids"]))
        basis = dict(
            job=job["id"], source="discovery_agent",
            choice=choice,
            reason=why,
            message_id=mid,
            parent_support={c["id"]: c.get("route_evidence", {}).get("support") for c in (a, b)},
        )
        basis["query_keys"] = [
            r["exchange"]["query_key"] for r in self.s.all("RouteExchange") if r["job"] == job["id"]
        ] + list(split_keys)
        if choice == "split_overlap":
            basis["allocation_rules"] = {side: req.get(side + "_rule") for side in ("left", "right")}
        if choice == "merge":
            basis["support"] = support_score(req.get("support"))
            if "proposed_label" in req:
                name = label(req["proposed_label"]) or "unknown"
            else:
                name = (
                    a["proposed_label"]
                    if self.key(a["proposed_label"]) == self.key(b["proposed_label"])
                    else "unknown"
                )
            chosen = [
                self.draft(
                    set(scopes["left"]) | set(scopes["right"]),
                    [left, right],
                    name,
                    "lineage",
                    basis,
                    evidence,
                )
            ]
        elif choice in ("left", "right"):
            chosen = [left if choice == "left" else right]
        elif choice == "separate" and not scopes["shared"]:
            chosen = [left, right]
        elif choice == "split_overlap":
            intersection = set(scopes["shared"])
            chosen = []
            # Restrict the decision to disputed cells; never remove exclusive flanks.
            for side, other, c in (("left", "right", a), ("right", "left", b)):
                own = (splits[side] - splits[other]) & intersection
                members = set(scopes[side + "_only"]) | own
                if members:
                    chosen.append(
                        c["id"]
                        if members == set(c["members"])
                        else self.draft(
                            members,
                            [c["id"], right if side == "left" else left],
                            c["proposed_label"],
                            "lineage",
                            dict(
                                basis,
                                side=side,
                                rule=req.get(side + "_rule"),
                                support=support_score(req.get(side + "_support")),
                            ),
                            evidence,
                        )
                    )
            if not chosen and intersection and set(scopes["left"]) == set(scopes["right"]):
                # Identical submissions agree on membership even if both names are uncertain.
                basis["alternative_labels"] = [a["proposed_label"], b["proposed_label"]]
                basis["support"] = 0  # No agreed identity for these otherwise identical members.
                chosen = [self.draft(intersection, [left, right], "unknown", "lineage", basis, evidence)]
            else:
                basis["unresolved_members"] = sorted(intersection - (splits["left"] ^ splits["right"]))
        else:
            raise ValueError("Choose a registered resolution; separate requires disjoint membership")
        self.active_change([left, right], chosen)
        self.s.put("RouteResolution", job["id"], dict(basis, retained=chosen))
        return dict(action="resolve", retained=chosen)

    def apply(self, job, ctx, scopes, req, mid):
        if isinstance(req, list):
            self.apply_query_batch(job, ctx, scopes, req, mid)
            return
        op = req.get("action", req.get("op"))
        if not isinstance(op, str) or op not in ctx["actions"]:
            raise ValueError("Choose a listed action for the current task")
        req = dict(req, action=op)
        if op == "route_seeds":
            reviews = req.get("reviews")
            if not isinstance(reviews, list) or not reviews or not all(isinstance(r, dict) for r in reviews):
                raise ValueError("Expected nonempty reviews array")
            expected = {r["candidate_id"] for r in ctx["seeds"]}
            ids = [r.get("candidate_id") for r in reviews]
            if not all(isinstance(i, str) for i in ids) or len(ids) != len(set(ids)):
                raise ValueError("Invalid or repeated seed assessment")
            with self.s.transaction():
                for row in reviews:
                    self.record_route(row, expected)
                self.complete_turn(mid, {"action": op, "reviewed": ids})
            return
        if op in TOOLS:
            result, key = self.measured(req, scopes)
            exchange = dict(request=req, response=result, query_key=key)
            with self.s.transaction():
                self.s.put("RouteExchange", mid, dict(job=job["id"], exchange=exchange))
                self.complete_turn(mid, {"action": op, "query_key": key, "summary": result})
        elif op == "freeze":
            if "rule" not in req:
                raise ValueError("Supply rule or null for abstention")
            members, key = [], None
            if req["rule"] is not None:
                result, key = self.measured(req, scopes)
                members = result["members"]
                if members and job["kind"] == "seed" and set(members).isdisjoint(scopes.get("retrieved", scopes["seed"])):
                    raise ValueError("Anchored claim must overlap the seed or its retrieved cells")
            with self.s.transaction():
                receipt = self.freeze(job, req, members, key, mid)
                self.complete_turn(mid, receipt)
        else:
            splits = {}
            split_keys = []
            for side in ("left", "right"):
                rule = req.get(side + "_rule")
                splits[side] = set()
                if rule is not None and req.get("choice") == "split_overlap":
                    out, key = self.measured(dict(action="freeze", rule=rule), scopes)
                    split_keys.append(key)
                    splits[side] = set(out["members"])
            with self.s.transaction():
                self.complete_turn(mid, self.resolve(job, req, scopes, splits, mid, split_keys))

    def apply_query_batch(self, job, ctx, scopes, items, mid):
        ops = [item.get("action", item.get("op")) for item in items]
        if not all(isinstance(op, str) and op in ctx["actions"] for op in ops):
            raise ValueError("Every item must use a listed action for this task")
        queries = items[:-1] if ops[-1] in {"resolve", "freeze"} else items
        if not queries or any(op not in TOOLS for op in ops[:len(queries)]):
            raise ValueError("Send one decision per reply. Only read-only queries can be batched.")
        room = min(ctx["remaining_queries"], self.ctrl.budget["max_actions_per_turn"])
        if not room:
            raise ValueError("No queries remaining; send one decision")
        exchanges = []
        for item, op in zip(queries[:room], ops):
            request = dict(item, action=op)
            result, key = self.measured(request, scopes)
            exchanges.append(dict(request=request, response=result, query_key=key))
        deferred = ops[len(exchanges):]
        if deferred:
            exchanges[-1]["not_executed"] = deferred
        with self.s.transaction():
            for i, exchange in enumerate(exchanges):
                self.s.put("RouteExchange", f"{mid}/{i}", dict(job=job["id"], exchange=exchange))
            self.complete_turn(mid, dict(action="query_batch", query_keys=[e["query_key"] for e in exchanges],
                                         not_executed=deferred))

    def complete_turn(self, mid, receipt):
        self.s.put(
            "Receipt",
            mid + "/0",
            dict(id=mid + "/0", message_id=mid, result=receipt, progress="error" not in receipt, metadata={}),
        )
        self.s.set_meta("finished_turn", int(mid.rsplit("T", 1)[1]))
        self.s.set_meta("route_last_error", None)
        self.s.event(
            "action",
            {"message_id": mid, "action": receipt.get("action"), "status": receipt.get("status", "ok")},
        )

    def run(self, obtain, pause_after, retry_unknown):
        self.s.set_meta("policy", "live")
        budget = self.ctrl.budget
        while True:
            job, ctx = self.next_task()
            if job is None and ctx is None:
                if self.ctrl.cfg.get("routed_conflict_policy") == "single":
                    from .forced_routing import finish
                    return finish(self, obtain, pause_after, retry_unknown)
                with self.s.transaction():
                    deferred = any(r["choice"] == "defer_to_submission" for r in self.s.all("RouteResolution"))
                    args = dict(candidate_ids=self.s.get_meta("route_active", []),
                                reason="Publish frozen lineages after boundary review.")
                    if deferred:
                        args.update(status="provisional", reason="Apply the existing overlap policy to surviving drafts, including pairs deferred by the controller.")
                    self.ctrl.decide("submit", args, None, source="controller_routed" if deferred else "agent")
                return "complete", "submitted"
            turn = self.s.get_meta("finished_turn", 0) + 1
            mid = f"{self.s.get_meta('run_id')}:T{turn:04d}"
            stopping = None
            if turn > budget["max_agent_turns"]:
                stopping = "routed_turn_budget"
            elif budget.get("max_reported_tokens") and self.s.get_optional("Message", mid) is None:
                counts = [m["usage"].get("total_tokens") for m in self.s.all("Message")]
                if any(v is None for v in counts) or sum(counts) >= budget["max_reported_tokens"]:
                    stopping = "token_usage_unreported" if None in counts else "max_reported_tokens"
            if stopping:
                if job and job["kind"] == "pair":
                    with self.s.transaction():
                        for pair in self.pending_pairs():
                            self.defer_pair(pair, stopping)
                    continue
                self.s.set_meta("current_issue", "Routed tasks pending: " + stopping)
                return "partial", stopping
            scopes = {}
            if job:
                try:
                    scopes = self.scope_members(job)
                except (BudgetExhausted, QueryBlocked) as exc:
                    self.s.set_meta("current_issue", str(exc))
                    return "partial", "routed_query_budget" if isinstance(
                        exc, BudgetExhausted
                    ) else "routed_retrieval_blocked"
                history = [r["exchange"] for r in self.s.all("RouteExchange") if r["job"] == job["id"]]
                ctx = self.search_context(job, scopes, history)
            ctx.update(tissue=self.ctrl.manifest.get("tissue"), species=self.ctrl.manifest.get("species"))
            task = job["id"] if job else "triage"
            last = self.s.get_meta("route_last_error")
            if last and last["task"] == task:
                ctx["previous_error"] = last["error"]
            while len(json_bytes(ctx)) > budget["context_bytes"] and ctx.get("history"):
                ctx["history"].pop(0)
            if len(json_bytes(ctx)) > budget["context_bytes"]:
                if job and job["kind"] == "pair":
                    self.defer_pair(job, "Pair context exceeds context_bytes")
                    continue
                self.s.set_meta("current_issue", "One routed task exceeds budget.context_bytes")
                return "partial", "context_blocked"
            try:
                msg = obtain(self.ctrl, "live", turn, retry_unknown, context=ctx)
            except providers.ProviderError as exc:
                self.s.set_meta("current_issue", str(exc))
                return "partial", "provider_blocked" if exc.status == "blocked" else exc.status
            try:
                if not msg.get("complete", True):
                    raise ValueError("Incomplete reply; send one complete JSON action")
                self.apply(job, ctx, scopes, response(msg["raw"]), mid)
            except BudgetExhausted as exc:
                if job and job["kind"] == "pair":
                    with self.s.transaction():
                        self.defer_pair(job, str(exc), mid)
                        self.complete_turn(mid, dict(action="defer_to_submission", reason=str(exc)))
                else:
                    self.s.set_meta("current_issue", str(exc))
                    return "partial", "routed_query_budget"
            except IntegrityError:
                raise
            except (ValueError, KeyError, TypeError, RecursionError) as exc:
                self.ctrl._load_candidates()  # Roll back any speculative in-memory revision.
                error = str(exc)[:300]
                failures = last["attempts"] + 1 if last and last["task"] == task else 1
                exhausted = failures >= max(2, budget["stall_turns"])
                defer = exhausted and job is not None and job["kind"] == "pair"
                with self.s.transaction():
                    self.complete_turn(mid, dict(error=error, status="deferred" if defer else "needs_clarification"))
                    self.s.set_meta("action_errors", self.s.get_meta("action_errors", 0) + 1)
                    if defer:
                        self.defer_pair(job, error, mid)
                    else:
                        self.s.set_meta("route_last_error", dict(task=task, error=error, attempts=failures))
                if exhausted and not defer:
                    self.s.set_meta("current_issue", f"Routed task {task}: {error}")
                    return "partial", "routed_action_blocked"
            if pause_after is not None and turn >= pause_after:
                return "partial", "paused"
