"""Produce a versioned total partition from immutable membership claims."""
import csv
import io
from collections import defaultdict

from . import exhaustive, pool_followup
from .config import atomic, fingerprint, hash_file, json_bytes, write_json
from .discovery import rows


def tsv(path, columns, records):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(records)
    atomic(path, stream.getvalue())


def jsonlines(path, records):
    atomic(path, b"".join(json_bytes(row) + b"\n" for row in records))


def membership_claims(ctrl, selected):
    claims = defaultdict(set)
    for cid in selected:
        c = ctrl.candidate(cid)
        if c["set_hash"] != fingerprint(sorted(set(c["members"]))):
            raise ValueError("Selected candidate version failed " + cid)
        for cell in c["members"]:
            claims[cell].add(cid)
    return claims


def assignments(ctrl):
    selected = ctrl.state.get_meta("selected", [])
    resolution = None
    if ctrl.adaptive:
        from .claim_resolution import plan
        resolution = plan(ctrl, selected)
        claims = {cell: {cid} for cell, cid in resolution["owners"].items()}
        claims.update(resolution["unresolved"])
    else:
        claims = membership_claims(ctrl, selected)
    baselines = rows(ctrl.snapshot / "baselines.tsv")
    recipes = {(b["view"], b["resolution"]) for b in baselines}
    priority = ["V_HVG", "V_LINEAR", "V_GINI", "V_TENSOR2"]
    base = min(recipes, key=lambda k: (priority.index(k[0]), float(k[1]), k[1])) if recipes else None
    fallback = {b["cell_id"]: b["group_id"] for b in baselines if (b["view"], b["resolution"]) == base}
    if base and set(fallback) != ctrl.cell_set:
        raise ValueError("Incomplete fallback partition")
    provisional = set(resolution["provisional_cells"]) if resolution else set()
    partition, conflicts = [], []
    for cell in sorted(ctrl.cells):
        rivals = sorted(claims.get(cell, ()))
        if len(rivals) == 1:
            cid = rivals[0]
            c = ctrl.candidate(cid)
            gid = "claim/" + cid
            label = c["proposed_label"]
            if label == "unknown" and c["initial"]:
                label = ctrl.current_label(cid).get("label", "unknown")
            status, reason = c["status"], "selected"
            if cell in provisional:
                status, reason = "provisional", "overlap_priority"
            if resolution and resolution["alternative_labels"].get(cid) and cid not in resolution["preferred"]:
                label, status = "unknown", "provisional"
        else:
            cid, label, status = "", "unknown", "unresolved"
            reason = "competing_claims" if rivals else "unclaimed"
            gid = ("residual/%s/%s/%s" % (base[0], base[1], fallback[cell])) if base else "unresolved/" + cell
            if rivals:
                conflicts.append({"cell_id": cell, "candidate_ids": rivals, "decision": "unresolved_baseline",
                                  "evidence_ids": sorted({e for id_ in rivals for e in ctrl.candidate(id_)["evidence_ids"]})})
        partition.append(dict(cell_id=cell, group_id=gid, candidate_id=cid, proposed_label=label,
                              status=status, assignment_reason=reason))
    ids = [r["cell_id"] for r in partition]
    if len(ids) != len(set(ids)) or set(ids) != ctrl.cell_set:
        raise ValueError("Partition violates complete/disjoint assignment")
    groups = defaultdict(list)
    for row in partition:
        groups[row["group_id"]].append(row["cell_id"])
        if row["candidate_id"] and row["candidate_id"] not in claims[row["cell_id"]]:
            raise ValueError("Unexpected partition membership")
    return partition, conflicts, [{"group_id": k, "member_count": len(v), "cell_set_hash": fingerprint(sorted(v))}
                                   for k, v in sorted(groups.items())], base


def export(ctrl, status, detail, policy, annotation_status):
    run, s = ctrl.run, ctrl.state
    part, conflicts, groups, base = assignments(ctrl)
    tsv(run / "partition.tsv", ["cell_id", "group_id", "candidate_id", "proposed_label", "status", "assignment_reason"], part)
    tsv(run / "groups.tsv", ["group_id", "member_count", "cell_set_hash"], groups)
    jsonlines(run / "conflicts.jsonl", conflicts)
    records = sorted(ctrl.candidates.values(), key=lambda c: c["id"])
    tsv(run / "candidate_members.tsv", ["candidate_id", "cell_id"],
        (dict(candidate_id=c["id"], cell_id=cell) for c in records for cell in c["members"]))
    jsonlines(run / "candidates.jsonl", ({k: v for k, v in c.items() if k != "members"} for c in records))
    annotations = s.all("Annotation")
    tsv(run / "annotations.tsv", ["candidate_id", "seed_id", "background_type", "n_background", "label", "status", "reason",
                                  "genes", "ontology_id", "partition_id", "marker_hash", "cache_key", "evidence_id"],
        (dict(row, genes=",".join(row.get("selection_history", {}).get("genes", []))) for row in annotations))
    jsonlines(run / "decisions.jsonl", s.all("Decision"))
    jsonlines(run / "evidence.jsonl", s.all("Evidence"))
    seeds = s.all("Seed")
    selected = set(s.get_meta("selected", []))
    depths = s.get_meta("depths", {})
    final = {r["cell_id"]: r for r in part}
    screened = {r["candidate_id"]: r for r in ctrl.screened_seeds} if ctrl.active and s.get_meta("policy") == "live" else {}
    followups = {r["anchor_id"]: r for r in s.all("SeedFollowup")}
    seed_rows = []
    for seed in seeds:
        cid = seed["candidate_id"]
        activity = "investigated" if depths.get(cid, 0) else "selected" if cid in selected else "pending"
        if ctrl.routed and s.get_optional("SeedRoute", cid):
            activity = "reviewed"
        seed_rows.append(dict(seed, seed_id=seed["id"], barcode=seed.get("barcode", seed["cell_id"]), task_status=activity,
                              screened=cid in screened, priority_peer=screened.get(cid, {}).get("peer"),
                              shared_markers=",".join(screened.get(cid, {}).get("shared_markers", [])),
                              controller_followup=followups.get(cid, {}).get("result", {}).get("status"), final_group=final[seed["cell_id"]]["group_id"]))
    tsv(run / "seeds.tsv", ["seed_id", "cell_id", "barcode", "candidate_id", "task_status", "final_group", "screened", "priority_peer", "shared_markers", "controller_followup"], seed_rows)
    # One durable queue covers all base candidates, including seeds left untouched
    # when the policy submits early or reaches its budget.
    investigated = {row["candidate_id"] for row in s.all("Investigation")}
    cited_parents = {cid for item in records if not item["initial"] for cid in item["parents"]}
    def queue_rows():
        for item in records:
            if not item["initial"]:
                continue
            cid = item["id"]
            if cid in selected:
                status = "selected"
            elif cid in investigated:
                status = "expansion_checked"
            elif cid in cited_parents:
                status = "revised"
            else:
                status = "pending"
            yield {"task_id": "review/" + cid, "candidate_id": cid, "cell_set_hash": item["set_hash"], "status": status}
    tsv(run / "candidate_tasks.tsv", ["task_id", "candidate_id", "cell_set_hash", "status"], queue_rows())
    # All baseline memberships come from the immutable prepare snapshot.
    atomic(run / "artifacts/baselines.tsv", (ctrl.snapshot / "baselines.tsv").read_bytes())
    messages = s.all("Message")
    dispatch = s.all("ProviderDispatch")
    counters = s.usage()
    counters.update(turns=s.get_meta("finished_turn", 0), provider_dispatches=len(dispatch),
                    protocol_errors=s.get_meta("protocol_errors", 0), action_errors=s.get_meta("action_errors", 0),
                    deterministic_format_repairs=s.get_meta("format_repairs", 0), semantic_parser_calls=0,
                    resume_count=s.get_meta("resume_count", 0), discovery_elapsed_seconds=s.get_meta("elapsed_seconds", 0.0),
                    api_usage=[m["usage"] for m in messages if "usage" in m], provider_cost_currency=None)
    for k in ("input_tokens", "output_tokens", "total_tokens"):
        values = [m["usage"][k] for m in messages if k in m.get("usage", {})]
        counters["reported_" + k] = sum(values) if values else None
        counters[k] = sum(values) if dispatch and len(values) == len(dispatch) else None
    usage_path = run / "usage.json"
    write_json(usage_path, counters)
    summary = {"mode": ctrl.mode, "policy": policy, "status": status, "detail": detail,
               "submitted": s.get_meta("submitted", False), "submission_status": s.get_meta("submission_status"),
               "annotation_status": annotation_status, "prepare_status": ctrl.manifest.get("status"),
               "stages": ctrl.manifest.get("stages", {}), "snapshot_hash": ctrl.snapshot_hash,
               "data_version": ctrl.version, "n_cells": len(part), "n_candidates": len(records), "n_groups": len(groups),
               "conflicted_cells": len(conflicts), "unresolved_cells": sum(r["status"] == "unresolved" for r in part),
               "issue": s.get_meta("current_issue"), "n_pending_seeds": sum(r["task_status"] == "pending" for r in seed_rows),
               "fallback_view_resolution": list(base) if base else None,
               "replication_id": s.get_meta("replication_id"), "run_id": s.get_meta("run_id"),
               "resumable": detail in ("paused", "provider_blocked", "annotation_blocked", "unknown_execution", "temporary_failure", "incomplete")}
    if policy == "live" and ctrl.active:
        summary.update(n_screened_seeds=len(screened), n_controller_seed_followups=len(followups),
                       n_expanded_seeds=sum(seed["candidate_id"] in investigated for seed in seeds))
    if s.get_meta("code_upgrade"):
        summary["code_upgrade"] = s.get_meta("code_upgrade")
        summary["code_upgrades"] = s.get_meta("code_upgrades", [summary["code_upgrade"]])
    if s.get_meta("closeout_trigger"):
        summary["closeout_trigger"] = s.get_meta("closeout_trigger")
    if ctrl.adaptive or policy == "live":
        summary["submission_source"] = s.get_meta("submission_source")
    if ctrl.exhaustive:
        summary["seed_review_coverage"] = exhaustive.coverage(ctrl)
        measurements = {r["candidate_id"]: r for r in s.all("SeedMeasurement")}
        notes = {r["candidate_id"]: r for r in s.all("SeedAssessment")}
        def review_rows():
            for seed in seed_rows:
                cid = seed["candidate_id"]
                measurement, note = measurements.get(cid, {}), notes.get(cid, {})
                yield {**seed, "measurement_status": measurement.get("status", "pending"),
                       "query_key": measurement.get("query_key"), "evidence_id": measurement.get("evidence_id"),
                       "measurement_error": measurement.get("error"), "proposal_id": measurement.get("proposal_id"),
                       "message_id": note.get("message_id"), "disposition": note.get("disposition"),
                       "guide_label": note.get("guide_label"), "label_reliability": note.get("label_reliability"),
                       "reason": note.get("reason")}
        tsv(run / "seed_reviews.tsv", ["seed_id", "cell_id", "candidate_id", "measurement_status", "query_key",
            "evidence_id", "measurement_error", "proposal_id", "message_id", "disposition", "guide_label", "label_reliability", "reason", "final_group"], review_rows())
    if ctrl.exhaustive and ctrl.labeled:
        summary["pool_followup_coverage"] = pool_followup.coverage(ctrl)
        pool_rows = []
        for c, notes in pool_followup.queue(ctrl):
            row = s.get_optional("PoolFollowup", c["id"]) or {"candidate_id": c["id"], "label": c["proposed_label"], "status": "pending"}
            pool_rows.append(dict(row, follow_up_seed_ids=[r["candidate_id"] for r in notes]))
        for cid, case in pool_followup.independent_cases(ctrl).items():
            row = s.get_optional("PoolFollowup", cid) or dict(case, status="pending")
            pool_rows.append(dict(row, follow_up_seed_ids=[case["anchor_id"]]))
        reviews = {r["candidate_id"]: r for r in s.all("ExpansionReview")}
        for row in pool_rows:
            if row["candidate_id"] in reviews:
                row["expansion_review"] = reviews[row["candidate_id"]]
        jsonlines(run / "pool_followups.jsonl", pool_rows)
    if ctrl.routed:
        from .routed import RoutedDiscovery
        route_files = {"seed_routes.jsonl": s.all("SeedRoute"), "route_jobs.jsonl": RoutedDiscovery(ctrl).jobs(),
                       "route_claims.jsonl": s.all("RouteFreeze"), "route_resolutions.jsonl": s.all("RouteResolution")}
        for name, data in route_files.items():
            jsonlines(run / name, data)
        forced_ref = s.get_meta("forced_choice_ref")
        if forced_ref:
            detail_record = s.read_artifact(forced_ref)
            summary["forced_choice"] = dict(artifact_ref=forced_ref, n_claimed=detail_record["n_claimed"],
                                             comparisons=len(detail_record["pair_votes"]),
                                             n_wins_tied=detail_record["n_wins_tied"])
        resolutions = route_files["route_resolutions.jsonl"]
        deferred = sum(row["choice"] == "defer_to_submission" for row in resolutions)
        summary["routing"] = {"reviewed_seeds": len(route_files["seed_routes.jsonl"]),
                              "name_searches": sum(j["kind"] == "name" for j in route_files["route_jobs.jsonl"]),
                              "seed_searches": sum(j["kind"] == "seed" for j in route_files["route_jobs.jsonl"]),
                              "frozen": len(route_files["route_claims.jsonl"]),
                              "reviews": len(resolutions) - deferred, "deferred_pairs": deferred}
    write_json(run / "run_summary.json", summary)
    files = ["run_summary.json", "usage.json", "manifest.json", "config_resolved.yaml", "project_contract.yaml",
             "partition.tsv", "groups.tsv", "conflicts.jsonl", "candidate_tasks.tsv", "candidate_members.tsv", "candidates.jsonl", "seeds.tsv",
             "annotations.tsv", "decisions.jsonl", "evidence.jsonl"]
    if ctrl.exhaustive:
        files.append("seed_reviews.tsv")
    if ctrl.exhaustive and ctrl.labeled:
        files.append("pool_followups.jsonl")
    if ctrl.routed:
        files.extend(route_files)
    files.extend(str(p.relative_to(run)) for p in sorted((run / "artifacts").glob("*")) if p.is_file())
    freeze = {"version": 1, "snapshot_hash": ctrl.snapshot_hash, "files": {name: hash_file(run / name) for name in files}}
    write_json(run / "frozen.json", freeze)
    return summary
