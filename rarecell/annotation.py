"""Read-only OntoAnno projection of existing Seurat partitions."""
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
from collections import defaultdict

from .config import atomic as _write, fingerprint as _hash, hash_file as _digest
from .config import json_bytes as _encode, read_json as _read
from .providers import safe_url
from .runtime import rscript_command

PIN = "ceda74c7d74af67a02093b07dd87623d846a59fe"
PROMPT = "ontoanno-parent-cl-v1"
SCRIPT = Path(__file__).parent / "r" / "annotation.R"
CONTRACT = Path(__file__).parent / "resources" / "PROJECT_CONTRACT.yaml"


def _options(cfg):
    options = cfg.get("annotation", cfg)
    if options is False or options == "off":
        return {"enabled": False}
    if not isinstance(options, dict):
        raise ValueError("annotation must be a mapping or off")
    return options


def _resolve(cfg):
    options = _options(cfg)
    provider = options.get("provider", "openai")
    if provider not in ("openai", "anthropic", "gemini", "ollama", "vllm"):
        raise ValueError("unsupported OntoAnno HTTP provider")
    top = options.get("top_genes", 10)
    if type(top) is not int or not 1 <= top <= 100:
        raise ValueError("top_genes must be between 1 and 100")
    batch_size = options.get("batch_size", 1)
    stop_line = options.get("single_line_stop", False)
    if type(stop_line) is not bool or (stop_line and batch_size != 1):
        raise ValueError("single_line_stop requires batch_size=1")
    if type(batch_size) is not int or not 1 <= batch_size <= 50:
        raise ValueError("annotation.batch_size must be between 1 and 50")
    model = options.get("model")
    if not isinstance(model, str) or not model.strip() or model.startswith("<"):
        raise ValueError("configure annotation.model to an available deployed model")
    env = options.get("api_key_env", {"openai": "OPENAI_API_KEY", "anthropic": "ANTHROPIC_API_KEY",
                                      "gemini": "GEMINI_API_KEY", "vllm": "VLLM_API_KEY"}.get(provider))
    if env and (not isinstance(env, str) or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", env)):
        raise ValueError("annotation.api_key_env must name an environment variable")
    if provider in ("openai", "anthropic", "gemini") and not env:
        raise ValueError("configure the annotation credential environment variable")
    # ellmer 0.5.2 requires VLLM_API_KEY even on loopback. A null config
    # selects the upstream default name, never a hard-coded secret.
    if provider == "vllm" and not env:
        env = "VLLM_API_KEY"
    if provider == "ollama" and env:
        raise ValueError("OntoAnno ollama requires api_key_env: null")
    if provider == "ollama" and ":" in model:
        raise ValueError("this OntoAnno commit strips Ollama tags; use a unique bare model alias")
    url = options.get("base_url")
    if url:
        if provider not in ("ollama", "vllm"):
            raise ValueError("OntoAnno custom base_url is supported for vllm and ollama")
        url = safe_url(url)
    timeout = float(options.get("timeout", 180))
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("annotation.timeout must be positive and finite")
    return {"provider": provider, "model": model, "api_key_env": env, "base_url": url,
            "top_genes": top, "timeout": timeout, "batch_size": batch_size,
            "single_line_stop": stop_line}


def diagnose(cfg):
    options = _options(cfg)
    if not options.get("enabled", True):
        return {"status": "not_run", "detail": "annotation explicitly off"}
    msgs = []
    try:
        settings = _resolve(cfg)
        env = settings["api_key_env"]
        if env and not os.environ.get(env):
            msgs.append("missing credential environment variable " + env)
    except (ValueError, TypeError) as exc:
        msgs.append(str(exc))
    command = rscript_command(options)
    r = shutil.which(command)
    versions = None
    if not r:
        msgs.append("Rscript executable not found")
    else:
        try:
            probe = subprocess.run([r, "--vanilla", str(SCRIPT), "--probe"], capture_output=True,
                                   text=True, timeout=90)
            if probe.returncode:
                msgs.append("OntoAnno version/import probe failed: " + probe.stderr.strip()[:400])
            else:
                versions = probe.stdout.strip()
        except (OSError, subprocess.TimeoutExpired):
            msgs.append("Rscript probe failed or timed out")
    return {"status": "blocked" if msgs else "configured", "detail": " | ".join(msgs),
            "rscript": r or command, "versions": versions, "network_test": "not_run"}


def _groups(snapshot):
    candidates = {}
    for line in (snapshot / "initial.jsonl").read_text().splitlines():
        row = json.loads(line)
        for source in row["sources"]:
            if source.get("view") == "V_HVG" and source.get("kind") == "partition":
                pair = (float(source["resolution"]), row["id"])
                candidates[pair] = {"candidate_id": row["id"], "resolution": pair[0],
                                    "partition_id": source["partition_id"],
                                    "n_members": len(row["members"]),
                                    "cell_set_hash": _hash(sorted(row["members"]))}
    marker_rows = defaultdict(list)
    if candidates:
        with (snapshot / "seurat" / "markers.tsv").open(newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                pair = (float(row["resolution"]), row["group_id"])
                if pair not in candidates:
                    raise ValueError("Seurat marker row lacks matching partition membership")
                marker_rows[pair].append(row)
    return candidates, marker_rows


def _jobs(snapshot, manifest, settings):
    candidates, markers = _groups(snapshot)
    by_resolution, skipped = defaultdict(list), []
    for key, info in sorted(candidates.items()):
        rows = markers[key]
        ranked = sorted([r for r in rows if float(r["avg_log2FC"]) > 0],
                        key=lambda r: (float(r["p_val_adj"]), -float(r["avg_log2FC"]), r["gene_id"]))
        # Gene IDs remain in the marker hash. The model gets available gene symbols.
        genes = list(dict.fromkeys(r.get("symbol") or r["gene_id"] for r in ranked))
        group = dict(info, marker_hash=_hash(rows), genes=genes[:settings["top_genes"]])
        if info["n_members"] == 1 or not genes:
            skipped.append(dict(group, status="not_annotatable",
                                reason="singleton" if info["n_members"] == 1 else "no_positive_markers"))
        else:
            by_resolution[key[0]].append(group)
    system_prompt = ("Provide one tentative cell type name per marker row, in the same order. "
                     "Use unknown when the markers cannot support a name. Return only those names.\n"
                     + CONTRACT.read_text())
    head = {"commit": PIN, "prompt_version": PROMPT, "adapter_sha256": _digest(SCRIPT),
            "contract_sha256": _digest(CONTRACT), "snapshot_hash": manifest["snapshot_hash"],
            "data_version": manifest["data_version"], "settings": settings, "system_prompt": system_prompt,
            "tissue": " ".join(str(manifest.get(k) or "") for k in ("species", "tissue")).strip()}
    jobs = []
    for res, group_list in sorted(by_resolution.items()):
        # Exact batch composition/order affects the upstream prompt and cache key.
        size = settings["batch_size"]
        for start in range(0, len(group_list), size):
            batch = group_list[start:start + size]
            note = ("\nAnnotation task: each comma-separated list of genes describes a single cluster. "
                    "There are {} cluster(s) in this request. Output exactly {} line(s), one tentative "
                    "cell type per cluster. If evidence is insufficient, answer unknown.").format(len(batch), len(batch))
            jobs.append(dict(head, resolution=res, groups=batch, system_prompt=system_prompt + note))
    return jobs, skipped


def _validate_receipt(path, job):
    record = _read(path)
    expected = [g["candidate_id"] for g in job["groups"]]
    if (record.get("commit") != PIN or record.get("valid") is not True
            or [r.get("candidate_id") for r in record["rows"]] != expected
            or not all(isinstance(r["label"], str) and r["label"].strip() for r in record["rows"])
            or not isinstance(record.get("raw_response"), str) or not record["raw_response"].strip()):
        raise ValueError("OntoAnno response failed partition/line alignment")
    for group, row in zip(job["groups"], record["rows"]):
        tokens = set(re.findall(r"[a-z0-9_.-]+", row["label"].casefold()))
        genes = {g.casefold() for g in group["genes"]}
        if len(tokens) > 1 and tokens <= genes:
            raise ValueError("OntoAnno model echoed marker names as a label")
    return record


def _batch(out, job, options):
    key = _hash(job)
    dst = out / "batches" / key
    dst.mkdir(parents=True, exist_ok=True)
    req = dst / "request.json"
    receipt = dst / "response.json"
    attempt = dst / "dispatch.json"
    if req.exists() and _read(req) != job:
        raise ValueError("annotation request content mismatch")
    if not req.exists():
        _write(req, _encode(job))
    cache_hit = receipt.exists()
    if not cache_hit:
        previous = _read(attempt)["attempts"] if attempt.exists() else 0
        if previous and (not options.get("retry_unknown", False) or previous >= 2):
            raise RuntimeError("remote execution unknown; inspect batch; one explicit retry available")
        _write(attempt, _encode({"status": "unknown_execution", "attempts": previous + 1}))
        try:
            # Upstream errors can contain sensitive transport details. Only the
            # validated, structured output crosses the process boundary.
            result = subprocess.run([rscript_command(options), "--vanilla", str(SCRIPT),
                                     str(req), str(receipt)], stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL, timeout=job["settings"]["timeout"])
        except subprocess.TimeoutExpired:
            raise RuntimeError("annotation deadline expired; remote execution unknown") from None
        if result.returncode or not receipt.exists():
            failure = dst / "response.json.failure.json"
            reason = _read(failure)["reason"] if failure.exists() else "R process returned without a response"
            raise RuntimeError("OntoAnno execution failed: " + reason)
    response = _validate_receipt(receipt, job)
    rows = []
    for group, normalized in zip(job["groups"], response["rows"]):
        fields = {k: v for k, v in group.items() if k != "genes"}
        fields.update(normalized, proposed_label=normalized["label"], status="annotated",
                      source="OntoAnno::gptcelltype", kind="interpretation",
                      data_version=job["data_version"], author_commit=PIN, model=job["settings"]["model"],
                      prompt_version=PROMPT, raw_response_ref=str(receipt.relative_to(out)),
                      cache_key=key, selection_history={"method": "Seurat positive DE", "genes": group["genes"]})
        rows.append(fields)
    return rows, cache_hit, {str(p.relative_to(out)): _digest(p) for p in [req, receipt]}


def _run(snapshot, out, cfg, manifest):
    options = _options(cfg)
    base = {"snapshot_hash": manifest["snapshot_hash"], "source_commit": PIN,
            "provider_calls_observed_at_wire": None, "annotation_cost_usd": None}
    def finish(rows, status, **extra):
        payload = b"".join(_encode(row) + b"\n" for row in rows)
        _write(out / "annotations.jsonl", payload)
        report = dict(base, status=status, annotation_sha256=hashlib.sha256(payload).hexdigest(), **extra)
        _write(out / "annotation_status.json", _encode(report))
        return report

    cfg_error = None
    try:
        settings = _resolve(cfg)
    except (ValueError, TypeError) as exc:
        settings, cfg_error = {}, str(exc)
    old = out / "annotation_status.json"
    if old.exists():
        prior = _read(old)
        if prior["snapshot_hash"] != manifest["snapshot_hash"]:
            raise ValueError("use a separate annotation directory for the new prepare snapshot")
        if prior["status"] == "success":
            if not options.get("enabled", True) or prior["settings"] != settings:
                raise ValueError("use a separate annotation directory for changed annotation settings")
            verify_annotation(snapshot, out, _manifest=manifest)
            return prior
    if not options.get("enabled", True):
        return finish([], "not_run", detail="annotation explicitly off")
    if cfg_error:
        return finish([], "blocked", detail=cfg_error)
    if "V_HVG" not in manifest["views"]:
        return finish([], "blocked", detail="V_HVG Seurat partition unavailable")
    marker_stage = manifest.get("stages", {}).get("seurat_markers", {})
    if marker_stage.get("status") != "ok":
        return finish([], "blocked", detail=marker_stage.get("reason", "Seurat marker stage unavailable"))
    jobs, skipped = _jobs(snapshot, manifest, settings)
    if not jobs:
        return finish([], "blocked", detail="no eligible Seurat marker lists", skipped=skipped)
    system = diagnose(cfg)
    rows, errors, hashes, hits, dispatches = [], [], {}, 0, 0
    for job in jobs:
        key = _hash(job)
        receipt = out / "batches" / key / "response.json"
        if system["status"] != "configured" and not receipt.exists():
            errors.append({"key": key, "reason": system["detail"]})
            continue
        try:
            batch_rows, hit, files = _batch(out, job, options)
            rows.extend(batch_rows)
            hashes.update(files)
            hits += int(hit)
            dispatches += int(not hit)
        except (ValueError, OSError, RuntimeError) as exc:
            errors.append({"key": key, "reason": str(exc)})
    status = ("partial" if rows else "blocked") if errors else "success"
    return finish(rows, status, settings=settings, skipped=skipped, errors=errors,
                  cached_batches=hits, newly_completed_batches=dispatches,
                  files=hashes, annotated=len(rows), versions=system["versions"])


def annotate(snapshot, out, cfg):
    from .prepare import verify_snapshot
    snapshot, out = Path(snapshot).resolve(), Path(out).resolve()
    if out == snapshot or snapshot in out.parents or out in snapshot.parents:
        raise ValueError("keep the annotation directory outside the prepare snapshot")
    manifest = verify_snapshot(snapshot)
    out.mkdir(parents=True, exist_ok=True)
    with (out / ".lock").open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        return _run(snapshot, out, cfg, manifest)


def verify_annotation(snapshot, out, *, allow_partial=False, _manifest=None):
    from .prepare import verify_snapshot
    snapshot, out = Path(snapshot), Path(out)
    manifest = _manifest if _manifest is not None else verify_snapshot(snapshot)
    try:
        state = _read(out / "annotation_status.json")
    except OSError as exc:
        raise ValueError("missing annotation_status.json; run annotate") from exc
    status = state["status"]
    if status != "success" and not (allow_partial and status == "partial"):
        raise ValueError("annotation is " + status)
    if state["snapshot_hash"] != manifest["snapshot_hash"] or state["source_commit"] != PIN:
        raise ValueError("annotation source/snapshot version mismatch")
    jobs, _ = _jobs(snapshot, manifest, state["settings"])
    by_key = {_hash(job): job for job in jobs}
    errors = state.get("errors", [])
    failed = {error["key"] for error in errors}
    if (len(failed) != len(errors) or not failed.issubset(by_key)
            or bool(failed) != (status == "partial")):
        raise ValueError("inconsistent annotation batch outcome records")
    expected = {}
    files = {}
    for key, job in by_key.items():
        if key in failed:
            continue
        for filename in ("request.json", "response.json"):
            name = "batches/{}/{}".format(key, filename)
            path = (out / name).resolve()
            if out.resolve() not in path.parents:
                raise ValueError("annotation artifact outside annotation directory")
            files[name] = _digest(path)
        if _read(out / "batches" / key / "request.json") != job:
            raise ValueError("annotation request mismatch")
        record = _validate_receipt(out / "batches" / key / "response.json", job)
        for group, reply in zip(job["groups"], record["rows"]):
            expected[(key, group["candidate_id"])] = (group, reply)
    if state.get("files") != files:
        raise ValueError("annotation artifact modified or batch references inconsistent")
    path = out / "annotations.jsonl"
    if _digest(path) != state["annotation_sha256"]:
        raise ValueError("annotation projection modified")
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    if not rows or len(rows) != len(expected) or state.get("annotated") != len(rows):
        raise ValueError("incomplete annotation projection")
    for row in rows:
        entry = expected.pop((row["cache_key"], row["candidate_id"]), None)
        if entry is None:
            raise ValueError("unexpected annotation row")
        group, reply = entry
        if (any(row[k] != group[k] for k in ("cell_set_hash", "marker_hash", "partition_id"))
                or row["proposed_label"] != reply["label"] or row["label"] != reply["label"]
                or row["ontology_id"] != reply["ontology_id"]):
            raise ValueError("annotation projection/receipt mismatch")
    return rows
