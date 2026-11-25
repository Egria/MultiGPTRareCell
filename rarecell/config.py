"""Immutable run configuration and content checks for the packaged contract."""
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path

import yaml

ADAPTIVE_MODE = "A3_adaptive"
ROUTED_MODE = "A3_routed"
MODES = ("A0", "A1", "A2", "A3", ADAPTIVE_MODE, ROUTED_MODE)
ACTIVE_MODES = {"A2", "A3", ADAPTIVE_MODE, ROUTED_MODE}
ANNOTATED_MODES = {"A1", "A3", ADAPTIVE_MODE, ROUTED_MODE}

RESOURCES = Path(__file__).parent / "resources"


def merge_config(base, override):
    output = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(output.get(key), dict):
            output[key] = merge_config(output[key], value)
        else:
            output[key] = value
    return output


def validate_settings(cfg):
    from .providers import safe_url
    names = {
        "provider": {"kind", "base_url", "endpoint", "api_style", "model", "api_key_env", "timeout", "options"},
        "annotation": {"enabled", "provider", "model", "base_url", "api_key_env", "top_genes", "timeout", "rscript", "retry_unknown", "batch_size", "single_line_stop", "include_singletons", "dual_backgrounds", "panel_mode"},
    }
    annotation = cfg.get("annotation")
    if isinstance(annotation, dict) and "panel_mode" in annotation and annotation["panel_mode"] not in ("global", "dual"):
        raise ValueError("annotation.panel_mode must be global or dual")
    for flag in ("include_singletons", "dual_backgrounds"):
        if isinstance(cfg.get("annotation"), dict) and type(cfg["annotation"].get(flag, True)) is not bool:
            raise ValueError(f"annotation.{flag} must be true or false")
    if type(cfg.get("exhaustive_seed_review", False)) is not bool:
        raise ValueError("exhaustive_seed_review must be true or false")
    hypotheses = cfg.get("potential_names", [])
    if not isinstance(hypotheses, list) or not all(isinstance(n, str) and 0 < len(n.strip()) <= 160 for n in hypotheses):
        raise ValueError("potential_names must contain nonempty names of at most 160 characters")
    if set(cfg) - {"potential_names", "exhaustive_seed_review"} != {"mode", "prepare", "search", "budget", "provider", "annotation"}:
        raise ValueError("Unknown or missing top-level configuration section")
    for name, fields in names.items():
        section = cfg[name]
        if not isinstance(section, dict) or set(section) - fields:
            raise ValueError("Unrecognized setting in " + name + "; credentials belong in environment variables")
        if section.get("base_url"):
            safe_url(section["base_url"])
        key = section.get("api_key_env")
        if key is not None and not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", key):
            raise ValueError("Supply the credential environment variable name in " + name)
    options = cfg["provider"].get("options", {})
    allowed = {"temperature", "seed", "top_p", "max_tokens", "max_completion_tokens", "max_output_tokens", "reasoning", "reasoning_effort"}
    if not isinstance(options, dict) or set(options) - allowed:
        raise ValueError("Unrecognized provider options")
    if "reasoning" in options and (not isinstance(options["reasoning"], dict) or set(options["reasoning"]) - {"effort", "summary"}):
        raise ValueError("Unknown provider reasoning options")
    base = yaml.safe_load((RESOURCES / "default.yaml").read_text())
    for section in ("prepare", "search", "budget"):
        if not isinstance(cfg[section], dict) or set(cfg[section]) - set(base[section]):
            raise ValueError("Unknown settings in " + section)
    return cfg


def config(path=None, override=None):
    base = yaml.safe_load((RESOURCES / "default.yaml").read_text())
    user = yaml.safe_load(Path(path).read_text()) if path else {}
    if (override or {}).get("mode", (user or {}).get("mode", base["mode"])) == ROUTED_MODE:
        # Routed budgets cover triage, every independent search and reconciliation together.
        base["budget"].update(max_agent_turns=800, max_new_queries=800, seed_review_batch_size=8)
    combined = merge_config(merge_config(base, user or {}), override or {})
    validate_settings(combined)
    if combined["mode"] not in MODES:
        raise ValueError("mode must be one of " + ", ".join(MODES))
    if any(type(v) is not int or v < 0 for v in combined["budget"].values()):
        raise ValueError("budgets must be nonnegative integers")
    if combined["budget"]["stall_turns"] < 1 or combined["budget"]["max_actions_per_turn"] < 1:
        raise ValueError("stall_turns and max_actions_per_turn must be positive")
    if not 1 <= combined["budget"]["seed_review_batch_size"] <= 128:
        raise ValueError("seed_review_batch_size must be between 1 and 128")
    return combined


def contract():
    obj = yaml.safe_load((RESOURCES / "PROJECT_CONTRACT.yaml").read_text())
    if set(obj["requirements"]) != {"C%02d" % n for n in range(1, 18)}:
        raise ValueError("Incomplete discovery contract")
    if any(len(row) != 3 or not all(row) for row in obj["requirements"].values()):
        raise ValueError("Each contract clause requires text, modules and checks")
    return obj


def fingerprint(obj):
    return hashlib.sha256(json_bytes(obj)).hexdigest()


def json_bytes(obj):
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def hash_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as tmp:
        try:
            tmp.write(data if isinstance(data, bytes) else data.encode())
            tmp.flush()
            os.fsync(tmp.fileno())
        except BaseException:
            os.unlink(tmp.name)
            raise
    os.replace(tmp.name, path)


def write_json(path, obj):
    atomic(path, json_bytes(obj) + b"\n")


def read_json(path):
    return json.loads(Path(path).read_bytes())


def read_lines(path):
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)
