"""Command-line workflows for preparation, discovery and isolated scoring."""
import argparse
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

from . import annotation, providers
from .config import config, merge_config, read_json, write_json
from .runtime import checkout, load_keys, model_overrides, project_model, rscript_command, save_project_model
from .prepare import prepare, verify_snapshot
from .run import discover, resume


def doctor(cfg):
    command = rscript_command(cfg["prepare"])
    r = shutil.which(command)
    probe = {"status": "blocked", "detail": "Rscript not found"}
    if r:
        expr = ('cat(as.character(getRversion()),as.character(packageVersion("Seurat")),'
                'as.character(packageVersion("SeuratObject")),as.character(packageVersion("leidenbase")),sep="/")')
        try:
            result = subprocess.run([r, "--vanilla", "-e", expr], capture_output=True, text=True, timeout=90)
            versions = result.stdout.strip().split("/")
            compatible = result.returncode == 0 and len(versions) == 4 and versions[1:3] == ["5.3.0", "5.2.0"]
            probe = {"status": "available" if compatible else "blocked",
                     "R/Seurat/SeuratObject/leidenbase": result.stdout.strip(),
                     "detail": result.stderr.strip()[-500:] if compatible or result.returncode else
                     "Validated pairing requires Seurat 5.3.0 and SeuratObject 5.2.0"}
        except (OSError, subprocess.TimeoutExpired) as exc:
            probe = {"status": "blocked", "detail": str(exc)}
    probe["rscript"] = r or command
    return {"python": platform.python_version(), "packages": {p: importlib.metadata.version(p) for p in
             ("numpy", "scipy", "scikit-learn", "igraph", "PyYAML")}, "Seurat": probe,
            "discovery_provider": providers.diagnose(cfg), "OntoAnno": annotation.diagnose(cfg)}


def demo(root, cfg, mode, seed):
    from .evaluation import evaluate
    from .synthetic import generate
    root = Path(root).resolve()
    if root.exists() and any(root.iterdir()):
        raise ValueError("demo needs an empty destination")
    root.mkdir(exist_ok=True, parents=True)
    dataset = generate(root, seed)
    if mode == "core":
        cfg = merge_config(cfg, {"mode": "A2", "prepare": {"seed": seed, "graph_backend": "igraph",
                            "views": ["V_LINEAR", "V_TENSOR2"]}, "annotation": {"enabled": False},
                            "provider": {"kind": "scripted_demo"}})
    else:
        cfg = merge_config(cfg, {"mode": "A3", "prepare": {"seed": seed}})
    prepare(dataset, root / "snapshot", cfg)
    # Evaluation fixtures live outside both discovery input and run state.
    annotation.annotate(root / "snapshot", root / "annotation", cfg)
    report = discover(root / "snapshot", root / "run", root / "cache", cfg,
                      policy="scripted_demo" if mode == "core" else "live",
                      annotations=root / "annotation", replication_id="demo-%s-%s" % (mode, seed),
                      pause_after=1 if mode == "core" else None)
    if report["detail"] == "paused":
        report = resume(root / "run")
    metrics = evaluate(root / "run", root / "evaluation_fixture/gt.tsv", root / "evaluation")
    return {"demo": mode, "root": str(root), "summary": report, "evaluation": str(root / "evaluation/metrics.json"),
            "freeze_sha256": metrics["freeze_sha256"]}


def parser():
    p = argparse.ArgumentParser(prog="rarecell", description=__doc__)
    subs = p.add_subparsers(dest="cmd", required=True)
    def cmd(name, cfg=True):
        sub = subs.add_parser(name)
        if cfg:
            sub.add_argument("--model", help="model for both discovery and annotation in this invocation")
            sub.add_argument("--config", type=Path, help="YAML overrides for the bundled defaults")
        return sub
    ui = cmd("ui", False)
    ui.add_argument("--workspace", type=Path, default=Path("work/tasks"), help="persistent dashboard task directory")
    ui.add_argument("--port", type=int, default=8787, help="loopback HTTP port (0 chooses a free port)")
    model_cmd = cmd("model", False)
    model_cmd.add_argument("--set", dest="value", help="save the project model (empty string clears it)")
    cmd("doctor")
    prep = cmd("prepare")
    prep.add_argument("--dataset", type=Path, required=True)
    prep.add_argument("--out", type=Path, required=True)
    an = cmd("annotate")
    an.add_argument("--snapshot", type=Path, required=True)
    an.add_argument("--out", type=Path, required=True)
    d = cmd("discover")
    d.add_argument("--snapshot", type=Path, required=True)
    d.add_argument("--run", type=Path, required=True)
    d.add_argument("--cache", type=Path)
    d.add_argument("--annotations", type=Path)
    d.add_argument("--mode", choices=["A0", "A1", "A2", "A3"])
    d.add_argument("--policy", choices=["live", "fixed"], default="live")
    d.add_argument("--replication-id")
    d.add_argument("--pause-after", type=int, help="checkpoint after this completed turn")
    resume_cmd = cmd("resume", False)
    resume_cmd.add_argument("--run", type=Path, required=True)
    resume_cmd.add_argument("--retry-unknown", action="store_true", help="allow one extra attempt for an incomplete external call")
    e = cmd("evaluate", False)
    e.add_argument("--run", type=Path, required=True)
    e.add_argument("--gt", type=Path, required=True)
    e.add_argument("--out", type=Path, required=True)
    dm = cmd("demo")
    dm.add_argument("--mode", choices=["core", "integration"], default="core")
    dm.add_argument("--out", type=Path, required=True)
    dm.add_argument("--seed", type=int, default=1009)
    ab = cmd("ablate")
    ab.add_argument("--snapshot", type=Path, required=True)
    ab.add_argument("--out", type=Path, required=True)
    ab.add_argument("--annotations", type=Path)
    ab.add_argument("--replication-prefix", help="assign each arm a planned ID PREFIX-ARM")
    ab.add_argument("--arms", nargs="+", choices=["A0", "A1", "A2", "A3", "fixed"], default=["A0", "A1", "A2", "A3", "fixed"])
    return p


def dispatch(a, default_model=""):
    cmd = a.cmd
    if cmd == "model":
        if a.value is not None:
            save_project_model(Path.cwd(), a.value)
        return {"model": project_model(Path.cwd())}
    if cmd == "ui":
        from .dashboard import serve
        return serve(a.workspace, a.port)
    if cmd == "resume":
        return resume(a.run, a.retry_unknown)
    if cmd == "evaluate":
        from .evaluation import evaluate
        metrics = evaluate(a.run, a.gt, a.out)
        return {"output": str(a.out), "freeze_sha256": metrics["freeze_sha256"], "submitted_partition": metrics["submitted_partition"]}
    cfg = config(a.config, model_overrides(a.model or default_model))
    if cmd == "doctor":
        return doctor(cfg)
    if cmd == "prepare":
        prepare(a.dataset, a.out, cfg)
        m = verify_snapshot(a.out)
        return {k: m.get(k) for k in ("snapshot_hash", "status", "views", "stages")}
    if cmd == "annotate":
        return annotation.annotate(a.snapshot, a.out, cfg)
    if cmd == "discover":
        return discover(a.snapshot, a.run, a.cache or a.run.parent / "query_cache", cfg, mode=a.mode,
                        annotations=a.annotations, policy=a.policy, replication_id=a.replication_id, pause_after=a.pause_after)
    if cmd == "demo":
        return demo(a.out, cfg, a.mode, a.seed)
    if cmd == "ablate":
        arms = list(dict.fromkeys(a.arms))
        if any((a.out / arm).exists() for arm in arms):
            raise ValueError("Use a fresh ablation output or resume each existing arm")
        sidecar = a.annotations or a.out / "annotation"
        if any(arm in {"A1", "A3", "fixed"} for arm in arms) and not a.annotations:
            annotation.annotate(a.snapshot, sidecar, cfg)
        reports = {}
        for arm in arms:
            summary = discover(a.snapshot, a.out / arm, a.out / "query_cache", cfg,
                               mode="A3" if arm == "fixed" else arm, policy="fixed" if arm == "fixed" else "live",
                               annotations=sidecar, replication_id=(a.replication_prefix + "-" + arm) if a.replication_prefix else None)
            reports[arm] = {"status": summary["status"], "detail": summary["detail"],
                            "snapshot_hash": summary["snapshot_hash"]}
        write_json(a.out / "ablation_summary.json", reports)
        return reports


def main(argv=None):
    try:
        args = parser().parse_args(argv)
        project = checkout(Path.cwd())
        if project:
            load_keys(project / ".api_keys")
        default = project_model(Path.cwd()) if hasattr(args, "config") else ""
        answer = dispatch(args, default)
        if answer is None:
            return 0
        print(json.dumps(answer, indent=2, ensure_ascii=False, default=str))
        if isinstance(answer, dict):
            result = answer.get("summary", answer)
            if result.get("status") in ("blocked", "partial") and result.get("detail") != "paused":
                return 2
            if "cmd" not in result and any(isinstance(v, dict) and v.get("detail") and "snapshot_hash" in v
                                            and v.get("status") == "partial" for v in result.values()):
                return 2
        return 0
    except (OSError, ValueError, KeyError) as exc:
        print("%s: %s" % (type(exc).__name__, exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
