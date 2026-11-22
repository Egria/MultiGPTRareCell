"""Project runtime settings across the GUI queue and batch entry point."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from rarecell import job_queue as jobs
from rarecell.config import config
from rarecell.runtime import load_keys, project_model, save_project_model


def test_key_file_is_literal_and_invalid_lines_are_private(tmp_path, monkeypatch):
    path = tmp_path / ".api_keys"
    monkeypatch.setenv("TEST_PROVIDER_KEY", "inherited")
    path.write_text('export TEST_PROVIDER_KEY="space # quote=ok" # trailing comment\n')
    load_keys(path)
    assert os.environ["TEST_PROVIDER_KEY"] == "space # quote=ok"
    path.write_text("TEST_PROVIDER_KEY=first\nnot valid sensitive-value\n")
    with pytest.raises(ValueError, match="第 2 行") as failure:
        load_keys(path)
    assert "sensitive-value" not in str(failure.value)
    assert os.environ["TEST_PROVIDER_KEY"] == "space # quote=ok"


def test_frozen_gui_models_and_clear(tmp_path):
    root = tmp_path / "queue"
    jobs.initialize(root)
    yaml_file = tmp_path / "a.yaml"
    yaml_file.write_text("provider: {model: yaml-disc}\nannotation: {model: yaml-ann}\n")
    form = {"kind": "pipeline", "dataset": str(yaml_file), "config": str(yaml_file)}
    save_project_model(tmp_path, "shared")
    first = jobs.submit(root, tmp_path, dict(form, annotation_model="custom"))
    save_project_model(tmp_path, "new")
    frozen = config(root / first / "settings.yaml")
    assert [frozen[s]["model"] for s in ("provider", "annotation")] == ["shared", "custom"]
    save_project_model(tmp_path, "")
    next_id = jobs.submit(root, tmp_path, form)
    resolved = config(root / next_id / "settings.yaml")
    assert [resolved[s]["model"] for s in ("provider", "annotation")] == ["yaml-disc", "yaml-ann"]


def test_checkout_batch_entry_reads_key_and_global_model(tmp_path):
    repo = Path(__file__).resolve().parents[1]
    # Isolate the credentials and project settings while executing the real code.
    for name in ("pyproject.toml", "rarecell"):
        (tmp_path / name).symlink_to(repo / name)
    (tmp_path / ".venv").symlink_to(Path(sys.executable).parent.parent, target_is_directory=True)
    secret = "local-offline-test-only"
    (tmp_path / ".api_keys").write_text("OPENAI_API_KEY=" + secret)
    save_project_model(tmp_path, "shared")
    env = dict(os.environ, SLURM_SUBMIT_DIR=str(tmp_path), RARECELL_R_MODULE="")
    env.pop("OPENAI_API_KEY", None)
    for extra, expected in (([], "shared"), (["--model", "one-job"], "one-job")):
        proc = subprocess.run(["/bin/bash", str(repo / "scripts/rarecell.sbatch"), "doctor", *extra],
                              env=env, cwd="/", capture_output=True, text=True, check=True)
        discovery = json.loads(proc.stdout)["discovery_provider"]
        assert discovery["status"] == "configured" and discovery["model"] == expected
        assert secret not in proc.stdout + proc.stderr
    assert project_model(tmp_path) == "shared"
