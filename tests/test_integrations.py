"""Wire protocol fixtures; these tests do not measure a discovery policy."""
import json
import os
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
import threading

import pytest

from rarecell import annotation, providers


@pytest.fixture
def endpoint():
    calls = []
    answer = {}
    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            calls.append((self.path, json.loads(self.rfile.read(int(self.headers["Content-Length"])))))
            self.send_response(200)
            self.end_headers()
            self.wfile.write(json.dumps(answer).encode())

        def log_message(self, *_args):
            return

    httpd = HTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield "http://127.0.0.1:{}".format(httpd.server_port), calls, answer
    httpd.shutdown()
    thread.join()
    httpd.server_close()


def test_provider_envelope_and_interruption(endpoint):
    base, calls, answer = endpoint
    cfg = {"base_url": base, "model": "test-envelope", "api_key_env": None}
    answer.update(choices=[{"message": {"content": '{"action":"list_candidates"}'}, "finish_reason": "stop"}],
                  usage={"prompt_tokens": 30, "completion_tokens": 8})
    message, usage = providers.complete({"contract": "C01"}, cfg)
    assert json.loads(message)["action"] == "list_candidates"
    assert usage["input_tokens"] == 30 and "total_tokens" not in usage
    assert calls[0][0] == "/chat/completions"
    assert set(calls[0][1]) == {"model", "messages"}
    answer["choices"][0]["finish_reason"] = "length"
    with pytest.raises(providers.ProviderError) as err:
        providers.complete({}, cfg)
    assert err.value.status == "incomplete" and len(calls) == 2


def test_responses_dialect(endpoint):
    base, calls, answer = endpoint
    answer.update(status="completed", output=[{"type": "message", "role": "assistant",
                  "content": [{"type": "output_text", "text": "hello"}]}])
    cfg = {"base_url": base, "model": "test", "api_style": "responses", "endpoint": "custom/path",
           "api_key_env": None, "options": {"max_output_tokens": 12}}
    assert providers.complete({"current": 2}, cfg)[0] == "hello"
    assert calls[0][0] == "/custom/path" and calls[0][1]["max_output_tokens"] == 12
    answer["status"] = "incomplete"
    with pytest.raises(providers.ProviderError):
        providers.complete({}, cfg)
    assert "model" in providers.diagnose({"base_url": base})["detail"]


def _marker_snapshot(path):
    (path / "seurat").mkdir()
    ids = ["H1", "H2"]
    candidates = [{"id": k, "members": [k + "_a", k + "_b", k + "_c"] if k == "H1" else [k],
                   "sources": [{"view": "V_HVG", "kind": "partition", "resolution": 1.0,
                                "partition_id": "pi"}]} for k in ids]
    (path / "initial.jsonl").write_text("".join(json.dumps(c) + "\n" for c in candidates))
    (path / "seurat" / "markers.tsv").write_text(
        "resolution\tgroup_id\tgene_id\tsymbol\tavg_log2FC\tp_val_adj\n"
        "1\tH1\tg1\tCD3D\t2\t0.02\n"
        "1\tH1\tg2\tCD3E\t1\t0.01\n"
        "1\tH1\tg3\tA\t-3\t0.001\n"
    )


def test_annotation_contents_and_separate_singleton(tmp_path):
    _marker_snapshot(tmp_path)
    manifest = {"snapshot_hash": "snap", "data_version": "v1"}
    settings = annotation._resolve({"model": "m"})
    first, skipped = annotation._jobs(tmp_path, manifest, settings)
    assert len(first) == 1 and first[0]["groups"][0]["genes"] == ["CD3E", "CD3D"]
    assert len(skipped) == 1 and skipped[0]["reason"] == "singleton"
    key = annotation._hash(first[0])
    changed, _ = annotation._jobs(tmp_path, manifest, dict(settings, model="m2"))
    assert annotation._hash(changed[0]) != key
    path = tmp_path / "seurat" / "markers.tsv"
    path.write_text(path.read_text().replace("CD3D", "CD8A"))
    changed, _ = annotation._jobs(tmp_path, manifest, settings)
    assert annotation._hash(changed[0]) != key


def test_annotation_unknown_dispatch_and_receipt_recovery(tmp_path):
    _marker_snapshot(tmp_path)
    manifest = {"snapshot_hash": "snap", "data_version": "v1"}
    job = annotation._jobs(tmp_path, manifest, annotation._resolve({"model": "m"}))[0][0]
    key = annotation._hash(job)
    dst = tmp_path / "batches" / key
    dst.mkdir(parents=True)
    (dst / "dispatch.json").write_text('{"attempts":1}')
    with pytest.raises(RuntimeError, match="unknown"):
        annotation._batch(tmp_path, job, {})
    response = {"commit": annotation.PIN, "valid": True, "raw_response": "unknown",
                "rows": [{"candidate_id": "H1", "label": "unknown", "ontology_id": None}]}
    (dst / "response.json").write_text(json.dumps(response))
    rows, hit, refs = annotation._batch(tmp_path, job, {"rscript": "must-not-run"})
    assert hit and rows[0]["kind"] == "interpretation" and len(refs) == 2
    response["rows"][0]["candidate_id"] = "H2"
    (dst / "response.json").write_text(json.dumps(response))
    with pytest.raises(ValueError, match="alignment"):
        annotation._batch(tmp_path, job, {})


def test_doctor_reports_actual_external_requirements(monkeypatch):
    monkeypatch.delenv("NO_SUCH_TEST_KEY", raising=False)
    cfg = {"model": "m", "api_key_env": "NO_SUCH_TEST_KEY", "rscript": "/absent/Rscript"}
    state = annotation.diagnose(cfg)
    assert state["status"] == "blocked"
    assert "NO_SUCH_TEST_KEY" in state["detail"] and "Rscript" in state["detail"]
    assert annotation.diagnose({"annotation": False})["status"] == "not_run"


def test_real_ontoanno_interface_if_installed():
    import shutil
    import subprocess
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("requires external R/OntoAnno environment")
    present = subprocess.run([rscript, "--vanilla", "-e",
                              'quit(status=if(requireNamespace("OntoAnno",quietly=TRUE)) 0 else 1)'],
                             capture_output=True, timeout=90)
    if present.returncode:
        pytest.skip("OntoAnno or its R imports unavailable")
    probe = subprocess.run([rscript, "--vanilla", str(annotation.SCRIPT), "--probe"],
                           capture_output=True, text=True, check=True, timeout=90)
    assert annotation.PIN in probe.stdout
    cmd = ('m <- OntoAnno::map_celltypes_to_cl(c("T cell","not_a_real_CL_term"),verbose=FALSE); '
           'stopifnot(m$clid[[1]] == "CL:0000084", is.na(m$clid[[2]]))')
    subprocess.run([rscript, "--vanilla", "-e", cmd], check=True, timeout=90)


def test_gate_copied_markers():
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "reply.json"
        path.write_text(json.dumps({"valid": True, "commit": annotation.PIN,
                                    "rows": [{"candidate_id": "H1", "label": "M1, M2"}],
                                    "raw_response": "M1, M2"}))
        job = {"groups": [{"candidate_id": "H1", "genes": ["M1", "M2"]}]}
        with pytest.raises(ValueError, match="echoed marker"):
            annotation._validate_receipt(path, job)


def test_partial_annotation_exact_successful_subset(tmp_path):
    """A receipt fixture exercises partial reads without rerunning a provider."""
    import hashlib
    fixture = tmp_path / "snap"
    fixture.mkdir()
    _marker_snapshot(fixture)
    # Two independent eligible clusters, one succeeds, the other is recorded failed.
    initial = fixture / "initial.jsonl"
    groups = [json.loads(line) for line in initial.read_text().splitlines()]
    groups[1]["members"] = ["b1", "b2", "b3"]
    initial.write_text("".join(json.dumps(row) + "\n" for row in groups))
    with (fixture / "seurat" / "markers.tsv").open("a") as f:
        f.write("1\tH2\tg2\tCD3E\t1\t0.01\n")
    manifest = {"snapshot_hash": "fixture", "data_version": "data"}
    settings = annotation._resolve({"model": "test"})
    jobs, _ = annotation._jobs(fixture, manifest, settings)
    assert len(jobs) == 2
    key, missing = map(annotation._hash, jobs)
    out = tmp_path / "labels"
    dest = out / "batches" / key
    dest.mkdir(parents=True)
    (dest / "response.json").write_text(json.dumps(
        {"valid": True, "commit": annotation.PIN, "raw_response": "unknown",
         "rows": [{"candidate_id": "H1", "label": "unknown", "ontology_id": None}]}))
    rows, hit, refs = annotation._batch(out, jobs[0], {})
    assert hit
    encoded = json.dumps(rows[0]).encode() + b"\n"
    (out / "annotations.jsonl").write_bytes(encoded)
    state = {"status": "partial", "annotated": 1, "snapshot_hash": "fixture", "settings": settings,
             "annotation_sha256": hashlib.sha256(encoded).hexdigest(), "errors": [{"key": missing}],
             "source_commit": annotation.PIN, "files": refs}
    status_file = out / "annotation_status.json"
    status_file.write_text(json.dumps(state))
    with pytest.raises(ValueError, match="partial"):
        annotation.verify_annotation(fixture, out, _manifest=manifest)
    assert annotation.verify_annotation(fixture, out, allow_partial=True, _manifest=manifest) == rows
    state["errors"].append({"key": key})
    status_file.write_text(json.dumps(state))
    with pytest.raises(ValueError):
        annotation.verify_annotation(fixture, out, allow_partial=True, _manifest=manifest)
