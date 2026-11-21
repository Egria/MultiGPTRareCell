"""Durable, process-safe local experiment state. Only trusted controllers write here."""

import contextlib
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile
import time
import uuid


class IntegrityError(ValueError):
    """An immutable object changed, or a declared dependency is invalid."""


class BudgetExhausted(ValueError):
    """This run cannot admit another distinct scientific query."""


class QueryBlocked(ValueError):
    """A registered local operation failed permanently or reached its retry limit."""


def _dump(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(value):
    return hashlib.sha256(_dump(value).encode("utf-8")).hexdigest()


def key_for(canonical):
    """The caller resolves member sets, aliases, versioned parameters and scope first."""
    if not isinstance(canonical, dict) or not canonical:
        raise ValueError("canonical scientific request must be a nonempty object")
    return digest({"schema": "rarecell-local-query-v1", "request": canonical})


def _sync_dir(path):
    fd = os.open(str(path), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _install_json(path, value):
    """Create once, even when multiple independent writers install the same JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    data = _dump(value).encode("utf-8")
    fd, name = tempfile.mkstemp(dir=str(path.parent), prefix=".pending-")
    try:
        with os.fdopen(fd, "wb") as tmp:
            tmp.write(data)
            tmp.flush()
            os.fsync(tmp.fileno())
        try:
            os.link(name, path)
        except FileExistsError:
            if path.read_bytes() != data:
                raise IntegrityError("immutable artifact has different content: " + path.name)
        else:
            _sync_dir(path.parent)
    finally:
        os.unlink(name)


@contextlib.contextmanager
def _lock(path, *, nonblocking=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("ab") as stream:
        flags = fcntl.LOCK_EX | (fcntl.LOCK_NB if nonblocking else 0)
        try:
            fcntl.flock(stream, flags)
        except BlockingIOError:
            yield False
            return
        try:
            yield True
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


_SCHEMA = """
CREATE TABLE IF NOT EXISTS record (
    kind TEXT NOT NULL, id TEXT NOT NULL, body TEXT NOT NULL, PRIMARY KEY (kind,id)
);
CREATE TABLE IF NOT EXISTS evidence_content (
    content_hash TEXT PRIMARY KEY, evidence_id TEXT UNIQUE NOT NULL
);
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, body TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS event (
    seq INTEGER PRIMARY KEY AUTOINCREMENT, type TEXT NOT NULL, at REAL NOT NULL, payload TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS query (
    key TEXT PRIMARY KEY, canonical TEXT NOT NULL, status TEXT NOT NULL,
    result_ref TEXT, error TEXT, global_hit INTEGER NOT NULL DEFAULT 0,
    reuses INTEGER NOT NULL DEFAULT 0, recovered INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE IF NOT EXISTS attempt (
    id TEXT PRIMARY KEY, key TEXT NOT NULL REFERENCES query(key), status TEXT NOT NULL,
    started REAL NOT NULL, wall REAL, cpu REAL, error TEXT
);
CREATE INDEX IF NOT EXISTS attempt_by_key ON attempt(key);
"""


class State:
    def __init__(self, run, cache):
        self.run, self.cache = Path(run), Path(cache)
        self.run.mkdir(parents=True, exist_ok=True)
        self.cache.mkdir(parents=True, exist_ok=True)
        (self.run / "artifacts").mkdir(exist_ok=True)
        self.db = sqlite3.connect(str(self.run / "state.sqlite"), isolation_level=None, timeout=60)
        self.db.row_factory = sqlite3.Row
        # Rollback journals work on more shared HPC filesystems than WAL.
        self.db.execute("PRAGMA foreign_keys=ON")
        self.db.execute("PRAGMA synchronous=FULL")
        self.db.executescript(_SCHEMA)
        self._event_dirty = False
        self._mirror_events(full=True)
        self._recover()

    def close(self):
        self.db.close()

    @contextlib.contextmanager
    def transaction(self):
        """Atomic controller edits; nested uses are isolated by savepoints."""
        nested = self.db.in_transaction
        savepoint = "s" + uuid.uuid4().hex
        self.db.execute("SAVEPOINT " + savepoint if nested else "BEGIN IMMEDIATE")
        try:
            yield self
            self.db.execute("RELEASE " + savepoint if nested else "COMMIT")
        except BaseException:
            if self.db.in_transaction:
                if nested:
                    self.db.execute("ROLLBACK TO " + savepoint)
                    self.db.execute("RELEASE " + savepoint)
                else:
                    self.db.rollback()
            raise
        if not nested and self._event_dirty:
            self._mirror_events()
            self._event_dirty = False

    def _put(self, kind, identifier, obj):
        if not isinstance(identifier, str) or not identifier or not isinstance(obj, dict):
            raise ValueError("immutable objects need a nonempty string ID and JSON object")
        body = _dump(obj)
        self.db.execute("INSERT OR IGNORE INTO record VALUES (?,?,?)", (kind, identifier, body))
        stored = self.db.execute("SELECT body FROM record WHERE kind=? AND id=?", (kind, identifier)).fetchone()
        if stored["body"] != body:
            raise IntegrityError("immutable version differs for %s/%s" % (kind, identifier))
        return identifier

    def put(self, kind, identifier, obj):
        if kind in {"Evidence", "CanonicalQuery", "ExecutionAttempt"}:
            raise ValueError("use the dedicated writer for " + kind)
        with self.transaction():
            return self._put(kind, identifier, obj)

    def get(self, kind, identifier):
        if kind == "CanonicalQuery":
            row = self.db.execute("SELECT * FROM query WHERE key=?", (identifier,)).fetchone()
            if row:
                return self._query_dict(row)
        elif kind == "ExecutionAttempt":
            row = self.db.execute("SELECT * FROM attempt WHERE id=?", (identifier,)).fetchone()
            if row:
                return self._attempt_dict(row)
        else:
            row = self.db.execute("SELECT body FROM record WHERE kind=? AND id=?", (kind, identifier)).fetchone()
            if row:
                return json.loads(row["body"])
        raise KeyError((kind, identifier))

    @staticmethod
    def _query_dict(row):
        return {"canonical_query_key": row["key"], "canonical": json.loads(row["canonical"]),
                "status": row["status"], "result_ref": row["result_ref"], "error": row["error"]}

    @staticmethod
    def _attempt_dict(row):
        return {"execution_attempt_id": row["id"], "canonical_query_key": row["key"],
                "status": row["status"], "started_unix": row["started"], "wall_seconds": row["wall"],
                "cpu_seconds": row["cpu"], "error": row["error"]}

    def all(self, kind):
        if kind == "CanonicalQuery":
            return [self._query_dict(r) for r in self.db.execute("SELECT * FROM query ORDER BY rowid")]
        if kind == "ExecutionAttempt":
            return [self._attempt_dict(r) for r in self.db.execute("SELECT * FROM attempt ORDER BY rowid")]
        rows = self.db.execute("SELECT body FROM record WHERE kind=? ORDER BY rowid", (kind,))
        return [json.loads(r["body"]) for r in rows]

    def get_meta(self, key, default=None):
        row = self.db.execute("SELECT body FROM meta WHERE key=?", (key,)).fetchone()
        return default if row is None else json.loads(row["body"])

    def set_meta(self, key, value):
        self.db.execute("INSERT INTO meta VALUES (?,?) ON CONFLICT(key) DO UPDATE SET body=excluded.body",
                        (key, _dump(value)))

    def event(self, type_, payload):
        with self.transaction():
            cur = self.db.execute("INSERT INTO event(type,at,payload) VALUES (?,?,?)",
                                  (type_, time.time(), _dump(payload)))
            self._event_dirty = True
            return cur.lastrowid

    def _mirror_events(self, full=False):
        """Append committed SQLite events; repair a torn or missing journal on startup."""
        path = self.run / "events.jsonl"
        with _lock(self.run / ".events.lock"), path.open("a+b") as stream:
            seq = 0
            if not full and stream.tell():
                # One event per line. Seek backwards by a block for large payloads.
                end = stream.tell()
                stream.seek(end - 1)
                if stream.read(1) != b"\n":
                    full = True
                else:
                    pos = end - 1
                    while pos > 0:
                        start = max(0, pos - 4096)
                        stream.seek(start)
                        previous = stream.read(pos - start).rfind(b"\n")
                        if previous >= 0:
                            pos = start + previous + 1
                            break
                        pos = start
                    stream.seek(pos)
                    try:
                        seq = json.loads(stream.readline())["seq"]
                    except (ValueError, KeyError, TypeError, UnicodeError):
                        full = True
            if full:
                stream.seek(0)
                # A consistent database read plus the mirror lock orders writers.
                for row in self.db.execute("SELECT * FROM event ORDER BY seq"):
                    pos = stream.tell()
                    if stream.readline() != self._event_line(row):
                        stream.seek(pos)
                        break
                    seq = row["seq"]
                stream.truncate(stream.tell())
            stream.seek(0, os.SEEK_END)
            for row in self.db.execute("SELECT * FROM event WHERE seq>? ORDER BY seq", (seq,)):
                stream.write(self._event_line(row))
            stream.flush()
            os.fsync(stream.fileno())

    @staticmethod
    def _event_line(row):
        return (_dump({"seq": row["seq"], "type": row["type"], "at": row["at"],
                       "payload": json.loads(row["payload"])}) + "\n").encode("utf-8")

    def artifact(self, payload):
        sha = digest(payload)
        ref = "artifacts/" + sha + ".json"
        _install_json(self.run / ref, payload)
        self.put("Artifact", sha, {"id": sha, "ref": ref, "sha256": sha})
        return ref

    def read_artifact(self, ref):
        match = re.fullmatch(r"artifacts/([0-9a-f]{64})\.json", ref)
        if not match:
            raise IntegrityError("unregistered artifact reference format")
        sha = match.group(1)
        self.get("Artifact", sha)
        path = self.run / ref
        if path.is_symlink() or path.resolve().parent != (self.run / "artifacts").resolve():
            raise IntegrityError("invalid artifact reference")
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != sha:
            raise IntegrityError("local artifact checksum mismatch")
        return json.loads(data)

    def evidence(self, obj):
        """Register an immutable node. Only already registered parents may be linked."""
        required = {"kind", "parent_ids", "source_family", "cell_set_hash", "data_version",
                    "method_version", "config_version", "prompt_version", "selection_history",
                    "artifact_refs", "purpose", "object", "summary"}
        missing = required.difference(obj)
        if missing:
            raise IntegrityError("missing evidence fields: " + ", ".join(sorted(missing)))
        entry = dict(obj)
        explicit = entry.pop("evidence_id", None)
        if entry["kind"] not in {"observation", "interpretation", "hypothesis", "criticism", "validation"}:
            raise IntegrityError("unknown evidence kind")
        if entry["purpose"] == "independent_validation" or entry["kind"] == "validation":
            raise IntegrityError("no independent validation protocol registered in discovery V1")
        if entry["purpose"] not in {"exploration", "related_check", "selection_dependent"}:
            raise IntegrityError("unknown evidence purpose")
        parents = entry["parent_ids"]
        if not isinstance(parents, list) or not all(isinstance(p, str) and p for p in parents):
            raise IntegrityError("parent_ids must be a list of evidence IDs")
        entry["parent_ids"] = sorted(set(parents))
        if entry["kind"] != "observation" and not parents:
            raise IntegrityError("a derived claim requires parent evidence")
        if not isinstance(entry["artifact_refs"], list):
            raise IntegrityError("artifact_refs must be a list")
        # References are run-local content-addressed artifacts installed by the controller.
        for ref in entry["artifact_refs"]:
            if not isinstance(ref, str) or not re.fullmatch(r"artifacts/[0-9a-f]{64}\.json", ref):
                raise IntegrityError("invalid evidence artifact reference")
            try:
                self.get("Artifact", Path(ref).stem)
            except KeyError as exc:
                raise IntegrityError("unregistered evidence artifact") from exc
        content = digest(entry)
        eid = explicit if explicit is not None else "E" + content
        if eid in entry["parent_ids"]:
            raise IntegrityError("cyclic evidence dependency")
        with self.transaction():
            for parent in entry["parent_ids"]:
                try:
                    self.get("Evidence", parent)
                except KeyError as exc:
                    raise IntegrityError("missing parent evidence " + parent) from exc
            equivalent = self.db.execute("SELECT evidence_id FROM evidence_content WHERE content_hash=?",
                                         (content,)).fetchone()
            if equivalent:
                if explicit is not None and explicit != equivalent[0]:
                    raise IntegrityError("equivalent evidence already has ID " + equivalent[0])
                return equivalent[0]
            self._put("Evidence", eid, dict(entry, evidence_id=eid))
            self.db.execute("INSERT INTO evidence_content VALUES (?,?)", (content, eid))
        return eid

    def result_ref(self, key):
        record = self.get("CanonicalQuery", key)
        if record["status"] != "success":
            raise QueryBlocked("no committed result for " + key)
        return record["result_ref"]

    def _entry_path(self, key):
        return self.cache / "results" / key[:2] / (key + ".json")

    @contextlib.contextmanager
    def _query_lock(self, key, **options):
        # The run lock also protects a run resumed with a relocated shared cache.
        with _lock(self.run / ".query_locks" / (key + ".lock"), **options) as local:
            if not local:
                yield False
                return
            with _lock(self.cache / "locks" / (key + ".lock"), **options) as shared:
                yield shared

    def _admit(self, key, canonical, budget):
        if not isinstance(budget, int) or isinstance(budget, bool) or budget < 0:
            raise ValueError("logical query budget must be a nonnegative integer")
        body = _dump(canonical)
        with self.transaction():
            row = self.db.execute("SELECT canonical FROM query WHERE key=?", (key,)).fetchone()
            if row is not None:
                if row["canonical"] != body:
                    raise IntegrityError("canonical request hash collision")
                return
            n = self.db.execute("SELECT count(*) FROM query").fetchone()[0]
            if n >= budget:
                raise BudgetExhausted("max_new_queries reached")
            self.db.execute("INSERT INTO query(key,canonical,status) VALUES (?,?,'pending')", (key, body))
            self.event("query_admitted", {"canonical_query_key": key, "new_queries": n + 1})

    def _load_cache(self, key):
        path = self._entry_path(key)
        if not path.exists():
            return None
        try:
            entry = json.loads(path.read_bytes())
            if (entry["canonical_query_key"] != key or entry["version"] != 1 or
                    not isinstance(entry["result"], dict) or digest(entry["result"]) != entry["result_sha256"]):
                raise IntegrityError("cached result identity or checksum mismatch")
            return entry
        except (KeyError, TypeError, ValueError) as exc:
            raise IntegrityError("invalid shared query cache: " + key) from exc

    def _finish(self, key, entry, *, recovering=False):
        with self.transaction():
            own = self.db.execute("SELECT id FROM attempt WHERE id=? AND key=?",
                                  (entry["execution_attempt_id"], key)).fetchone()
            ref = self.artifact(entry["result"])
            if own is not None:
                self.db.execute("UPDATE attempt SET status='success',wall=?,cpu=?,error=NULL WHERE id=?",
                                (entry["wall_seconds"], entry["cpu_seconds"], own["id"]))
            self.db.execute("UPDATE attempt SET status='interrupted' WHERE key=? AND status='running'", (key,))
            self.db.execute("UPDATE query SET status='success',result_ref=?,error=NULL,global_hit=?,recovered=? WHERE key=?",
                            (ref, int(own is None), int(recovering), key))
            self.event("query_succeeded", {"canonical_query_key": key, "artifact_ref": ref,
                                           "global_cache_hit": own is None, "recovered": recovering})
        return entry["result"]

    def _recover(self):
        rows = self.db.execute("SELECT key FROM query WHERE status IN ('pending','running','retryable')").fetchall()
        for row in rows:
            key = row["key"]
            with self._query_lock(key, nonblocking=True) as locked:
                if not locked:
                    continue
                current = self.get("CanonicalQuery", key)
                if current["status"] == "success":
                    continue
                entry = self._load_cache(key)
                if entry is not None:
                    self._finish(key, entry, recovering=True)
                elif current["status"] == "running":
                    self._interrupt(key)

    def _interrupt(self, key):
        with self.transaction():
            self.db.execute("UPDATE attempt SET status='interrupted' WHERE key=? AND status='running'", (key,))
            self.db.execute("UPDATE query SET status='retryable',error='interrupted before publishing result' WHERE key=?",
                            (key,))
            self.event("query_interrupted", {"canonical_query_key": key})

    def query(self, canonical, budget, compute, max_attempts=2):
        """Admit, resolve cache, or run a deterministic local callable. Never call external APIs here."""
        if self.db.in_transaction:
            raise ValueError("run scientific query outside controller write transaction")
        if not isinstance(max_attempts, int) or isinstance(max_attempts, bool) or max_attempts < 1:
            raise ValueError("max_attempts must be a positive integer")
        key = key_for(canonical)
        self._admit(key, canonical, budget)
        with self._query_lock(key):
            current = self.get("CanonicalQuery", key)
            if current["status"] == "success":
                result = self.read_artifact(current["result_ref"])
                self.db.execute("UPDATE query SET reuses=reuses+1 WHERE key=?", (key,))
                return result, key
            entry = self._load_cache(key)
            if entry is not None:
                return self._finish(key, entry, recovering=current["status"] == "running"), key
            if current["status"] == "blocked":
                raise QueryBlocked(current["error"])
            if current["status"] == "running":
                self._interrupt(key)
            for _ in range(max_attempts):
                n = self.db.execute("SELECT count(*) FROM attempt WHERE key=?", (key,)).fetchone()[0]
                if n >= max_attempts:
                    break
                aid = "X" + uuid.uuid4().hex
                with self.transaction():
                    self.db.execute("INSERT INTO attempt(id,key,status,started) VALUES (?,?,'running',?)",
                                    (aid, key, time.time()))
                    self.db.execute("UPDATE query SET status='running',error=NULL WHERE key=?", (key,))
                    self.event("compute_started", {"execution_attempt_id": aid, "canonical_query_key": key})
                wall, cpu = time.monotonic(), time.process_time()
                try:
                    result = compute()
                    if not isinstance(result, dict):
                        raise TypeError("local tool must return a JSON object")
                    sha = digest(result)
                except Exception as exc:
                    retryable = isinstance(exc, (TimeoutError, OSError))
                    status = "retryable" if retryable else "blocked"
                    error = type(exc).__name__ + ": " + str(exc)
                    with self.transaction():
                        self.db.execute("UPDATE attempt SET status=?,error=?,wall=?,cpu=? WHERE id=?",
                                        (status, error, time.monotonic() - wall, time.process_time() - cpu, aid))
                        self.db.execute("UPDATE query SET status=?,error=? WHERE key=?", (status, error, key))
                        self.event("compute_failed", {"execution_attempt_id": aid, "canonical_query_key": key,
                                                      "status": status, "error": error})
                    if retryable:
                        continue
                    raise QueryBlocked(error) from exc
                entry = {"version": 1, "canonical_query_key": key, "execution_attempt_id": aid,
                         "result_sha256": sha, "result": result, "cpu_seconds": time.process_time() - cpu,
                         "wall_seconds": time.monotonic() - wall}
                _install_json(self._entry_path(key), entry)
                return self._finish(key, entry), key
            raise QueryBlocked("local attempt limit reached: " + str(self.get("CanonicalQuery", key)["error"]))

    def usage(self):
        q = self.db.execute("""SELECT count(*) AS new_queries, coalesce(sum(global_hit),0) AS global_cache_hits,
                           coalesce(sum(reuses),0) AS run_reuses, coalesce(sum(recovered),0) AS recovered_results
                           FROM query""").fetchone()
        a = self.db.execute("""SELECT count(*) AS compute_executions, count(*)-count(DISTINCT key) AS technical_retries,
                           coalesce(sum(cpu),0.0) AS cpu_seconds, coalesce(sum(wall),0.0) AS wall_seconds,
                           count(CASE WHEN status IN ('running','interrupted') THEN 1 END) AS unmeasured_attempts
                           FROM attempt""").fetchone()
        return dict(dict(q), **dict(a))
