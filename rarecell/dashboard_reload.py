"""Keep the loopback listener and login alive while replacing the GUI process."""
from pathlib import Path
import secrets
import socket
import subprocess
import sys
import tempfile
import time

PACKAGE = Path(__file__).resolve().parent


def source_stamp():
    paths = sorted([*PACKAGE.glob("*.py"), *(PACKAGE / "web").glob("*")])
    versions = []
    for path in paths:
        try:
            stat = path.stat()
        except FileNotFoundError:
            continue
        versions.append((path.name, stat.st_ino, stat.st_mtime_ns, stat.st_size))
    return versions


def start_backend(listener, root, project, token):
    control, child_end = socket.socketpair()
    proc = None
    try:
        # Read current source even for equal-size edits within a .pyc timestamp.
        with child_end, tempfile.TemporaryDirectory(prefix="rarecell-ui-imports-") as empty_cache:
            argv = [sys.executable, "-B", "-X", "pycache_prefix=" + empty_cache, "-m", "rarecell.dashboard",
                    str(root), str(project), str(listener.fileno()), str(child_end.fileno())]
            proc = subprocess.Popen(argv, pass_fds=(listener.fileno(), child_end.fileno()), start_new_session=True)
        control.settimeout(30)
        control.sendall(token.encode("ascii") + b"\n")
        if control.recv(1) != b"1":
            raise OSError("新版 GUI 未完成启动")
        return proc, control
    except BaseException:
        control.close()
        if proc is not None:
            proc.terminate()
            proc.wait()
        raise


def serve(root, port):
    root, project = Path(root).resolve(), Path.cwd()
    token = secrets.token_urlsafe(32)
    with socket.create_server(("127.0.0.1", port), backlog=64) as listener:
        observed = attempted = source_stamp()
        proc, control = start_backend(listener, root, project, token)
        print("RareCell 任务控制台（私人链接）", flush=True)
        print("http://127.0.0.1:%d/#%s" % (listener.getsockname()[1], token), flush=True)
        retired = []
        try:
            while True:
                time.sleep(1.5)
                retired = [p for p in retired if p.poll() is None]
                stamp = source_stamp()
                if stamp != observed:
                    observed = stamp
                    continue  # Wait for two equal scans, e.g. during git pull.
                if proc.poll() is not None:
                    raise RuntimeError("GUI 服务进程已退出")
                if stamp == attempted:
                    continue
                attempted = stamp
                try:
                    new_proc, new_control = start_backend(listener, root, project, token)
                except OSError as exc:
                    print("GUI 更新未能启动，原版本继续服务：", exc, file=sys.stderr, flush=True)
                    continue
                control.close()  # Old backend drains requests; detached task workers keep their lease.
                retired.append(proc)
                proc, control = new_proc, new_control
                print("GUI 已热更新。", flush=True)
        except KeyboardInterrupt:
            pass
        finally:
            control.close()
            proc.wait()
            for old in retired:
                old.wait()
