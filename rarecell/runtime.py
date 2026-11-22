"""Project defaults and process startup for the CLI, dashboard and batch jobs."""
import os
from pathlib import Path
import re
import shlex
import socket
import sys


def checkout(directory):
    """Only read project credentials from an identifiable RareCell checkout."""
    path = Path(directory).resolve()
    for parent in (path, *path.parents):
        if (parent / "pyproject.toml").is_file() and (parent / "rarecell/r/annotation.R").is_file():
            return parent
    return None


# The credentials file may also be sourced by an interactive shell. Its
# literal PATH (e.g. "$PATH:...") cannot replace the initialized module paths.
_RUNTIME_ENV = {"PATH", "LD_LIBRARY_PATH", "LD_PRELOAD", "R_HOME", "RHOME",
                "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "RARECELL_MODULE_RSCRIPT"}


def rscript_command(settings):
    configured = settings.get("rscript", "Rscript")
    if configured == "Rscript":
        return os.environ.get("RARECELL_MODULE_RSCRIPT") or configured
    return configured


def load_keys(path):
    if not path.exists():
        return
    values = {}
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = re.fullmatch(r"(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)=(.*)", line)
        if match:
            if match[1] in _RUNTIME_ENV:
                continue
            try:
                # Parse only the assignment value. Quotes and inline comments are
                # supported, without evaluating shell commands from a key file.
                lexer = shlex.shlex(match[2], posix=True)
                lexer.whitespace_split = True
                lexer.commenters = ""
                value = lexer.get_token()
                rest = lexer.instream.read().lstrip()
            except ValueError:
                pass
            else:
                if (not rest or rest.startswith("#")) and "\x00" not in (value or ""):
                    values[match[1]] = value or ""
                    continue
        raise ValueError(f"{path.name} 第 {number} 行需要 NAME=value 或 export NAME=value")
    os.environ.update(values)


def validate_model(value):
    if not isinstance(value, str) or len(value) > 256 or any(c.isspace() or not c.isprintable() for c in value):
        raise ValueError("模型 ID 最多 256 个字符，不能包含空白或控制字符")
    return value


def model_file(project):
    return (checkout(project) or Path(project)) / ".rarecell-model"


def project_model(project):
    try:
        value = model_file(project).read_text(encoding="utf-8").removesuffix("\n")
    except FileNotFoundError:
        return ""
    return validate_model(value)


def save_project_model(project, value):
    from .config import atomic
    atomic(model_file(project), validate_model(value) + "\n")


def model_overrides(model):
    if not model:
        return {}
    return {section: {"model": validate_model(model)} for section in ("provider", "annotation")}


def main():
    root = checkout(Path.cwd())
    module = os.environ.get("RARECELL_R_MODULE")
    if module is None and socket.gethostname().split(".")[0].startswith("dcc-"):
        module = "R/4.3.1-rhel8"
    if root and module:
        # Start a new interpreter after loading the module so the dynamic
        # linker, Rscript and detached workers all inherit the same environment.
        script = Path(__file__).parent / "resources" / "load-r.sh"
        os.execv("/bin/bash", ["/bin/bash", str(script), str(root), module, sys.executable, *sys.argv[1:]])
    from .cli import main as run
    return run()
