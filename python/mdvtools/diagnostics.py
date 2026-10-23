"""
Version and error-report formatting for the mdvtools CLI.

``format_version_info`` is what ``mdvtools --version`` prints. ``report_errors``
wraps any callable so that an unexpected exception is printed as a single
report (traceback, command line, versions, platform) suitable for pasting into
a bug report, and the process exits with status 1.
"""
import functools
import json
import os
import platform
import shlex
import sys
import traceback
from importlib import metadata
from typing import Any, Callable, Dict, Optional, ParamSpec, TypeVar

import click

from .build_info import get_build_info

ISSUE_TRACKER_URL = "https://github.com/Taylor-CCB-Group/MDV/issues"

# Libraries whose versions most often explain a conversion/serving failure.
KEY_DEPENDENCIES = (
    "anndata",
    "scanpy",
    "mudata",
    "spatialdata",
    "spatialdata-io",
    "zarr",
    "numpy",
    "pandas",
    "h5py",
    "click",
)

P = ParamSpec("P")
R = TypeVar("R")


def _package_version(name: str) -> Optional[str]:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def _install_kind() -> Optional[str]:
    """'editable' or 'installed', from the distribution's PEP 610 record."""
    try:
        dist = metadata.distribution("mdvtools")
    except metadata.PackageNotFoundError:
        return None
    direct_url = dist.read_text("direct_url.json")
    if direct_url:
        try:
            if json.loads(direct_url).get("dir_info", {}).get("editable"):
                return "editable"
        except ValueError:
            pass
    return "installed"


def get_version_info(include_dependencies: bool = False) -> Dict[str, Any]:
    """Collect mdvtools version, build, interpreter and platform details."""
    info: Dict[str, Any] = {
        "mdvtools": _package_version("mdvtools") or "unknown",
        "location": os.path.dirname(os.path.abspath(__file__)),
        "install": _install_kind(),
        "python": f"{platform.python_version()} ({platform.python_implementation()})",
        "executable": sys.executable,
        "platform": platform.platform(),
        "build": get_build_info(),
    }
    if include_dependencies:
        info["dependencies"] = {
            name: _package_version(name) for name in KEY_DEPENDENCIES
        }
    return info


def _format_build(build: Dict[str, Any]) -> Optional[str]:
    commit = build.get("git_commit_hash")
    if not commit:
        return None
    details = [d for d in (build.get("git_branch"), build.get("git_commit_date")) if d]
    if build.get("git_dirty"):
        details.append("uncommitted changes")
    suffix = f" ({', '.join(details)})" if details else ""
    return f"{commit[:12]}{suffix}"


def format_version_info(info: Optional[Dict[str, Any]] = None) -> str:
    """Render ``get_version_info()`` as an aligned, human-readable block."""
    if info is None:
        info = get_version_info()
    rows = []
    location = info["location"]
    if info.get("install"):
        location += f" ({info['install']})"
    rows.append(("location", location))
    build = _format_build(info["build"])
    if build:
        rows.append(("commit", build))
    rows.append(("python", f"{info['python']} {info['executable']}"))
    rows.append(("platform", info["platform"]))
    deps = info.get("dependencies")
    if deps:
        rows.append(("packages", ", ".join(
            f"{name} {version}" for name, version in deps.items() if version
        )))
    width = max(len(label) for label, _ in rows)
    lines = [f"mdvtools {info['mdvtools']}"]
    lines += [f"  {label.ljust(width)}  {value}" for label, value in rows]
    return "\n".join(lines)


def _command_line(argv: Optional[list] = None) -> str:
    argv = list(sys.argv if argv is None else argv)
    if argv:
        argv[0] = os.path.basename(argv[0])
    return shlex.join(argv)


def format_error_report(exc: BaseException, argv: Optional[list] = None) -> str:
    """Render an exception with the context needed to reproduce it."""
    rule = "=" * 72
    summary = f"{type(exc).__name__}: {exc}"
    try:
        environment = format_version_info(get_version_info(include_dependencies=True))
    except Exception as info_exc:  # the report must not itself fail
        environment = f"(version information unavailable: {info_exc!r})"
    tb = "".join(traceback.format_exception(exc)).rstrip()
    return "\n".join([
        rule,
        "mdvtools error report",
        rule,
        f"command:  {_command_line(argv)}",
        f"cwd:      {os.getcwd()}",
        "",
        tb,
        "",
        environment,
        rule,
        f"Error: {summary}",
        f"If this looks like a bug, please include the report above in an issue: {ISSUE_TRACKER_URL}",
    ])


def report_errors(func: Callable[P, R]) -> Callable[P, R]:
    """Print a formatted error report and exit(1) on unexpected exceptions.

    Click's own control-flow exceptions (usage errors, ``Abort``, ``Exit``)
    and ``SystemExit`` / ``KeyboardInterrupt`` pass through untouched.
    """
    @functools.wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        try:
            return func(*args, **kwargs)
        except (click.ClickException, click.exceptions.Abort, click.exceptions.Exit):
            raise
        except Exception as exc:
            # drop this wrapper's own frame from the reported traceback
            if exc.__traceback__ is not None and exc.__traceback__.tb_next is not None:
                exc.__traceback__ = exc.__traceback__.tb_next
            click.echo(format_error_report(exc), err=True)
            sys.exit(1)

    return wrapper
