#!/usr/bin/env python3
"""Snapshot and restore MDV project JSON files.

Copies ``datasources.json`` and, when present, ``views.json`` and ``state.json``
into ``history/<revision>/`` inside each project. Does not read or write HDF5,
Zarr, or the database, and does not import mdvtools.

A directory containing ``datasources.json`` is a project. The walk does not
descend into a project, so ``history/`` and other nested data are not treated
as further projects. ``snapshot``, ``log``, and ``restore`` accept either one
project directory or a parent directory of projects.

One ``snapshot`` invocation stamps the same ``run`` id on every project. Revision
ids differ per project. Restore a single project by revision id, or every
project from one snapshot by ``--run``.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import secrets
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

SNAPSHOT_FILES = ("datasources.json", "views.json", "state.json")
_REV_ID = re.compile(r"^[0-9]{8}T[0-9]{6}Z-[0-9a-f]{6}$")
_RUN_ID = re.compile(r"^[0-9a-f]{16}$")


class HistoryError(Exception):
    """A history command cannot proceed, and no files should change."""


def discover_projects(roots: list[Path]) -> tuple[list[Path], list[str]]:
    """Return project directories and discovery errors.

    A directory containing ``datasources.json`` is a project. Children of that
    directory are not walked.
    """
    projects: list[Path] = []
    seen: set[Path] = set()
    errors: list[str] = []
    for root in roots:
        if not root.exists():
            errors.append(f"{root}: not found")
            continue
        if not root.is_dir():
            errors.append(f"{root}: not a directory")
            continue
        try:
            for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
                if "datasources.json" not in filenames:
                    continue
                path = Path(dirpath).resolve()
                if path not in seen:
                    seen.add(path)
                    projects.append(path)
                dirnames.clear()
        except OSError as exc:
            errors.append(f"{root}: {exc}")
    projects.sort()
    return projects, errors


def _require_projects(roots: list[Path]) -> list[Path]:
    projects, errors = discover_projects(roots)
    if errors:
        raise HistoryError("\n".join(errors))
    if not projects:
        raise HistoryError("no projects found")
    return projects


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _stamp(when: datetime) -> str:
    return when.strftime("%Y%m%dT%H%M%SZ")


def _created_at(when: datetime) -> str:
    return when.strftime("%Y-%m-%dT%H:%M:%SZ")


def _new_run_id() -> str:
    return secrets.token_hex(8)


def _new_revision_id(when: datetime) -> str:
    return f"{_stamp(when)}-{secrets.token_hex(3)}"


def _history_dir(project: Path) -> Path:
    return project / "history"


def _read_head(project: Path) -> str | None:
    path = _history_dir(project) / "HEAD"
    if not path.is_file():
        return None
    text = path.read_text(encoding="utf-8").strip()
    if not _REV_ID.fullmatch(text):
        return None
    return text


def _replace_bytes(path: Path, data: bytes) -> None:
    """Write ``data`` onto ``path`` via a temp file in the same directory."""
    mode = path.stat().st_mode if path.exists() else 0o664
    directory = path.parent
    directory.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=directory, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(tmp_name, mode)
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def _write_json(path: Path, payload: dict[str, object]) -> None:
    _replace_bytes(path, (json.dumps(payload, indent=2) + "\n").encode("utf-8"))


def _read_manifest(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HistoryError(f"{path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise HistoryError(f"{path}: manifest must be a JSON object")
    return payload


def _manifest_files(manifest: dict[str, object], label: str) -> list[str]:
    files = manifest.get("files")
    if not isinstance(files, list) or not all(isinstance(name, str) for name in files):
        raise HistoryError(f"{label}: manifest files must be a list of names")
    unknown = [name for name in files if name not in SNAPSHOT_FILES]
    if unknown:
        raise HistoryError(f"{label}: manifest lists unknown files: {', '.join(unknown)}")
    ordered = [name for name in SNAPSHOT_FILES if name in files]
    if "datasources.json" not in ordered:
        raise HistoryError(f"{label}: revision is missing datasources.json")
    return ordered


def list_revisions(project: Path) -> list[dict[str, object]]:
    """Return manifests for ``project``, newest first."""
    history = _history_dir(project)
    if not history.is_dir():
        return []
    manifests: list[dict[str, object]] = []
    for child in history.iterdir():
        manifest_path = child / "manifest.json"
        if child.is_dir() and manifest_path.is_file():
            manifests.append(_read_manifest(manifest_path))
    manifests.sort(
        key=lambda item: (str(item.get("created_at", "")), str(item.get("id", ""))),
        reverse=True,
    )
    return manifests


def _revision_dir(project: Path, rev_id: str) -> Path:
    if not _REV_ID.fullmatch(rev_id):
        raise HistoryError(f"{project}: revision {rev_id} not found")
    return _history_dir(project) / rev_id


def snapshot_project(project: Path, *, run: str, message: str, when: datetime) -> str:
    """Copy the live JSON files into a new revision and move HEAD to it."""
    if not _RUN_ID.fullmatch(run):
        raise HistoryError(f"invalid run id: {run}")
    if not (project / "datasources.json").is_file():
        raise HistoryError(f"{project}: not a project")
    history = _history_dir(project)
    history.mkdir(parents=True, exist_ok=True)
    rev_id = _new_revision_id(when)
    dest = history / rev_id
    while dest.exists():
        rev_id = _new_revision_id(when)
        dest = history / rev_id
    dest.mkdir()
    copied: list[str] = []
    for name in SNAPSHOT_FILES:
        src = project / name
        if src.is_file():
            (dest / name).write_bytes(src.read_bytes())
            copied.append(name)
    if "datasources.json" not in copied:
        raise HistoryError(f"{project}: not a project")
    manifest: dict[str, object] = {
        "id": rev_id,
        "parent": _read_head(project),
        "run": run,
        "created_at": _created_at(when),
        "message": message,
        "files": copied,
    }
    _write_json(dest / "manifest.json", manifest)
    _replace_bytes(history / "HEAD", (rev_id + "\n").encode("utf-8"))
    return rev_id


def snapshot_projects(projects: list[Path], message: str) -> tuple[str, list[tuple[Path, str]]]:
    run = _new_run_id()
    when = _utc_now()
    written: list[tuple[Path, str]] = []
    for project in projects:
        written.append((project, snapshot_project(project, run=run, message=message, when=when)))
    return run, written


def _format_revision(manifest: dict[str, object]) -> str:
    parent = manifest.get("parent")
    parent_text = parent if isinstance(parent, str) and parent else "-"
    run = manifest.get("run")
    run_text = run if isinstance(run, str) else "-"
    created = manifest.get("created_at")
    created_text = created if isinstance(created, str) else "-"
    message = manifest.get("message")
    message_text = message if isinstance(message, str) else ""
    rev = manifest.get("id")
    rev_text = rev if isinstance(rev, str) else "-"
    return f"{rev_text}  {created_text}  parent={parent_text}  run={run_text}  {message_text}"


def format_log(projects: list[Path]) -> str:
    blocks: list[str] = []
    show_path = len(projects) != 1
    for project in projects:
        lines = [_format_revision(manifest) for manifest in list_revisions(project)]
        if show_path:
            body = "\n".join(lines)
            blocks.append(f"{project}\n{body}" if body else str(project))
        else:
            blocks.append("\n".join(lines))
    text = "\n\n".join(blocks)
    if text:
        text += "\n"
    return text


def _load_revision(project: Path, rev_id: str) -> tuple[dict[str, object], list[str]]:
    rev_dir = _revision_dir(project, rev_id)
    manifest_path = rev_dir / "manifest.json"
    if not manifest_path.is_file():
        raise HistoryError(f"{project}: revision {rev_id} not found")
    manifest = _read_manifest(manifest_path)
    files = _manifest_files(manifest, f"{project}: revision {rev_id}")
    missing = [name for name in files if not (rev_dir / name).is_file()]
    if missing:
        raise HistoryError(
            f"{project}: revision {rev_id} is missing {', '.join(missing)}"
        )
    return manifest, files


def _find_run(project: Path, run: str) -> str | None:
    matches: list[str] = []
    for manifest in list_revisions(project):
        if manifest.get("run") != run:
            continue
        rev = manifest.get("id")
        if isinstance(rev, str):
            matches.append(rev)
    if not matches:
        return None
    if len(matches) > 1:
        raise HistoryError(f"{project}: run {run} matches more than one revision")
    return matches[0]


def restore_project(project: Path, rev_id: str, *, run: str, message: str) -> str:
    """Copy ``rev_id`` onto the live JSON files, then snapshot that state."""
    _manifest, files = _load_revision(project, rev_id)
    rev_dir = _revision_dir(project, rev_id)
    for name in files:
        _replace_bytes(project / name, (rev_dir / name).read_bytes())
    return snapshot_project(project, run=run, message=message, when=_utc_now())


def _direct_project(root: Path) -> Path | None:
    """Return ``root`` when it is itself a project directory."""
    if not root.exists():
        raise HistoryError(f"{root}: not found")
    if not root.is_dir():
        raise HistoryError(f"{root}: not a directory")
    resolved = root.resolve()
    if (resolved / "datasources.json").is_file():
        return resolved
    return None


def restore_by_revision(root: Path, rev_id: str) -> list[tuple[Path, str]]:
    project = _direct_project(root)
    if project is None:
        raise HistoryError(
            f"{root}: revision ids belong to one project; use --run to restore a directory of projects"
        )
    _load_revision(project, rev_id)
    run = _new_run_id()
    restore_project(project, rev_id, run=run, message=f"restore {rev_id}")
    return [(project, rev_id)]


def restore_by_run(roots: list[Path], run: str) -> tuple[list[tuple[Path, str]], list[Path]]:
    if not _RUN_ID.fullmatch(run):
        raise HistoryError(f"invalid run id: {run}")
    projects = _require_projects(roots)
    planned: list[tuple[Path, str]] = []
    skipped: list[Path] = []
    for project in projects:
        rev_id = _find_run(project, run)
        if rev_id is None:
            skipped.append(project)
            continue
        _load_revision(project, rev_id)
        planned.append((project, rev_id))
    if not planned:
        raise HistoryError(f"no project has run {run}")
    restore_run = _new_run_id()
    for project, rev_id in planned:
        restore_project(project, rev_id, run=restore_run, message=f"restore {rev_id}")
    return planned, skipped


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Snapshot and restore datasources.json, views.json, and state.json "
            "for MDV projects. HDF5 and Zarr files are left unchanged."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    snapshot = subparsers.add_parser("snapshot", help="Copy project JSON files into history/")
    snapshot.add_argument(
        "roots",
        nargs="+",
        help="Project directory, or a parent directory whose children are projects",
    )
    snapshot.add_argument(
        "--message",
        default="snapshot",
        help="Text stored on each new revision (default: snapshot)",
    )

    log = subparsers.add_parser("log", help="List revisions, newest first")
    log.add_argument(
        "roots",
        nargs="+",
        help="Project directory, or a parent directory whose children are projects",
    )

    restore = subparsers.add_parser(
        "restore",
        help="Put JSON files back from a revision id or a shared run id",
    )
    restore.add_argument(
        "root",
        help="One project directory, or a parent directory when using --run",
    )
    restore.add_argument(
        "revision",
        nargs="?",
        help="Revision id to restore. Only valid when root is a single project",
    )
    restore.add_argument(
        "--run",
        help="Restore every project under root that was written by this snapshot run",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "snapshot":
            projects = _require_projects([Path(root) for root in args.roots])
            run, written = snapshot_projects(projects, args.message)
            print(f"run {run}")
            for project, rev_id in written:
                print(f"{project} {rev_id}")
            return 0
        if args.command == "log":
            projects = _require_projects([Path(root) for root in args.roots])
            sys.stdout.write(format_log(projects))
            return 0
        if args.command == "restore":
            if args.revision and args.run:
                parser.error("pass a revision id or --run, not both")
            if not args.revision and not args.run:
                parser.error("pass a revision id or --run")
            if args.run:
                restored, skipped = restore_by_run([Path(args.root)], args.run)
                for project in skipped:
                    print(f"skip {project}: no revision for run {args.run}", file=sys.stderr)
                for project, rev_id in restored:
                    print(f"{project} {rev_id}")
                return 0
            restored = restore_by_revision(Path(args.root), args.revision)
            for project, rev_id in restored:
                print(f"{project} {rev_id}")
            return 0
    except HistoryError as exc:
        print(exc, file=sys.stderr)
        return 1
    parser.error(f"unknown command {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
