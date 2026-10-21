"""Changelog helpers for releasing mdvtools, used by .github/workflows/release*.yml.

    python .github/scripts/release.py current
    python .github/scripts/release.py changelog --previous 1.3.0 [--date 2026-10-21] [--fallback-notes notes.md]
    python .github/scripts/release.py notes 1.4.0

Bump the version with `uv version --bump minor` (from python/) before running `changelog`.
`changelog` turns the CHANGELOG.md `[Unreleased]` section into a dated section for the
version now in python/pyproject.toml. Paths are relative to the repository root, so run it
from there.
"""

import argparse
import datetime
import re
import sys
import tomllib
from pathlib import Path

PYPROJECT = Path("python/pyproject.toml")
CHANGELOG = Path("CHANGELOG.md")
REPO_URL = "https://github.com/Taylor-CCB-Group/MDV"


def current_version() -> str:
    return tomllib.loads(PYPROJECT.read_text())["project"]["version"]


def replace_once(pattern: str, replacement: str, text: str) -> str:
    new_text, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        sys.exit(f"could not find {pattern!r} in {CHANGELOG}")
    return new_text


def unreleased_body(changelog: str) -> str:
    match = re.search(r"^## \[Unreleased\]\n(.*?)(?=^## \[)", changelog, re.MULTILINE | re.DOTALL)
    return match.group(1).strip() if match else ""


def date_changelog(previous: str, date: str, fallback_notes: Path | None) -> str:
    new = current_version()
    if new == previous:
        sys.exit(f"python/pyproject.toml is still at {previous}; run `uv version --bump` first")

    changelog = CHANGELOG.read_text()
    heading = f"## [Unreleased]\n\n## [{new}] - {date}\n"
    if not unreleased_body(changelog):
        if fallback_notes is None:
            sys.exit("CHANGELOG.md has nothing under [Unreleased]; add entries before releasing")
        # the release PR reviewer rewrites these generated PR titles into changelog entries
        heading += f"\n{fallback_notes.read_text().strip()}\n"
    changelog = replace_once(r"^## \[Unreleased\]\n", heading, changelog)
    changelog = replace_once(
        r"^\[Unreleased\]: .*$",
        f"[Unreleased]: {REPO_URL}/compare/v{new}...HEAD\n"
        f"[{new}]: {REPO_URL}/compare/v{previous}...v{new}",
        changelog,
    )
    CHANGELOG.write_text(changelog)
    return new


def notes(version: str) -> str:
    changelog = CHANGELOG.read_text()
    match = re.search(
        rf"^## \[{re.escape(version)}\][^\n]*\n(.*?)(?=^## \[|^\[Unreleased\]: )",
        changelog,
        re.MULTILINE | re.DOTALL,
    )
    if not match:
        sys.exit(f"no CHANGELOG.md section for {version}")
    return match.group(1).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("current")
    changelog_parser = commands.add_parser("changelog")
    changelog_parser.add_argument("--previous", required=True, help="the version before the bump")
    changelog_parser.add_argument("--date", default=datetime.date.today().isoformat())
    changelog_parser.add_argument(
        "--fallback-notes", type=Path, help="notes to use when [Unreleased] is empty"
    )
    notes_parser = commands.add_parser("notes")
    notes_parser.add_argument("version")
    args = parser.parse_args()

    if args.command == "current":
        print(current_version())
    elif args.command == "changelog":
        print(date_changelog(args.previous, args.date, args.fallback_notes))
    else:
        print(notes(args.version))


if __name__ == "__main__":
    main()
