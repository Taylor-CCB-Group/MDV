"""Changelog helpers for releasing mdvtools, used by .github/workflows/release*.yml.

    python .github/scripts/release.py current
    python .github/scripts/release.py changelog --previous 1.3.0 [--date 2026-10-21] [--fallback-notes notes.md]
    python .github/scripts/release.py to-rewrite 1.4.0 entries.md
    python .github/scripts/release.py notes 1.4.0

Bump the version with `uv version --bump minor` (from python/) before running `changelog`.
`changelog` turns the CHANGELOG.md `[Unreleased]` section into a dated section for the
version now in python/pyproject.toml. `to-rewrite` adds lines to that section's
`### To rewrite` list, and `notes` refuses to print a section that still has one. Paths are
relative to the repository root, so run it from there.
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
# heading for PR titles a person still has to rewrite into changelog entries
TO_REWRITE = "### To rewrite"


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


def section(changelog: str, version: str) -> re.Match[str]:
    match = re.search(
        rf"^## \[{re.escape(version)}\][^\n]*\n(.*?)(?=^## \[|^\[Unreleased\]: )",
        changelog,
        re.MULTILINE | re.DOTALL,
    )
    if not match:
        sys.exit(f"no CHANGELOG.md section for {version}")
    return match


def add_to_rewrite(version: str, entries: Path) -> None:
    changelog = CHANGELOG.read_text()
    match = section(changelog, version)
    body = match.group(1).rstrip()
    lines = entries.read_text().strip()
    if re.search(rf"^{TO_REWRITE}$", body, re.MULTILINE):
        # add to the end of the existing list, wherever the reviewer has left it
        body = re.sub(
            rf"(^{TO_REWRITE}$.*?)(?=\n\n|\n#|\Z)",
            lambda list_match: f"{list_match.group(1)}\n{lines}",
            body,
            count=1,
            flags=re.MULTILINE | re.DOTALL,
        )
    else:
        body += f"\n\n{TO_REWRITE}\n{lines}"
    CHANGELOG.write_text(f"{changelog[: match.start(1)]}{body}\n\n{changelog[match.end(1) :]}")


def notes(version: str) -> str:
    body = section(CHANGELOG.read_text(), version).group(1)
    if re.search(rf"^{TO_REWRITE}$", body, re.MULTILINE):
        sys.exit(
            f"the [{version}] section of CHANGELOG.md still has a '{TO_REWRITE}' list; "
            "rewrite it into changelog entries and delete the heading"
        )
    return body.strip()


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
    to_rewrite_parser = commands.add_parser("to-rewrite")
    to_rewrite_parser.add_argument("version")
    to_rewrite_parser.add_argument("entries", type=Path, help="lines to add to the list")
    notes_parser = commands.add_parser("notes")
    notes_parser.add_argument("version")
    args = parser.parse_args()

    if args.command == "current":
        print(current_version())
    elif args.command == "changelog":
        print(date_changelog(args.previous, args.date, args.fallback_notes))
    elif args.command == "to-rewrite":
        add_to_rewrite(args.version, args.entries)
    else:
        print(notes(args.version))


if __name__ == "__main__":
    main()
