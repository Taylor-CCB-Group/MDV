#!/usr/bin/env python3
"""
Optional pre-merge check: list git-changed files outside the ChatMDV boundary.

See CHATMDV_BOUNDARY.md. Exit code 1 if any disallowed path is changed.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

# Repo root: .../python/mdvtools/llm -> parents[3] = repo root if structure is repo/python/mdvtools/llm
_REPO_ROOT = Path(__file__).resolve().parents[3]

ALLOWED_PREFIXES: tuple[str, ...] = ("python/mdvtools/llm/",)

ALLOWED_FILES: frozenset[str] = frozenset(
    {
        "python/mdvtools/tests/test_chat_first_text_table_policy.py",
        "python/mdvtools/tests/test_datasource_roles.py",
        "python/mdvtools/tests/test_column_field_resolve.py",
        "python/mdvtools/tests/test_code_execution.py",
    }
)


def _normalize(path: str) -> str:
    return path.replace("\\", "/").lstrip("./")


def is_allowed(relative_path: str) -> bool:
    p = _normalize(relative_path)
    if any(p.startswith(prefix) for prefix in ALLOWED_PREFIXES):
        return True
    return p in ALLOWED_FILES


def _git_changed_files(git_range: str | None, staged_only: bool) -> list[str]:
    repo = _REPO_ROOT
    if git_range:
        cmd = ["git", "-C", str(repo), "diff", "--name-only", git_range]
    elif staged_only:
        cmd = ["git", "-C", str(repo), "diff", "--name-only", "--cached"]
    else:
        cmd = ["git", "-C", str(repo), "diff", "--name-only", "HEAD"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if out.returncode != 0:
        print(out.stderr or out.stdout or "git diff failed", file=sys.stderr)
        sys.exit(2)
    lines = [ln.strip() for ln in (out.stdout or "").splitlines() if ln.strip()]
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--git-range",
        default=None,
        metavar="REV",
        help="Compare with git range (e.g. origin/main...HEAD). If set, overrides --staged.",
    )
    parser.add_argument(
        "--staged",
        action="store_true",
        help="Only check staged files (ignored if --git-range is set).",
    )
    args = parser.parse_args()

    try:
        changed = _git_changed_files(args.git_range, args.staged)
    except FileNotFoundError:
        print("git not found; skipping boundary check", file=sys.stderr)
        sys.exit(0)

    if not changed:
        print("ChatMDV boundary: no changed files detected.")
        sys.exit(0)

    bad = [p for p in changed if not is_allowed(p)]
    if not bad:
        print("ChatMDV boundary: all changed files are allowed.")
        sys.exit(0)

    print("ChatMDV boundary: disallowed changed files:", file=sys.stderr)
    for p in bad:
        print(f"  {p}", file=sys.stderr)
    print(
        "\nAllowed: prefixes "
        + ", ".join(repr(p) for p in ALLOWED_PREFIXES)
        + " or "
        + ", ".join(sorted(ALLOWED_FILES)),
        file=sys.stderr,
    )
    sys.exit(1)


if __name__ == "__main__":
    main()
