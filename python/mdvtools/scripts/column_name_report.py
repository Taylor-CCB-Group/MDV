#!/usr/bin/env python3
"""Report column-name inconsistencies across MDV projects on disk.

Reads ``datasources.json`` only (no HDF5, no database, no mdvtools import).
Projects are directories that contain ``datasources.json``. The walk does not
descend into a project, so nested data files are not treated as further projects.

Columns are compared within a datasource name. ``field`` values that match
after trimming, lowercasing, and dropping every character outside ``[a-z0-9]``
are one spelling cluster. ``cell_type``, ``Cell Type``, and ``cell-type`` all
join ``celltype``. Every group is reported, including a column that uses one
spelling everywhere. A group is inconsistent when it contains more than one
raw ``field`` string.

``--fuzzy RATIO`` additionally merges clusters in the same datasource whose
normalized keys have a ``difflib`` ratio at or above ``RATIO``. Keys shorter
than ``--min-fuzzy-length`` (default 6) are left out of fuzzy merges, so short
tokens such as ``pc1`` / ``pc2`` stay separate. Fuzzy matches are suggestions
and can still be wrong.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from pathlib import Path
from typing import TextIO


def normalize_field(value: str) -> str:
    """Collapse case and separators so near-duplicate field ids compare equal."""
    lowered = value.strip().lower()
    return "".join(ch for ch in lowered if ch.isascii() and ch.isalnum())


@dataclass(frozen=True)
class ColumnHit:
    project: str
    datasource: str
    field_name: str
    display_name: str


@dataclass(frozen=True)
class Occurrence:
    project: str
    display_name: str


@dataclass(frozen=True)
class Spelling:
    field_name: str
    occurrences: tuple[Occurrence, ...]

    @property
    def projects(self) -> tuple[str, ...]:
        return tuple(sorted({item.project for item in self.occurrences}))

    @property
    def display_names(self) -> tuple[str, ...]:
        return tuple(sorted({item.display_name for item in self.occurrences}))


@dataclass(frozen=True)
class Cluster:
    datasource: str
    normalized: tuple[str, ...]
    spellings: tuple[Spelling, ...]

    @property
    def inconsistent(self) -> bool:
        return len(self.spellings) > 1


@dataclass
class Scan:
    roots: list[str]
    projects: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    clusters: list[Cluster] = field(default_factory=list)


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


def _string_field(column: dict[str, object], key: str) -> str | None:
    value = column.get(key)
    if isinstance(value, str) and value != "":
        return value
    return None


def load_project_columns(
    project: Path, *, include_internal: bool
) -> tuple[list[ColumnHit], str | None]:
    """Read column hits from one project's ``datasources.json``."""
    path = project / "datasources.json"
    try:
        text = path.read_text(encoding="utf-8")
        parsed: object = json.loads(text)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return [], f"{path}: {exc}"
    if not isinstance(parsed, list):
        return [], f"{path}: expected a list of datasources"
    hits: list[ColumnHit] = []
    project_name = str(project)
    for item in parsed:
        if not isinstance(item, dict):
            continue
        raw_name = item.get("name")
        datasource = raw_name if isinstance(raw_name, str) else ""
        columns = item.get("columns")
        if not isinstance(columns, list):
            continue
        for column in columns:
            if not isinstance(column, dict):
                continue
            field_name = _string_field(column, "field") or _string_field(column, "name")
            if field_name is None:
                continue
            if not include_internal and field_name.startswith("__"):
                continue
            display_name = _string_field(column, "name") or field_name
            hits.append(
                ColumnHit(
                    project=project_name,
                    datasource=datasource,
                    field_name=field_name,
                    display_name=display_name,
                )
            )
    return hits, None


def _union(parent: list[int], left: int, right: int) -> None:
    left_root = _find(parent, left)
    right_root = _find(parent, right)
    if left_root != right_root:
        parent[right_root] = left_root


def _find(parent: list[int], index: int) -> int:
    while parent[index] != index:
        parent[index] = parent[parent[index]]
        index = parent[index]
    return index


def find_inconsistencies(
    hits: list[ColumnHit],
    *,
    fuzzy: float | None = None,
    min_fuzzy_length: int = 6,
) -> list[Cluster]:
    """Group every column by datasource and normalized field spelling."""
    by_datasource: dict[str, list[ColumnHit]] = defaultdict(list)
    for hit in hits:
        by_datasource[hit.datasource].append(hit)

    clusters: list[Cluster] = []
    for datasource in sorted(by_datasource):
        by_norm: dict[str, list[ColumnHit]] = defaultdict(list)
        for hit in by_datasource[datasource]:
            by_norm[normalize_field(hit.field_name)].append(hit)
        norms = sorted(by_norm)
        parent = list(range(len(norms)))
        if fuzzy is not None:
            for i, left in enumerate(norms):
                if len(left) < min_fuzzy_length:
                    continue
                for j in range(i + 1, len(norms)):
                    right = norms[j]
                    if len(right) < min_fuzzy_length:
                        continue
                    if SequenceMatcher(None, left, right).ratio() >= fuzzy:
                        _union(parent, i, j)
        groups: dict[int, list[int]] = defaultdict(list)
        for index in range(len(norms)):
            groups[_find(parent, index)].append(index)
        for indexes in groups.values():
            grouped_hits: list[ColumnHit] = []
            grouped_norms: list[str] = []
            for index in indexes:
                grouped_norms.append(norms[index])
                grouped_hits.extend(by_norm[norms[index]])
            by_field: dict[str, list[ColumnHit]] = defaultdict(list)
            for hit in grouped_hits:
                by_field[hit.field_name].append(hit)
            spellings: list[Spelling] = []
            for field_name in sorted(by_field):
                field_hits = by_field[field_name]
                seen: set[tuple[str, str]] = set()
                occurrences: list[Occurrence] = []
                for hit in sorted(field_hits, key=lambda item: (item.project, item.display_name)):
                    key = (hit.project, hit.display_name)
                    if key in seen:
                        continue
                    seen.add(key)
                    occurrences.append(Occurrence(project=hit.project, display_name=hit.display_name))
                spellings.append(Spelling(field_name=field_name, occurrences=tuple(occurrences)))
            clusters.append(
                Cluster(
                    datasource=datasource,
                    normalized=tuple(sorted(grouped_norms)),
                    spellings=tuple(spellings),
                )
            )
    clusters.sort(key=lambda cluster: (cluster.datasource, cluster.normalized))
    return clusters


def scan_roots(
    roots: list[Path],
    *,
    include_internal: bool = False,
    fuzzy: float | None = None,
    min_fuzzy_length: int = 6,
) -> Scan:
    projects, errors = discover_projects(roots)
    scan = Scan(roots=[str(root) for root in roots], projects=[str(p) for p in projects], errors=errors)
    hits: list[ColumnHit] = []
    for project in projects:
        project_hits, error = load_project_columns(project, include_internal=include_internal)
        if error is not None:
            scan.errors.append(error)
            continue
        hits.extend(project_hits)
    scan.clusters = find_inconsistencies(
        hits, fuzzy=fuzzy, min_fuzzy_length=min_fuzzy_length
    )
    return scan


def _normalized_label(cluster: Cluster) -> str:
    return " | ".join(cluster.normalized)


def _spelling_count_label(count: int) -> str:
    if count == 1:
        return "1 spelling"
    return f"{count} spellings"


def format_markdown(
    scan: Scan,
    *,
    fuzzy: float | None = None,
    min_fuzzy_length: int = 6,
) -> str:
    inconsistent_count = sum(1 for cluster in scan.clusters if cluster.inconsistent)
    lines: list[str] = ["# Column naming report", ""]
    lines.append(f"- Roots: {', '.join(scan.roots) if scan.roots else '(none)'}")
    lines.append(f"- Projects scanned: {len(scan.projects)}")
    lines.append(f"- Column groups: {len(scan.clusters)}")
    lines.append(f"- Inconsistent groups: {inconsistent_count}")
    lines.append(f"- Unreadable files: {len(scan.errors)}")
    for error in scan.errors:
        lines.append(f"  - {error}")
    lines.append("")
    if fuzzy is None:
        lines.append("Fuzzy matching is off. Only case and separator differences are grouped.")
    else:
        lines.append(
            f"Fuzzy ratio: {fuzzy} (suggestions only; matches can be wrong). "
            f"Minimum key length: {min_fuzzy_length}."
        )
    lines.append("")
    if not scan.clusters:
        lines.append("No columns found.")
        lines.append("")
        return "\n".join(lines)

    current = None
    for cluster in scan.clusters:
        if cluster.datasource != current:
            current = cluster.datasource
            heading = current if current else "(unnamed datasource)"
            lines.append(f"## {heading}")
            lines.append("")
        spelling_label = _spelling_count_label(len(cluster.spellings))
        lines.append(f"### {_normalized_label(cluster)} ({spelling_label})")
        lines.append("")
        lines.append("| field | projects | display names |")
        lines.append("| --- | ---: | --- |")
        for spelling in cluster.spellings:
            differing = [
                name for name in spelling.display_names if name != spelling.field_name
            ]
            display = ", ".join(differing) if differing else "(same as field)"
            lines.append(
                f"| `{spelling.field_name}` | {len(spelling.projects)} | {display} |"
            )
        lines.append("")
        for spelling in cluster.spellings:
            lines.append(f"- `{spelling.field_name}` ({len(spelling.projects)})")
            for project in spelling.projects:
                lines.append(f"  - {project}")
        lines.append("")
    return "\n".join(lines)


def iter_csv_rows(scan: Scan) -> Iterator[dict[str, str]]:
    for cluster in scan.clusters:
        label = _normalized_label(cluster)
        for spelling in cluster.spellings:
            for occurrence in spelling.occurrences:
                yield {
                    "datasource": cluster.datasource,
                    "normalized": label,
                    "field": spelling.field_name,
                    "display_name": occurrence.display_name,
                    "project": occurrence.project,
                    "inconsistent": "true" if cluster.inconsistent else "false",
                }


def write_csv(scan: Scan, destination: TextIO) -> None:
    writer = csv.DictWriter(
        destination,
        fieldnames=[
            "datasource",
            "normalized",
            "field",
            "display_name",
            "project",
            "inconsistent",
        ],
    )
    writer.writeheader()
    for row in iter_csv_rows(scan):
        writer.writerow(row)


def cluster_to_json(cluster: Cluster) -> dict[str, object]:
    return {
        "datasource": cluster.datasource,
        "normalized": list(cluster.normalized),
        "inconsistent": cluster.inconsistent,
        "spellings": [
            {
                "field": spelling.field_name,
                "project_count": len(spelling.projects),
                "projects": list(spelling.projects),
                "display_names": list(spelling.display_names),
                "occurrences": [
                    {"project": item.project, "display_name": item.display_name}
                    for item in spelling.occurrences
                ],
            }
            for spelling in cluster.spellings
        ],
    }


def scan_to_json(scan: Scan, *, fuzzy: float | None) -> dict[str, object]:
    return {
        "roots": scan.roots,
        "project_count": len(scan.projects),
        "projects": scan.projects,
        "unreadable": scan.errors,
        "fuzzy": fuzzy,
        "clusters": [cluster_to_json(cluster) for cluster in scan.clusters],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Report every column field across MDV projects, grouped by datasource "
            "and normalized spelling. Groups with more than one raw field are "
            "marked inconsistent. Fuzzy matches are optional suggestions and "
            "can be wrong."
        )
    )
    parser.add_argument(
        "roots",
        nargs="+",
        help="Project directory, or a parent directory whose children are projects",
    )
    parser.add_argument("--csv", metavar="PATH", help="Also write a CSV report to PATH")
    parser.add_argument("--json", metavar="PATH", help="Also write a JSON report to PATH")
    parser.add_argument(
        "--fuzzy",
        type=float,
        default=None,
        metavar="RATIO",
        help=(
            "Also merge normalized keys in the same datasource with difflib "
            "ratio >= RATIO (0-1). Off by default. Suggestions can be wrong."
        ),
    )
    parser.add_argument(
        "--min-fuzzy-length",
        type=int,
        default=6,
        help="Skip fuzzy merges for normalized keys shorter than this (default: 6)",
    )
    parser.add_argument(
        "--include-internal",
        action="store_true",
        help="Include field names that start with __ (omitted by default)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.fuzzy is not None and not 0.0 <= args.fuzzy <= 1.0:
        parser.error("--fuzzy must be between 0 and 1")
    if args.min_fuzzy_length < 1:
        parser.error("--min-fuzzy-length must be >= 1")
    scan = scan_roots(
        [Path(root) for root in args.roots],
        include_internal=args.include_internal,
        fuzzy=args.fuzzy,
        min_fuzzy_length=args.min_fuzzy_length,
    )
    sys.stdout.write(
        format_markdown(scan, fuzzy=args.fuzzy, min_fuzzy_length=args.min_fuzzy_length)
    )
    if args.csv:
        with open(args.csv, "w", encoding="utf-8", newline="") as handle:
            write_csv(scan, handle)
    if args.json:
        payload = scan_to_json(scan, fuzzy=args.fuzzy)
        Path(args.json).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
