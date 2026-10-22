#!/usr/bin/env python3
"""Rename or soft-delete MDV project fields using the Data Cleaning sheet."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


EXPECTED_HEADERS = {
    "displayed name": {"displayed name", "displayed_name"},
    "new_display_name": {"new_display_name", "new display name", "new-display-name"},
    "Final Column Name": {"final column name", "final_column_name", "final-column-name"},
    "original name": {"original name", "original_name", "original-name"},
    "Project Id": {"project id", "project_id", "project-id", "project_name", "project name", "mdv_project", "mdv-project"},
    "datasource": {"datasource"},
}
CELL_RENAME_OVERRIDES = {
    "bucket_caf": "Cell Compartment",
    "sub_bucket_caf": "Cell Subgroup",
    "annotation_caf": "Cell State",
}
CELL_RENAME_EXCLUSIONS = {"bucket", "sub_bucket", "annotation"}


def normalize_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"nan", "none", "null"} else text


def normalize_key(value: Any) -> str:
    return normalize_text(value).casefold()


def project_matches(value: Any, project_name: str) -> bool:
    text = normalize_text(value)
    if not text:
        return False
    target = normalize_key(project_name)
    return target in {
        normalize_key(text),
        normalize_key(Path(text).name),
        normalize_key(Path(text).stem),
    }


def duplicated_name_keys(values: list[str]) -> set[str]:
    """Case-folded names that occur more than once. Blanks are ignored."""
    counts = Counter(normalize_key(value) for value in values if normalize_text(value))
    return {key for key, count in counts.items() if count > 1}


def sheet_target_name(
    final_name: str,
    new_name: str,
    displayed_name: str,
    duplicate_finals: set[str],
    duplicate_news: set[str],
) -> tuple[str, str]:
    """Pick the display name for one sheet row.

    A unique Final Column Name wins. A unique new_display_name is used when Final
    is blank or duplicated. When the name that would have been used is duplicated
    within the project and datasource, the row's displayed name is used instead.
    """
    final_is_unique = bool(final_name) and normalize_key(final_name) not in duplicate_finals
    new_is_unique = bool(new_name) and normalize_key(new_name) not in duplicate_news
    if final_is_unique:
        return final_name, "workbook"
    if new_is_unique:
        return new_name, "workbook"
    if displayed_name and (final_name or new_name):
        return displayed_name, "displayed_name"
    if final_name or new_name:
        return "", "displayed_name_blank"
    return "", "blank"


def find_header(df: pd.DataFrame, candidates: set[str]) -> str | None:
    for actual in df.columns:
        if normalize_key(actual) in candidates:
            return str(actual)
    return None


def read_rows_from_excel(
    excel_path: Path, project_name: str | None = None
) -> tuple[str, pd.DataFrame, str, str]:
    workbook = pd.ExcelFile(excel_path)
    sheet_name = next(
        (name for name in workbook.sheet_names if normalize_key(name) == "data cleaning"),
        None,
    )
    if sheet_name is None:
        raise ValueError(f"Worksheet 'Data Cleaning' was not found in {excel_path}.")

    df = pd.read_excel(excel_path, sheet_name=sheet_name, dtype=object)
    project_header = find_header(df, EXPECTED_HEADERS["Project Id"])
    datasource_header = str(df.columns[0]) if len(df.columns) else None
    displayed_header = find_header(df, EXPECTED_HEADERS["displayed name"])
    final_header = find_header(df, EXPECTED_HEADERS["Final Column Name"])
    new_header = find_header(df, EXPECTED_HEADERS["new_display_name"])
    if (
        not all((project_header, datasource_header, displayed_header))
        or normalize_key(datasource_header) != "datasource"
        or not (final_header or new_header)
    ):
        raise ValueError(
            f"Worksheet '{sheet_name}' is missing required columns. Found: {list(df.columns)}"
        )

    if project_name is not None and not any(
        project_matches(value, project_name) for value in df[project_header].tolist()
    ):
        raise ValueError(
            f"Project '{project_name}' was not found in column '{project_header}' "
            f"of worksheet '{sheet_name}' in {excel_path}."
        )
    return sheet_name, df, project_header, datasource_header


def resolve_datasources(project_dir: Path, excel_path: Path, project_name: str) -> list[str]:
    with (project_dir / "datasources.json").open("r", encoding="utf-8") as handle:
        project_datasources = json.load(handle)
    available = {
        normalize_key(item.get("name")): normalize_text(item.get("name"))
        for item in project_datasources
        if normalize_text(item.get("name"))
    }

    _, df, project_header, datasource_header = read_rows_from_excel(excel_path, project_name)
    requested: list[str] = []
    missing: set[str] = set()
    for _, row in df.iterrows():
        if not project_matches(row.get(project_header), project_name):
            continue
        workbook_name = normalize_text(row.get(datasource_header))
        if not workbook_name:
            continue
        datasource = available.get(normalize_key(workbook_name))
        if datasource is None:
            missing.add(workbook_name)
        elif datasource not in requested:
            requested.append(datasource)

    if missing:
        raise ValueError(
            f"Datasource(s) {sorted(missing)} from the workbook are absent from {project_dir}. "
            f"Available datasources: {sorted(available.values())}"
        )
    if not requested:
        raise ValueError(f"No datasource mappings found for project '{project_name}'.")
    return requested


def view_references(views: dict[str, Any], datasource: str, field: str) -> list[tuple[str, str]]:
    """Charts in saved views whose config contains this field id."""

    def references(value: Any) -> bool:
        if isinstance(value, str):
            return value == field
        if isinstance(value, dict):
            return any(references(item) for item in value.values())
        if isinstance(value, (list, tuple)):
            return any(references(item) for item in value)
        return False

    hits: list[tuple[str, str]] = []
    for view_name, view in views.items():
        if not isinstance(view, dict):
            continue
        charts = (view.get("initialCharts") or {}).get(datasource) or []
        for chart in charts:
            if not isinstance(chart, dict) or not references(chart):
                continue
            title = chart.get("title")
            if isinstance(title, list):
                title = " ".join(str(part) for part in title)
            hits.append((view_name, str(title or chart.get("type") or chart.get("id") or "untitled chart")))
    return hits


def load_views(project_dir: Path) -> dict[str, Any]:
    views_path = project_dir / "views.json"
    if not views_path.is_file():
        return {}
    views = json.loads(views_path.read_text(encoding="utf-8"))
    if not isinstance(views, dict):
        raise TypeError(f"{views_path} must be an object")
    return views


def backup_project_jsons(project_dir: Path, backup_dir: Path) -> None:
    backup_dir.mkdir(parents=True, exist_ok=False)
    for path in project_dir.glob("*.json"):
        shutil.copy2(path, backup_dir / f"{path.stem}.orig.json")


def copy_project_for_mutation(
    source_dir: Path, target_dir: Path, backup_dir: Path | None = None
) -> Path:
    if target_dir.exists():
        raise FileExistsError(f"Target directory already exists: {target_dir}")
    shutil.copytree(source_dir, target_dir)
    final_backup_dir = backup_dir or target_dir.parent / (
        f"{target_dir.name}_json_backup_"
        f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
    )
    backup_project_jsons(target_dir, final_backup_dir)
    print(f"Copied source project: {source_dir} -> {target_dir}")
    print(f"Backup created at: {final_backup_dir}")
    return target_dir


def build_field_lookup(datasource_metadata: dict[str, Any]) -> dict[str, str]:
    return {
        normalize_key(column.get("name")): normalize_text(column.get("field"))
        for column in datasource_metadata.get("columns", [])
        if normalize_text(column.get("field")) and normalize_text(column.get("name"))
    }


def build_plan(
    project_dir: Path,
    datasource: str,
    excel_path: Path,
    project_name: str | None = None,
    source_project_dir: Path | None = None,
) -> tuple[list[dict[str, Any]], list[str], list[str], set[str], set[str]]:
    source_project_dir = source_project_dir or project_dir
    with (source_project_dir / "datasources.json").open("r", encoding="utf-8") as handle:
        source_datasources = json.load(handle)
    with (project_dir / "datasources.json").open("r", encoding="utf-8") as handle:
        target_datasources = json.load(handle)
    source_metadata = next((item for item in source_datasources if item.get("name") == datasource), None)
    target_metadata = next((item for item in target_datasources if item.get("name") == datasource), None)
    if source_metadata is None or target_metadata is None:
        raise ValueError(f"Datasource '{datasource}' is missing from the source or target project.")
    sheet_name, df, project_header, datasource_header = read_rows_from_excel(
        excel_path, project_name
    )
    displayed_header = find_header(df, EXPECTED_HEADERS["displayed name"])
    final_header = find_header(df, EXPECTED_HEADERS["Final Column Name"])
    new_header = find_header(df, EXPECTED_HEADERS["new_display_name"])
    original_header = find_header(df, EXPECTED_HEADERS["original name"])
    raw_rows: list[dict[str, Any]] = []
    for index, row in df.iterrows():
        if project_name is not None and not project_matches(row.get(project_header), project_name):
            continue
        if normalize_text(row.get(datasource_header)) != datasource:
            continue
        displayed_name = normalize_text(row.get(displayed_header))
        original_name = normalize_text(row.get(original_header)) if original_header else ""
        if not displayed_name and not original_name:
            continue
        raw_rows.append({
            "row": index + 2,
            "displayed_name": displayed_name,
            "original_name": original_name,
            "final_name": normalize_text(row.get(final_header)) if final_header else "",
            "new_name": normalize_text(row.get(new_header)) if new_header else "",
        })
    duplicate_finals = duplicated_name_keys([item["final_name"] for item in raw_rows])
    duplicate_news = duplicated_name_keys([item["new_name"] for item in raw_rows])
    sheet_rows: list[dict[str, Any]] = []
    for item in raw_rows:
        target_name, rule = sheet_target_name(
            item["final_name"],
            item["new_name"],
            item["displayed_name"],
            duplicate_finals,
            duplicate_news,
        )
        sheet_rows.append({**item, "target_name": target_name, "rule": rule})

    target_columns = {
        normalize_text(column.get("field")): column
        for column in target_metadata.get("columns", [])
        if normalize_text(column.get("field"))
    }
    plan: list[dict[str, Any]] = []
    kept_fields: set[str] = set()
    ambiguous_fields: set[str] = set()
    warnings: list[str] = []
    errors: list[str] = []
    for source_column in source_metadata.get("columns", []):
        display_name = normalize_text(source_column.get("name"))
        field = normalize_text(source_column.get("field"))
        if not display_name or not field or field not in target_columns:
            continue
        current_name = normalize_text(target_columns[field].get("name"))
        # Field id only. rna:bucket is displayed as "bucket" but is a different column, so the sheet still renames it.
        if datasource == "cells" and field in CELL_RENAME_EXCLUSIONS:
            kept_fields.add(field)
            plan.append({
                "row": None,
                "field": field,
                "display_name": display_name,
                "current_name": current_name,
                "target_name": current_name,
                "rule": "exclusion",
                "action": "keep",
                "status": "noop",
            })
            continue
        if datasource == "cells" and (
            display_name in CELL_RENAME_OVERRIDES or field in CELL_RENAME_OVERRIDES
        ):
            target_name = CELL_RENAME_OVERRIDES.get(display_name) or CELL_RENAME_OVERRIDES[field]
            rule = "override"
            matched_rows: list[dict[str, Any]] = []
        else:
            matched_rows = [
                item for item in sheet_rows
                if item["original_name"] == field
                or item["displayed_name"] in {field, display_name}
            ]
            targets = {item["target_name"] for item in matched_rows if item["target_name"]}
            if len(targets) != 1:
                blank_fallback = any(item["rule"] == "displayed_name_blank" for item in matched_rows)
                if len(targets) > 1:
                    warnings.append(
                        f"Ambiguous workbook targets for field {field!r}: {sorted(targets)}"
                    )
                    ambiguous_fields.add(field)
                elif blank_fallback:
                    warnings.append(
                        f"Displayed name is blank for duplicated sheet names on field {field!r}"
                    )
                    ambiguous_fields.add(field)
                continue
            target_name = next(iter(targets))
            matched_rule = next(
                item["rule"] for item in matched_rows if item["target_name"] == target_name
            )
            rule = matched_rule if matched_rule == "displayed_name" else "workbook"
        if not target_name:
            continue
        kept_fields.add(field)
        if current_name == target_name:
            plan.append({
                "row": matched_rows[0]["row"] if matched_rows else None,
                "field": field,
                "display_name": display_name,
                "current_name": current_name,
                "target_name": target_name,
                "rule": rule,
                "action": "keep",
                "status": "noop",
            })
            continue
        plan.append({
            "row": matched_rows[0]["row"] if matched_rows else None,
            "field": field,
            "display_name": display_name,
            "current_name": current_name,
            "target_name": target_name,
            "rule": rule,
            "action": "rename",
        })

    if project_name is not None and not plan and not kept_fields and not sheet_rows:
        raise ValueError(
            f"No applicable rename rows matched project '{project_name}' and datasource '{datasource}' "
            f"in worksheet '{sheet_name}'."
        )
    return plan, warnings, errors, kept_fields, ambiguous_fields


def build_soft_delete_plan(
    project_dir: Path,
    datasource: str,
    kept_fields: set[str],
    ambiguous_fields: set[str],
    views: dict[str, Any],
) -> list[dict[str, Any]]:
    """Columns the sheet does not keep, unless a saved view still uses them."""
    with (project_dir / "datasources.json").open("r", encoding="utf-8") as handle:
        datasources = json.load(handle)
    metadata = next((item for item in datasources if item.get("name") == datasource), None)
    if metadata is None:
        raise ValueError(f"Datasource '{datasource}' not found in {project_dir}")
    plan: list[dict[str, Any]] = []
    for column in metadata.get("columns", []):
        field = normalize_text(column.get("field"))
        display_name = normalize_text(column.get("name"))
        if not field or field in kept_fields or column.get("deleted"):
            continue
        item: dict[str, Any] = {
            "row": None,
            "field": field,
            "display_name": display_name or field,
            "action": "soft_delete",
        }
        if field in ambiguous_fields:
            item["status"] = "blocked"
            item["error"] = "ambiguous workbook targets"
        elif column.get("sgindex") is not None or column.get("sgtype") is not None:
            spatial_key = "sgindex" if column.get("sgindex") is not None else "sgtype"
            item["status"] = "blocked"
            item["error"] = f"spatial ({spatial_key}) column"
        else:
            hits = view_references(views, datasource, field)
            if hits:
                used_by = "; ".join(f"{view_name}: {label}" for view_name, label in hits)
                item["status"] = "blocked"
                item["error"] = f"used by {used_by}"
            else:
                item["status"] = "planned"
        plan.append(item)
    return plan


def apply_plan(
    project_dir: Path, datasource: str, plan: list[dict[str, Any]], dry_run: bool
) -> tuple[Counter, list[dict[str, Any]], list[str]]:
    summary: Counter = Counter()
    results: list[dict[str, Any]] = []
    log_lines: list[str] = []
    metadata_path = project_dir / "datasources.json"
    project_datasources = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata = next((item for item in project_datasources if item.get("name") == datasource), None)
    if metadata is None:
        raise ValueError(f"Datasource '{datasource}' not found in {metadata_path}")
    columns = metadata.get("columns", [])
    columns_by_field = {normalize_text(column.get("field")): column for column in columns}
    plans_by_field = {item["field"]: item for item in plan}

    final_names: dict[str, str] = {}
    for field, column in columns_by_field.items():
        item = plans_by_field.get(field)
        final_names[field] = item["target_name"] if item else normalize_text(column.get("name"))
    duplicate_fields: dict[str, list[str]] = {}
    fields_by_name: dict[str, list[str]] = {}
    for field, name in final_names.items():
        fields_by_name.setdefault(normalize_key(name), []).append(field)
    for fields in fields_by_name.values():
        if len(fields) > 1:
            for field in fields:
                if field in plans_by_field:
                    duplicate_fields[field] = [other for other in fields if other != field]

    writes = False
    for item in plan:
        result = {**item}
        field = item["field"]
        if item.get("action") == "keep":
            result["status"] = "noop"
            summary["keep"] += 1
            log_lines.append(f"keep: field={field} name={item['current_name']} ({item['rule']})")
            results.append(result)
            continue
        column = columns_by_field.get(field)
        if column is None:
            result["status"] = "blocked"
            result["error"] = "field is missing from target datasource"
        elif field in duplicate_fields:
            result["status"] = "blocked"
            result["error"] = f"target name collides with field(s) {duplicate_fields[field]}"
        else:
            result["status"] = "planned" if dry_run else "renamed"
            summary["rename"] += 1
            if not dry_run:
                column["name"] = item["target_name"]
                writes = True
        if result["status"] == "blocked":
            summary["blocked"] += 1
            log_lines.append(
                f"blocked: field={field} row={item['row']} reason={result['error']}"
            )
        else:
            log_lines.append(
                f"{'dry-run rename' if dry_run else 'renamed'}: field={field} "
                f"{item['display_name']} -> {item['target_name']}"
            )
        results.append(result)

    if writes:
        fd, temporary_path = tempfile.mkstemp(
            prefix="datasources.", suffix=".json.tmp", dir=project_dir
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(project_datasources, handle, indent=2)
                handle.write("\n")
            os.replace(temporary_path, metadata_path)
        except Exception:
            if os.path.exists(temporary_path):
                os.unlink(temporary_path)
            raise
    return summary, results, log_lines


def apply_soft_delete_plan(
    project_dir: Path, datasource: str, plan: list[dict[str, Any]], dry_run: bool
) -> tuple[Counter, list[dict[str, Any]], list[str]]:
    summary: Counter = Counter()
    results: list[dict[str, Any]] = []
    log_lines: list[str] = []
    metadata_path = project_dir / "datasources.json"
    project_datasources = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata = next((item for item in project_datasources if item.get("name") == datasource), None)
    if metadata is None:
        raise ValueError(f"Datasource '{datasource}' not found in {metadata_path}")
    columns_by_field = {
        normalize_text(column.get("field")): column for column in metadata.get("columns", [])
    }
    views = load_views(project_dir)
    writes = False
    for item in plan:
        result = {**item}
        field = item["field"]
        column = columns_by_field.get(field)
        if column is None:
            result["status"] = "blocked"
            result["error"] = "field is missing from target datasource"
        elif result.get("status") == "blocked":
            pass
        else:
            hits = view_references(views, datasource, field)
            if hits:
                used_by = "; ".join(f"{view_name}: {label}" for view_name, label in hits)
                result["status"] = "blocked"
                result["error"] = f"used by {used_by}"
            else:
                result["status"] = "planned" if dry_run else "soft_deleted"
                summary["soft_delete"] += 1
                if not dry_run:
                    column["deleted"] = True
                    writes = True
        if result["status"] == "blocked":
            summary["blocked"] += 1
            log_lines.append(f"blocked soft-delete: field={field} reason={result.get('error', '')}")
        else:
            log_lines.append(
                f"{'dry-run soft-delete' if dry_run else 'soft-deleted'}: "
                f"field={field} name={item['display_name']}"
            )
        results.append(result)

    if writes:
        fd, temporary_path = tempfile.mkstemp(
            prefix="datasources.", suffix=".json.tmp", dir=project_dir
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(project_datasources, handle, indent=2)
                handle.write("\n")
            os.replace(temporary_path, metadata_path)
        except Exception:
            if os.path.exists(temporary_path):
                os.unlink(temporary_path)
            raise
    return summary, results, log_lines


def build_report(
    project_dir: Path,
    project_name: str,
    outcomes: dict[str, list[dict[str, Any]]],
    warnings_by_datasource: dict[str, list[str]],
    errors_by_datasource: dict[str, list[str]],
) -> dict[str, Any]:
    report: dict[str, Any] = {
        "project": project_name,
        "target_path": str(project_dir),
        "datasources": {},
    }
    for datasource, results in outcomes.items():
        renamed = []
        kept = []
        soft_deleted = []
        blocked = []
        for item in results:
            entry = {"displayed_name": item["display_name"], "field": item["field"]}
            if item["action"] == "rename":
                entry["target_name"] = item["target_name"]
                entry["current_name"] = item.get("current_name", "")
                if item["status"] in {"renamed", "planned", "noop"}:
                    renamed.append(entry)
            elif item["action"] == "keep":
                entry["name"] = item.get("current_name", item["display_name"])
                entry["rule"] = item.get("rule", "")
                kept.append(entry)
            elif item["status"] in {"soft_deleted", "planned", "noop"}:
                soft_deleted.append(entry)
            if item["status"] == "blocked":
                entry["error"] = item.get("error", "")
                blocked.append(entry)
        report["datasources"][datasource] = {
            "renamed": renamed,
            "kept": kept,
            "soft_deleted": soft_deleted,
            "blocked": blocked,
            "warnings": warnings_by_datasource.get(datasource, []),
            "errors": errors_by_datasource.get(datasource, []),
        }
    return report


def print_report(report: dict[str, Any]) -> None:
    print("\nField curation report")
    print("========================")
    for datasource, payload in report["datasources"].items():
        print(f"\nDatasource: {datasource}")
        print(f"  Renamed/planned: {len(payload['renamed'])}")
        for item in payload["renamed"]:
            print(
                f"    - {item['current_name']!r} -> {item['target_name']!r} (field={item['field']})"
            )
        print(f"  Kept: {len(payload['kept'])}")
        for item in payload["kept"]:
            print(f"    - {item['name']!r} (field={item['field']}, {item['rule']})")
        print(f"  Soft-deleted/planned: {len(payload['soft_deleted'])}")
        for item in payload["soft_deleted"]:
            print(f"    - {item['displayed_name']!r} (field={item['field']})")
        print(f"  Blocked: {len(payload['blocked'])}")
        for item in payload["blocked"]:
            print(f"    - {item['displayed_name']!r} (field={item['field']}): {item['error']}")
        for warning in payload["warnings"]:
            print(f"  Warning: {warning}")
        for error in payload["errors"]:
            print(f"  Error: {error}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path("/well/cartography/shared/MDV/ondemand/projects"))
    parser.add_argument(
        "--source-project",
        type=Path,
        default=None,
        help="Project folder to read. Defaults to SOURCE_DIR/PROJECT_NAME. --project-name remains the Excel Project Id.",
    )
    parser.add_argument("--destination-dir", type=Path, default=Path("/well/cartography/shared/MDV"))
    parser.add_argument("--project-name", default="myeloid_v2_testing", help="Excel Project Id.")
    parser.add_argument("--target-name", default=None)
    parser.add_argument("--xlsx", type=Path, default=Path("/well/cartography/shared/MDV/MIM _ Nov25.xlsx"))
    parser.add_argument("--datasource", default=None, help="Optional single datasource; defaults to all project mappings.")
    parser.add_argument("--backup-dir", type=Path, default=None)
    parser.add_argument("--source-backup-dir", type=Path, default=None)
    parser.add_argument("--dry-run", action="store_true", help="Preview planned renames without changing project metadata.")
    parser.add_argument("--write", action="store_true", help="Apply changes; default is dry-run.")
    parser.add_argument(
        "--copy-project",
        action="store_true",
        help="Copy source project before writing. Ignored during a dry run.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="Where to write the JSON report. A dry run never writes inside the project.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    target_name = args.target_name or f"{args.project_name}_fields_mutated"
    source_project_dir = (args.source_project or (args.source_dir / args.project_name)).resolve()
    project_dir = (args.destination_dir / target_name).resolve()
    excel_path = args.xlsx.resolve()
    dry_run = not args.write
    if not excel_path.is_file():
        raise FileNotFoundError(f"Spreadsheet does not exist: {excel_path}")
    if not source_project_dir.is_dir():
        raise FileNotFoundError(f"Source project does not exist: {source_project_dir}")

    if args.copy_project and not dry_run:
        project_dir = copy_project_for_mutation(source_project_dir, project_dir, args.backup_dir)
    elif dry_run:
        if args.copy_project:
            print("Dry run: not copying the project.")
        project_dir = source_project_dir
    elif not project_dir.is_dir():
        raise FileNotFoundError(f"Project directory does not exist: {project_dir}")

    datasources = [args.datasource] if args.datasource else resolve_datasources(project_dir, excel_path, args.project_name)
    views = load_views(project_dir)
    outcomes: dict[str, list[dict[str, Any]]] = {}
    warnings_by_datasource: dict[str, list[str]] = {}
    errors_by_datasource: dict[str, list[str]] = {}
    plans: dict[str, list[dict[str, Any]]] = {}
    soft_delete_plans: dict[str, list[dict[str, Any]]] = {}
    for datasource in datasources:
        plan, warnings, errors, kept_fields, ambiguous_fields = build_plan(
            project_dir,
            datasource,
            excel_path,
            args.project_name,
            source_project_dir=source_project_dir,
        )
        plans[datasource] = plan
        soft_delete_plans[datasource] = build_soft_delete_plan(
            project_dir, datasource, kept_fields, ambiguous_fields, views
        )
        warnings_by_datasource[datasource] = warnings
        errors_by_datasource[datasource] = errors

    if args.write:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        source_backup = args.source_backup_dir or project_dir.parent / (
            f"{args.project_name}_original_json_backup_{timestamp}"
        )
        backup_project_jsons(source_project_dir, source_backup)
        print(f"Original source JSON backup created at: {source_backup}")
        if not args.copy_project:
            target_backup = args.backup_dir or project_dir.parent / (
                f"{project_dir.name}_prewrite_json_backup_{timestamp}"
            )
            backup_project_jsons(project_dir, target_backup)
            print(f"Target pre-write JSON backup created at: {target_backup}")

    for datasource, plan in plans.items():
        summary, results, actions = apply_plan(project_dir, datasource, plan, dry_run)
        delete_summary, delete_results, delete_actions = apply_soft_delete_plan(
            project_dir, datasource, soft_delete_plans[datasource], dry_run
        )
        summary.update(delete_summary)
        outcomes[datasource] = results + delete_results
        print(
            f"\nDatasource: {datasource}; dry_run={dry_run}; "
            f"renames={len(plan)}; soft_deletes={len(soft_delete_plans[datasource])}; summary={dict(summary)}"
        )
        for action in actions + delete_actions:
            print(f"  - {action}")

    report = build_report(project_dir, args.project_name, outcomes, warnings_by_datasource, errors_by_datasource)
    report["source_path"] = str(source_project_dir)
    report["dry_run"] = dry_run
    print_report(report)
    if args.report is not None:
        report_path = args.report.resolve()
    elif dry_run:
        report_path = source_project_dir.parent / f"{source_project_dir.name}_dry_run_report.json"
    else:
        report_path = project_dir.parent / f"{project_dir.name}_report.json"
    if report_path.parent == project_dir:
        raise ValueError(f"Refusing to write the report inside the project: {report_path}")
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nSaved field report to: {report_path}")
    print("\nDry run complete. No project fields were changed." if dry_run else "\nWrite complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
