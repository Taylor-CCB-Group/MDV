#!/usr/bin/env python3
"""Rename or soft-delete MDV project fields using the Data Cleaning sheet."""

from __future__ import annotations

import argparse
import json
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from mdvtools.mdvproject import MDVProject


EXPECTED_HEADERS = {
    "displayed name": {"displayed name", "displayed_name"},
    "new_display_name": {"new_display_name", "new display name", "new-display-name"},
    "Final Column Name": {"final column name", "final_column_name", "final-column-name"},
    "Project Id": {"project id", "project_id", "project-id", "project_name", "project name", "mdv_project", "mdv-project"},
    "datasource": {"datasource"},
}


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
    datasource_header = find_header(df, EXPECTED_HEADERS["datasource"])
    displayed_header = find_header(df, EXPECTED_HEADERS["displayed name"])
    final_header = find_header(df, EXPECTED_HEADERS["Final Column Name"])
    new_header = find_header(df, EXPECTED_HEADERS["new_display_name"])
    if not all((project_header, datasource_header, displayed_header)) or not (final_header or new_header):
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


def backup_project_jsons(project_dir: Path, backup_dir: Path) -> None:
    backup_dir.mkdir(parents=True, exist_ok=True)
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
) -> tuple[list[dict[str, Any]], list[str], list[str]]:
    with (project_dir / "datasources.json").open("r", encoding="utf-8") as handle:
        project_datasources = json.load(handle)
    metadata = next(
        (item for item in project_datasources if normalize_key(item.get("name")) == normalize_key(datasource)),
        None,
    )
    if metadata is None:
        raise ValueError(f"Datasource '{datasource}' not found in {project_dir / 'datasources.json'}")

    sheet_name, df, project_header, datasource_header = read_rows_from_excel(
        excel_path, project_name
    )
    displayed_header = find_header(df, EXPECTED_HEADERS["displayed name"])
    final_header = find_header(df, EXPECTED_HEADERS["Final Column Name"])
    new_header = find_header(df, EXPECTED_HEADERS["new_display_name"])
    field_lookup = build_field_lookup(metadata)
    plan: list[dict[str, Any]] = []
    warnings: list[str] = []
    errors: list[str] = []

    for index, row in df.iterrows():
        if project_name is not None and not project_matches(row.get(project_header), project_name):
            continue
        if normalize_key(row.get(datasource_header)) != normalize_key(datasource):
            continue
        display_name = normalize_text(row.get(displayed_header))
        if not display_name:
            continue
        field = field_lookup.get(normalize_key(display_name))
        if field is None:
            warnings.append(
                f"Row {index + 2}: display name '{display_name}' not found in datasource '{datasource}'"
            )
            continue

        final_name = normalize_text(row.get(final_header)) if final_header else ""
        new_name = normalize_text(row.get(new_header)) if new_header else ""
        target_name = final_name or new_name
        item = {
            "row": index + 2,
            "field": field,
            "display_name": display_name,
            "action": "rename" if target_name else "soft_delete",
        }
        if target_name:
            item["target_name"] = target_name
        plan.append(item)

    if project_name is not None and not plan:
        raise ValueError(
            f"No rows matched project '{project_name}' and datasource '{datasource}' "
            f"in worksheet '{sheet_name}'."
        )
    return plan, warnings, errors


def apply_plan(
    project_dir: Path, datasource: str, plan: list[dict[str, Any]], dry_run: bool
) -> tuple[Counter, list[dict[str, Any]], list[str]]:
    summary: Counter = Counter()
    results: list[dict[str, Any]] = []
    log_lines: list[str] = []
    project = None if dry_run else MDVProject(str(project_dir))

    for item in plan:
        result = {**item, "status": "planned" if dry_run else "blocked"}
        try:
            if item["action"] == "rename":
                if dry_run:
                    changed = True
                else:
                    changed = project.rename_column(datasource, item["field"], item["target_name"])
                if changed:
                    summary["rename"] += 1
                    result["status"] = "planned" if dry_run else "renamed"
                    log_lines.append(
                        f"{'dry-run rename' if dry_run else 'renamed'}: field={item['field']} "
                        f"display_name={item['display_name']} -> {item['target_name']}"
                    )
                else:
                    summary["rename_noop"] += 1
                    result["status"] = "noop"
                    log_lines.append(f"rename_noop: field={item['field']} already has the target name")
            else:
                if dry_run:
                    changed = True
                else:
                    changed = project.soft_delete_column(datasource, item["field"])
                if changed:
                    summary["soft_delete"] += 1
                    result["status"] = "planned" if dry_run else "soft_deleted"
                    log_lines.append(
                        f"{'dry-run soft-delete' if dry_run else 'soft_deleted'}: field={item['field']}"
                    )
                else:
                    summary["soft_delete_noop"] += 1
                    result["status"] = "noop"
                    log_lines.append(f"soft_delete_noop: field={item['field']} already hidden")
        except (ValueError, AttributeError) as exc:
            summary["blocked"] += 1
            result["status"] = "blocked"
            result["error"] = str(exc)
            log_lines.append(f"blocked: field={item['field']} row={item['row']} reason={exc}")
        results.append(result)
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
        soft_deleted = []
        blocked = []
        for item in results:
            entry = {"displayed_name": item["display_name"], "field": item["field"]}
            if item["action"] == "rename":
                entry["target_name"] = item["target_name"]
                if item["status"] in {"renamed", "planned", "noop"}:
                    renamed.append(entry)
            elif item["status"] in {"soft_deleted", "planned", "noop"}:
                soft_deleted.append(entry)
            if item["status"] == "blocked":
                entry["error"] = item.get("error", "")
                blocked.append(entry)
        report["datasources"][datasource] = {
            "renamed": renamed,
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
            print(f"    - {item['displayed_name']!r} -> {item['target_name']!r} (field={item['field']})")
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
    parser.add_argument("--destination-dir", type=Path, default=Path("/well/cartography/shared/MDV"))
    parser.add_argument("--project-name", default="myeloid_v2_testing")
    parser.add_argument("--target-name", default=None)
    parser.add_argument("--xlsx", type=Path, default=Path("/well/cartography/shared/MDV/MIM _ Nov25.xlsx"))
    parser.add_argument("--datasource", default=None, help="Optional single datasource; defaults to all project mappings.")
    parser.add_argument("--backup-dir", type=Path, default=None)
    parser.add_argument("--write", action="store_true", help="Apply changes; default is dry-run.")
    parser.add_argument("--copy-project", action="store_true", help="Copy source project before planning or writing.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    target_name = args.target_name or f"{args.project_name}_fields_mutated"
    source_project_dir = (args.source_dir / args.project_name).resolve()
    project_dir = (args.destination_dir / target_name).resolve()
    excel_path = args.xlsx.resolve()
    if not excel_path.is_file():
        raise FileNotFoundError(f"Spreadsheet does not exist: {excel_path}")

    if args.copy_project:
        if not source_project_dir.is_dir():
            raise FileNotFoundError(f"Source project does not exist: {source_project_dir}")
        project_dir = copy_project_for_mutation(source_project_dir, project_dir, args.backup_dir)
    elif not project_dir.is_dir():
        raise FileNotFoundError(f"Project directory does not exist: {project_dir}")

    dry_run = not args.write
    if args.write and not args.copy_project:
        backup_dir = args.backup_dir or project_dir.parent / (
            f"{project_dir.name}_json_backup_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
        )
        backup_project_jsons(project_dir, backup_dir)
        print(f"Backup created at: {backup_dir}")

    datasources = [args.datasource] if args.datasource else resolve_datasources(project_dir, excel_path, args.project_name)
    outcomes: dict[str, list[dict[str, Any]]] = {}
    warnings_by_datasource: dict[str, list[str]] = {}
    errors_by_datasource: dict[str, list[str]] = {}
    plans: dict[str, list[dict[str, Any]]] = {}
    for datasource in datasources:
        plan, warnings, errors = build_plan(project_dir, datasource, excel_path, args.project_name)
        plans[datasource] = plan
        warnings_by_datasource[datasource] = warnings
        errors_by_datasource[datasource] = errors

    for datasource, plan in plans.items():
        summary, results, actions = apply_plan(project_dir, datasource, plan, dry_run)
        outcomes[datasource] = results
        print(f"\nDatasource: {datasource}; dry_run={dry_run}; rows={len(plan)}; summary={dict(summary)}")
        for action in actions:
            print(f"  - {action}")

    report = build_report(project_dir, args.project_name, outcomes, warnings_by_datasource, errors_by_datasource)
    print_report(report)
    report_path = project_dir.parent / f"{project_dir.name}_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nSaved field report to: {report_path}")
    print("\nDry run complete. No project fields were changed." if dry_run else "\nWrite complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
