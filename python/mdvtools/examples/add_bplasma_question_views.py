#!/usr/bin/env python3
"""Add the two B-cell question views to carth_BPlasma_mu2mdv.

Backs up views.json and state.json, then adds the views that can be built from
the gene list below. datafile.h5 and column metadata are not modified. Display
names from the Data Cleaning sheet are used only to choose Cell Type and Tissue;
they are not written back.

XPB1 is not in this project's gene table. The list uses XBP1.

The app reads the project under ondemand/projects. Write there to see the views
in Select View. Omit --write to print the plan without changing files.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from mdvtools.examples.curate_mdv_project import apply_planned_names, plan_question_views
from mdvtools.examples.rename_softdelete_mdv_project_fields import (
    apply_plan,
    build_plan,
    load_views,
    resolve_datasources,
)
from mdvtools.template_views.create_question_views import atomic_write_json

# XPB1 is absent from the rna table; XBP1 is the gene that is present.
BPLASMA_GENES = [
    "CD19",
    "MS4A1",
    "TCL1A",
    "IGHD",
    "IGHM",
    "CD27",
    "IGHG1",
    "IGHA1",
    "TBX21",
    "FCRL5",
    "MZB1",
    "DERL3",
    "XBP1",
]

DEFAULT_PROJECT = Path("/well/cartography/shared/MDV/ondemand/projects/carth_BPlasma_mu2mdv")
DEFAULT_XLSX = Path("/well/cartography/shared/MDV/MIM_Nov25.xlsx")
DEFAULT_PROJECT_ID = "carth_BPlasma_mu2mdv"


def backup_file(path: Path) -> Path:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    backup = path.with_name(f"{path.name}.bak.{timestamp}")
    backup.write_bytes(path.read_bytes())
    return backup


def planned_datasources(project: Path, excel_path: Path, project_id: str) -> list[dict[str, Any]]:
    source = json.loads((project / "datasources.json").read_text(encoding="utf-8"))
    outcomes: dict[str, list[dict[str, Any]]] = {}
    for datasource in resolve_datasources(project, excel_path, project_id):
        plan, _warnings, _errors, _kept, _ambiguous = build_plan(
            project,
            datasource,
            excel_path,
            project_id,
            source_project_dir=project,
        )
        _summary, results, _actions = apply_plan(project, datasource, plan, True)
        outcomes[datasource] = results
    return apply_planned_names(source, outcomes)


def copy_project(source: Path, output: Path) -> None:
    """Copy the project. datafile.h5 is hard-linked when it is on the same filesystem."""
    if output.exists():
        raise FileExistsError(f"Output project already exists: {output}")
    output.mkdir(parents=True)
    for path in source.iterdir():
        target = output / path.name
        if path.is_dir():
            shutil.copytree(path, target, copy_function=shutil.copy2)
        elif path.name == "datafile.h5" and path.stat().st_dev == output.stat().st_dev:
            os.link(path, target)
        else:
            shutil.copy2(path, target)


def add_views(project: Path, question_plan: dict[str, Any]) -> list[str]:
    views = load_views(project)
    added = list(question_plan["added"])
    for name in added:
        views[name] = question_plan["views"][name]
    atomic_write_json(project / "views.json", views)
    state_path = project / "state.json"
    state = json.loads(state_path.read_text(encoding="utf-8"))
    all_views = list(state.get("all_views") or [])
    for name in added:
        if name not in all_views:
            all_views.append(name)
    state["all_views"] = all_views
    atomic_write_json(state_path, state)
    return added


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT, help="MDV project to read.")
    parser.add_argument("--project-id", default=DEFAULT_PROJECT_ID, help="Excel Project id.")
    parser.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX, help="Workbook with the Data Cleaning sheet.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Write a new project here. Omit to update --project in place.",
    )
    parser.add_argument("--write", action="store_true", help="Back up views and state, then add the views.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project = args.project.expanduser().resolve()
    excel_path = args.xlsx.expanduser().resolve()
    datasources = planned_datasources(project, excel_path, args.project_id)
    question_plan = plan_question_views(
        project,
        datasources,
        load_views(project),
        marker_genes=BPLASMA_GENES,
    )
    print("Genes:", ", ".join(question_plan.get("markers") or []))
    print("Views to add:")
    for name in question_plan["added"]:
        print(f"  add: {name}")
    if not question_plan["added"]:
        print("  (none)")
    print("Skipped:")
    for item in question_plan["skipped"]:
        print(f"  skip: {item['question']}: {item['reason']}")
    if not question_plan["added"]:
        print("No new views to add.")
        return 0
    if not args.write:
        print("Dry run. Pass --write to back up views and state, then add the views.")
        return 0

    target = project
    if args.output is not None:
        target = args.output.expanduser().resolve()
        copy_project(project, target)
    views_backup = backup_file(target / "views.json")
    state_backup = backup_file(target / "state.json")
    added = add_views(target, question_plan)
    print(f"Views backup: {views_backup}")
    print(f"State backup: {state_backup}")
    print(f"Wrote {len(added)} view(s) to {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
