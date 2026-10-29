#!/usr/bin/env python3
"""Add an All Cells view to an MDV project."""

from __future__ import annotations

import argparse
import json
import os
import random
import string
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

VIEW_NAME = "All Cells"
CATEGORICAL_TYPES = {"text", "text16", "multitext"}


def chart_id() -> str:
    return "".join(random.choices(string.ascii_letters + string.digits, k=6))


def find_field(columns: list[dict[str, Any]], display_name: str) -> str:
    matches = [
        str(column["field"])
        for column in columns
        if column.get("name") == display_name and column.get("field")
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one cells column named {display_name!r}; found {len(matches)}."
        )
    return matches[0]


def build_view(datasources: list[dict[str, Any]]) -> dict[str, Any]:
    cells = next((item for item in datasources if item.get("name") == "cells"), None)
    if cells is None:
        raise ValueError("Project has no 'cells' datasource.")
    columns = cells.get("columns", [])
    categorical_fields = [
        str(column["field"])
        for column in columns
        if column.get("datatype") in CATEGORICAL_TYPES and column.get("field")
    ]
    if not categorical_fields:
        raise ValueError("The cells datasource contains no categorical columns.")

    x_field = find_field(columns, "UMAP 2_1")
    y_field = find_field(columns, "UMAP 2_2")
    compartment_field = find_field(columns, "Cell Compartment")
    subgroup_field = find_field(columns, "Cell Subgroup")

    selection_dialog = {
        "type": "selection_dialog",
        "title": "Cell filters",
        "param": categorical_fields,
        "size": [360, 720],
        "position": [0, 0],
        "gsposition": [0, 0],
        "gssize": [3, 8],
        "id": chart_id(),
        "legend": "",
    }
    umap = {
        "type": "wgl_scatter_plot_dev",
        "title": "Cell Compartment UMAP",
        "param": [x_field, y_field],
        "color_by": compartment_field,
        "axis": {
            "x": {"rotate_labels": False, "size": 30, "tickfont": 10},
            "y": {"rotate_labels": False, "size": 45, "tickfont": 10},
        },
        "color_legend": {"display": True, "pos": [12, 12]},
        "field_legend": {"display": True},
        "dimension": "2d",
        "point_shape": "circle",
        "radius": 10,
        "course_radius": 1,
        "opacity": 0.8,
        "on_filter": "hide",
        "zoom_on_filter": False,
        "tooltip": {"show": True, "column": [compartment_field]},
        "size": [720, 720],
        "position": [0, 0],
        "gsposition": [3, 0],
        "gssize": [6, 8],
        "id": chart_id(),
        "legend": "",
    }
    subgroup_chart = {
        "type": "row_chart",
        "title": "Cell Subgroup",
        "param": [subgroup_field],
        "axis": {"x": {"label": "", "size": 25, "textSize": 13, "tickfont": 10}},
        "size": [360, 720],
        "position": [0, 0],
        "gsposition": [9, 0],
        "gssize": [3, 8],
        "id": chart_id(),
        "legend": "",
    }
    return {
        "dataSources": {"cells": {"layout": "gridstack", "panelWidth": 100}},
        "initialCharts": {"cells": [selection_dialog, umap, subgroup_chart]},
    }


def atomic_write_json(path: Path, data: Any) -> None:
    fd, temporary_path = tempfile.mkstemp(prefix=f"{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2)
            handle.write("\n")
        os.replace(temporary_path, path)
    except Exception:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise


def add_all_cells_view(project_path: Path, *, dry_run: bool = False) -> dict[str, Any]:
    project_path = project_path.expanduser().resolve()
    views_path = project_path / "views.json"
    state_path = project_path / "state.json"
    datasources_path = project_path / "datasources.json"
    for path in (views_path, state_path, datasources_path):
        if not path.is_file():
            raise FileNotFoundError(f"Not a complete MDV project; missing {path}")

    datasources = json.loads(datasources_path.read_text(encoding="utf-8"))
    view = build_view(datasources)
    views = json.loads(views_path.read_text(encoding="utf-8"))
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if not isinstance(views, dict):
        raise TypeError("views.json must contain an object.")
    if not isinstance(state, dict):
        raise TypeError("state.json must contain an object.")

    all_views = list(state.get("all_views") or [])
    if VIEW_NAME not in all_views:
        all_views.append(VIEW_NAME)
    views[VIEW_NAME] = view
    state["all_views"] = all_views

    summary = {
        "project": str(project_path),
        "view": VIEW_NAME,
        "categorical_filters": len(view["initialCharts"]["cells"][0]["param"]),
        "charts": [chart["type"] for chart in view["initialCharts"]["cells"]],
        "dry_run": dry_run,
    }
    if dry_run:
        return summary

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    for path in (views_path, state_path):
        backup_path = path.with_name(f"{path.name}.bak.{timestamp}")
        backup_path.write_bytes(path.read_bytes())
    atomic_write_json(views_path, views)
    atomic_write_json(state_path, state)
    summary["backups"] = [
        str(views_path.with_name(f"{views_path.name}.bak.{timestamp}")),
        str(state_path.with_name(f"{state_path.name}.bak.{timestamp}")),
    ]
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True, help="MDV project directory.")
    parser.add_argument("--dry-run", action="store_true", help="Validate and summarize the view without writing files.")
    args = parser.parse_args()
    summary = add_all_cells_view(args.project, dry_run=args.dry_run)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
