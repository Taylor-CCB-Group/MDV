#!/usr/bin/env python3
"""Add a project overview view and put the project's views in a fixed order.

The overview copies the layout of the "Cell atlas" view in myeloid_v2_testing_latest:
the project title and abstract across the top, UMAPs coloured by cell type, cell
group and tissue, a UMAP coloured by disease, and composition charts below. Charts
use columns that already exist in the cells datasource. datafile.h5 is not opened.

The views are then ordered: the overview first, the question views in question
order, then every other listed view in its current order. ``initial_view`` is set to
the overview. Views that are not in ``all_views`` stay out of it.

The default is a dry run that prints the plan. ``--write`` backs up views.json and
state.json beside them, then writes both.
"""

from __future__ import annotations

import argparse
import json
import random
import string
import sys
from dataclasses import dataclass, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from mdvtools.examples import curate_mdv_project as curate
from mdvtools.template_views import create_question_views as questions
from mdvtools.template_views import create_synovial_question_views as synovial
from mdvtools.template_views.create_question_views import atomic_write_json

VIEW_NAME = "Cell atlas"
CELLS = "cells"
TEXT_TYPES = {"text", "text16", "multitext"}
PLACEHOLDER_TITLE = "Project title"
PLACEHOLDER_ABSTRACT = "Project abstract."

# Question views in question order. Q8 and Q9 have different names per project.
QUESTION_VIEW_ORDER: tuple[str, ...] = tuple(
    dict.fromkeys(
        [
            questions.Q1_VIEW_NAME,
            questions.Q2_VIEW_NAME,
            questions.Q3_VIEW_NAME,
            questions.Q6_VIEW_NAME,
            questions.Q7_VIEW_NAME,
            questions.Q8_VIEW_NAME,
            synovial.Q8_VIEW_NAME,
            curate.Q8_VIEW_NAME,
            questions.Q9_VIEW_NAME,
            synovial.Q9_VIEW_NAME,
            curate.Q9_VIEW_NAME,
        ]
    )
)

# Display names tried in order for each role. The first non-deleted text column
# whose display name matches (ignoring case) is used. --<role> overrides this.
ROLE_NAMES: dict[str, tuple[str, ...]] = {
    "cell_type": ("Cell State", "Cell Type"),
    "cell_group": ("Cell Subgroup", "Cell Group"),
    "disease": ("Diagnosis", "Disease"),
    "tissue": ("Tissue Simple", "Tissue", "Biopsy Site"),
}
UMAP_FIELD_SUFFIXES = ("X_umap_mindist_0.25_1", "X_umap_mindist_0.25_2")
UMAP_DISPLAY_NAMES = ("UMAP 1", "UMAP 2")

# Grid placement and pixel size copied from the myeloid Cell atlas view (12 columns).
WIDE = {"gssize": [12, 2], "size": [1839, 274]}
THIRD = {"gssize": [4, 5], "size": [609, 694]}
HALF = {"gssize": [6, 4], "size": [917, 554]}


@dataclass(frozen=True)
class Column:
    field: str
    name: str


@dataclass(frozen=True)
class Roles:
    umap_x: Column
    umap_y: Column
    cell_type: Column
    cell_group: Column
    disease: Column
    tissue: Column


def chart_id(used: set[str]) -> str:
    while True:
        value = "".join(random.choices(string.ascii_letters + string.digits, k=6))
        if value not in used:
            used.add(value)
            return value


def live_columns(datasource: dict[str, Any]) -> list[Column]:
    columns: list[Column] = []
    for column in datasource.get("columns", []):
        field = column.get("field")
        if not isinstance(field, str) or not field or column.get("deleted"):
            continue
        name = column.get("name")
        columns.append(Column(field=field, name=name if isinstance(name, str) and name else field))
    return columns


def text_columns(datasource: dict[str, Any]) -> list[Column]:
    fields = {
        column.get("field")
        for column in datasource.get("columns", [])
        if column.get("datatype") in TEXT_TYPES
    }
    return [column for column in live_columns(datasource) if column.field in fields]


def find_column(columns: list[Column], wanted: str, role: str) -> Column:
    """Match ``wanted`` against field ids first, then display names (ignoring case)."""
    for column in columns:
        if column.field == wanted:
            return column
    matches = [column for column in columns if column.name.casefold() == wanted.casefold()]
    if len(matches) == 1:
        return matches[0]
    if matches:
        fields = ", ".join(column.field for column in matches)
        raise ValueError(f"{role}: {wanted!r} matches more than one column ({fields}); pass a field id")
    raise ValueError(f"{role}: no visible column has the field id or display name {wanted!r}")


def choose_role(columns: list[Column], role: str, override: str | None) -> Column:
    if override:
        return find_column(columns, override, role)
    for candidate in ROLE_NAMES[role]:
        matches = [column for column in columns if column.name.casefold() == candidate.casefold()]
        if len(matches) == 1:
            return matches[0]
        if matches:
            fields = ", ".join(column.field for column in matches)
            raise ValueError(f"{role}: several columns are named {candidate!r} ({fields}); pass --{role.replace('_', '-')}")
    tried = ", ".join(repr(name) for name in ROLE_NAMES[role])
    raise ValueError(
        f"{role}: no visible text column is named {tried}; pass --{role.replace('_', '-')} with a field id or display name"
    )


def choose_umap(columns: list[Column], overrides: tuple[str | None, str | None]) -> tuple[Column, Column]:
    if overrides[0] or overrides[1]:
        if not (overrides[0] and overrides[1]):
            raise ValueError("umap: pass both --umap-x and --umap-y")
        return find_column(columns, overrides[0], "umap-x"), find_column(columns, overrides[1], "umap-y")
    picked: list[Column] = []
    for suffix, display in zip(UMAP_FIELD_SUFFIXES, UMAP_DISPLAY_NAMES):
        matches = [column for column in columns if column.field.endswith(suffix)]
        if not matches:
            matches = [column for column in columns if column.name.casefold() == display.casefold()]
        if len(matches) != 1:
            found = ", ".join(column.field for column in matches) or "none"
            raise ValueError(
                f"umap: expected one visible column ending in {suffix!r} or named {display!r}; "
                f"found {found}. Pass --umap-x and --umap-y"
            )
        picked.append(matches[0])
    return picked[0], picked[1]


def choose_roles(datasource: dict[str, Any], overrides: dict[str, str | None]) -> Roles:
    texts = text_columns(datasource)
    umap_x, umap_y = choose_umap(live_columns(datasource), (overrides.get("umap_x"), overrides.get("umap_y")))
    return Roles(
        umap_x=umap_x,
        umap_y=umap_y,
        cell_type=choose_role(texts, "cell_type", overrides.get("cell_type")),
        cell_group=choose_role(texts, "cell_group", overrides.get("cell_group")),
        disease=choose_role(texts, "disease", overrides.get("disease")),
        tissue=choose_role(texts, "tissue", overrides.get("tissue")),
    )


def overview_text(title: str, abstract: str) -> str:
    return f"# {title.strip()}\n{abstract.strip()}"


def bar_axes(stacked: bool) -> dict[str, Any]:
    axis = {"label": "", "size": 25, "textSize": 13, "tickfont": 10}
    if not stacked:
        return {"x": dict(axis)}
    return {"ry": dict(axis), "x": dict(axis), "y": {**axis, "size": 40}}


def umap_chart(roles: Roles, color_by: Column, grid: dict[str, Any], position: list[int]) -> dict[str, Any]:
    chart = questions.scatter(
        color_by.name, [roles.umap_x.field, roles.umap_y.field], color_by.field, color_by.field
    )
    chart.update({
        "on_filter": "hide",
        "opacity": 0.8,
        "position": [0, 0],
        "gsposition": position,
        **grid,
    })
    return chart


def stacked_chart(rows: Column, stacks: Column, position: list[int]) -> dict[str, Any]:
    return {
        "type": "stacked_row_chart",
        "title": f"{stacks.name} by {rows.name}",
        "param": [rows.field, stacks.field],
        "axis": bar_axes(stacked=True),
        "color_legend": {"display": True, "pos": [8, 8]},
        "position": [0, 0],
        "gsposition": position,
        **HALF,
    }


def build_overview_view(roles: Roles, title: str, abstract: str, datasource: str = CELLS) -> dict[str, Any]:
    """Return a view with the same chart set and grid as the myeloid Cell atlas."""
    charts: list[dict[str, Any]] = [
        {
            "type": "text_box_chart",
            "title": "",
            "legend": "",
            "param": [],
            "text": overview_text(title, abstract),
            "gsposition": [0, 0],
            **WIDE,
        },
        umap_chart(roles, roles.cell_type, THIRD, [0, 2]),
        umap_chart(roles, roles.cell_group, THIRD, [4, 2]),
        umap_chart(roles, roles.tissue, THIRD, [8, 2]),
        umap_chart(roles, roles.disease, HALF, [0, 7]),
        stacked_chart(roles.disease, roles.cell_group, [6, 7]),
        stacked_chart(roles.tissue, roles.cell_type, [0, 11]),
        {
            "type": "row_chart",
            "title": roles.cell_type.name,
            "param": [roles.cell_type.field],
            "axis": bar_axes(stacked=False),
            "position": [0, 0],
            "gsposition": [6, 11],
            **HALF,
        },
    ]
    used: set[str] = set()
    for chart in charts:
        chart["id"] = chart_id(used)
    return {
        "dataSources": {datasource: {"layout": "gridstack", "panelWidth": 100}},
        "initialCharts": {datasource: charts},
    }


def order_views(all_views: list[str], overview: str) -> list[str]:
    """Overview first, question views in question order, then the rest as they were."""
    listed = [name for name in dict.fromkeys(all_views) if name != overview]
    asked = [name for name in QUESTION_VIEW_ORDER if name in listed]
    rest = [name for name in listed if name not in asked]
    return [overview, *asked, *rest]


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def backup_file(path: Path, timestamp: str) -> Path:
    backup = path.with_name(f"{path.name}.bak.{timestamp}")
    backup.write_bytes(path.read_bytes())
    return backup


def add_overview_view(
    project: Path,
    *,
    title: str | None,
    abstract: str | None,
    overrides: dict[str, str | None] | None = None,
    view_name: str = VIEW_NAME,
    datasource: str = CELLS,
    replace: bool = False,
    write: bool = False,
) -> dict[str, Any]:
    """Plan the overview view and view order, and write them when ``write`` is true."""
    project = project.expanduser().resolve()
    paths = {name: project / f"{name}.json" for name in ("datasources", "views", "state")}
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(f"Not a complete MDV project; missing {path}")
    if write and not (title and abstract):
        raise ValueError("--write needs --title and --abstract-file")

    datasources = read_json(paths["datasources"])
    views: dict[str, Any] = read_json(paths["views"])
    state: dict[str, Any] = read_json(paths["state"])
    cells = next((item for item in datasources if item.get("name") == datasource), None)
    if cells is None:
        raise ValueError(f"Project has no {datasource!r} datasource")
    if view_name in views and not replace:
        raise ValueError(f"View {view_name!r} already exists; pass --replace to rebuild it")

    roles = choose_roles(cells, overrides or {})
    view = build_overview_view(
        roles, title or PLACEHOLDER_TITLE, abstract or PLACEHOLDER_ABSTRACT, datasource
    )
    current = [name for name in state.get("all_views") or [] if isinstance(name, str)]
    missing = [name for name in current if name not in views and name != view_name]
    new_order = order_views([name for name in current if name not in missing], view_name)

    summary: dict[str, Any] = {
        "project": str(project),
        "view_name": view_name,
        "replaced": view_name in views,
        "columns": {
            role: {"field": column.field, "name": column.name}
            for role, column in ((item.name, getattr(roles, item.name)) for item in fields(roles))
        },
        "charts": [
            {"type": chart["type"], "title": chart["title"], "param": chart["param"], "color_by": chart.get("color_by")}
            for chart in view["initialCharts"][datasource]
        ],
        "old_order": current,
        "new_order": new_order,
        "dropped_missing_views": missing,
        "initial_view": {"old": state.get("initial_view"), "new": view_name},
        "show_gallery_on_open": state.get("show_gallery_on_open"),
        "written": False,
    }
    if not write:
        return summary

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    summary["backups"] = [str(backup_file(paths[name], timestamp)) for name in ("views", "state")]
    views[view_name] = view
    atomic_write_json(paths["views"], views)
    state["all_views"] = new_order
    state["initial_view"] = view_name
    atomic_write_json(paths["state"], state)
    summary["written"] = True
    return summary


def print_summary(summary: dict[str, Any], *, has_text: bool) -> None:
    print(f"Project: {summary['project']}")
    action = "Rebuild" if summary["replaced"] else "Add"
    print(f"{action} view: {summary['view_name']}")
    print("\nColumns")
    for role, column in summary["columns"].items():
        print(f"  {role:10} {column['field']}  ({column['name']})")
    print("\nCharts")
    for chart in summary["charts"]:
        colour = f"  colour: {chart['color_by']}" if chart["color_by"] else ""
        print(f"  {chart['type']:18} {chart['title']!r:40} {chart['param']}{colour}")
    print("\nView order")
    for index, name in enumerate(summary["new_order"], start=1):
        print(f"  {index:2}. {name}")
    if summary["dropped_missing_views"]:
        print("  dropped (in all_views but not in views.json):")
        for name in summary["dropped_missing_views"]:
            print(f"    - {name}")
    initial = summary["initial_view"]
    print(f"\ninitial_view: {initial['old']!r} -> {initial['new']!r}")
    if summary["show_gallery_on_open"]:
        print(
            "Note: show_gallery_on_open is true, so the project opens on the view gallery "
            "rather than on the initial view. This script leaves it unchanged."
        )
    if not has_text:
        print("Note: no --title/--abstract-file given; the plan uses placeholder text.")
    if summary["written"]:
        print("\nBackups:")
        for path in summary["backups"]:
            print(f"  {path}")
        print("Write complete. datafile.h5 was not opened.")
    else:
        print("\nDry run. Nothing was written. Pass --write to apply.")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project", type=Path, required=True, help="MDV project directory.")
    parser.add_argument("--title", help="Project title, shown as the heading of the text box.")
    parser.add_argument("--abstract-file", type=Path, help="Text file holding the abstract (Markdown allowed).")
    parser.add_argument("--view-name", default=VIEW_NAME, help=f"Name of the overview view (default: {VIEW_NAME}).")
    parser.add_argument("--datasource", default=CELLS, help=f"Cells datasource name (default: {CELLS}).")
    for role in ("cell_type", "cell_group", "disease", "tissue"):
        parser.add_argument(
            f"--{role.replace('_', '-')}",
            dest=role,
            help=f"Field id or display name for {role.replace('_', ' ')} (default: first of {', '.join(ROLE_NAMES[role])}).",
        )
    parser.add_argument("--umap-x", dest="umap_x", help="Field id or display name of the UMAP x column.")
    parser.add_argument("--umap-y", dest="umap_y", help="Field id or display name of the UMAP y column.")
    parser.add_argument("--replace", action="store_true", help="Rebuild the view when it already exists.")
    parser.add_argument("--write", action="store_true", help="Back up views.json and state.json, then write them.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    abstract = None
    if args.abstract_file is not None:
        abstract = args.abstract_file.expanduser().read_text(encoding="utf-8")
    overrides = {
        key: getattr(args, key)
        for key in ("cell_type", "cell_group", "disease", "tissue", "umap_x", "umap_y")
    }
    try:
        summary = add_overview_view(
            args.project,
            title=args.title,
            abstract=abstract,
            overrides=overrides,
            view_name=args.view_name,
            datasource=args.datasource,
            replace=args.replace,
            write=args.write,
        )
    except (ValueError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print_summary(summary, has_text=bool(args.title and abstract))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
