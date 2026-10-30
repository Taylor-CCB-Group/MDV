"""Overview view and view order written by create_overview_view."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from mdvtools.template_views import create_question_views as questions
from mdvtools.template_views.create_overview_view import (
    VIEW_NAME,
    add_overview_view,
    main,
    order_views,
)


def column(field: str, name: str, datatype: str = "text", deleted: bool = False) -> dict[str, Any]:
    item: dict[str, Any] = {"field": field, "name": name, "datatype": datatype}
    if deleted:
        item["deleted"] = True
    return item


def write_project(root: Path, *, columns: list[dict[str, Any]] | None = None, state: dict[str, Any] | None = None) -> Path:
    project = root / "project"
    project.mkdir()
    cells = columns or [
        column("rna:X_umap_mindist_0.25_1", "UMAP 1", "double"),
        column("rna:X_umap_mindist_0.25_2", "UMAP 2", "double"),
        column("rna:X_umap_1", "rna:X_umap_1", "double", deleted=True),
        column("rna:annotations", "Cell Type"),
        column("rna:bucket", "Cell Group"),
        column("rna:diagnosis", "Diagnosis"),
        column("rna:site", "Biopsy Site"),
    ]
    datasources = [{"name": "cells", "size": 10, "columns": cells}]
    views = {"Summary": {"initialCharts": {}}, "Stroma": {"initialCharts": {}}, questions.Q1_VIEW_NAME: {"initialCharts": {}}}
    project_state = state or {
        "all_views": ["Summary", questions.Q1_VIEW_NAME, "Stroma"],
        "permission": "view",
        "name": "Synovial",
    }
    (project / "datasources.json").write_text(json.dumps(datasources), encoding="utf-8")
    (project / "views.json").write_text(json.dumps(views), encoding="utf-8")
    (project / "state.json").write_text(json.dumps(project_state), encoding="utf-8")
    return project


def read(project: Path, name: str) -> Any:
    return json.loads((project / f"{name}.json").read_text(encoding="utf-8"))


def test_order_puts_overview_then_questions_then_the_rest():
    order = order_views(
        ["Summary", questions.Q2_VIEW_NAME, "metadata", questions.Q1_VIEW_NAME, VIEW_NAME],
        VIEW_NAME,
    )
    assert order == [VIEW_NAME, questions.Q1_VIEW_NAME, questions.Q2_VIEW_NAME, "Summary", "metadata"]


def test_dry_run_plans_columns_and_leaves_files_alone(tmp_path: Path):
    project = write_project(tmp_path)
    before = {name: (project / f"{name}.json").read_bytes() for name in ("views", "state")}

    summary = add_overview_view(project, title=None, abstract=None)

    assert summary["written"] is False
    assert summary["columns"]["cell_type"]["field"] == "rna:annotations"
    assert summary["columns"]["tissue"]["field"] == "rna:site"
    assert summary["columns"]["umap_x"]["field"] == "rna:X_umap_mindist_0.25_1"
    assert summary["new_order"] == [VIEW_NAME, questions.Q1_VIEW_NAME, "Summary", "Stroma"]
    assert {name: (project / f"{name}.json").read_bytes() for name in ("views", "state")} == before


def test_write_adds_view_sets_order_and_keeps_other_state(tmp_path: Path):
    project = write_project(tmp_path)

    summary = add_overview_view(project, title="Synovial atlas", abstract="About the cells.", write=True)

    view = read(project, "views")[VIEW_NAME]
    charts = view["initialCharts"]["cells"]
    assert view["dataSources"] == {"cells": {"panelWidth": 100}}
    assert len(charts) == 8
    assert max(chart["position"][0] + chart["size"][0] for chart in charts) == 1837
    umap_row = [chart for chart in charts if chart["position"][1] == charts[1]["position"][1]]
    gaps = [b["position"][0] - (a["position"][0] + a["size"][0]) for a, b in zip(umap_row, umap_row[1:])]
    assert len(gaps) == 2 and all(23 <= gap <= 25 for gap in gaps)
    assert charts[1]["position"][1] - charts[0]["size"][1] == 25
    assert len({chart["id"] for chart in charts}) == 8
    assert charts[0]["text"] == "# Synovial atlas\nAbout the cells."
    used = {field for chart in charts for field in chart["param"]}
    assert used == {"rna:X_umap_mindist_0.25_1", "rna:X_umap_mindist_0.25_2", "rna:annotations", "rna:bucket", "rna:diagnosis", "rna:site"}
    state = read(project, "state")
    assert state["all_views"] == [VIEW_NAME, questions.Q1_VIEW_NAME, "Summary", "Stroma"]
    assert state["initial_view"] == VIEW_NAME
    assert state["permission"] == "view"
    assert state["name"] == "Synovial"
    assert all(Path(path).is_file() for path in summary["backups"])


def test_existing_view_needs_replace(tmp_path: Path):
    project = write_project(tmp_path)
    add_overview_view(project, title="T", abstract="A", write=True)

    with pytest.raises(ValueError, match="--replace"):
        add_overview_view(project, title="T", abstract="A", write=True)

    add_overview_view(project, title="New", abstract="A", replace=True, write=True)
    assert read(project, "views")[VIEW_NAME]["initialCharts"]["cells"][0]["text"].startswith("# New")
    assert read(project, "state")["all_views"].count(VIEW_NAME) == 1


def test_write_needs_title_and_abstract(tmp_path: Path):
    project = write_project(tmp_path)
    with pytest.raises(ValueError, match="--title"):
        add_overview_view(project, title="T", abstract=None, write=True)


def test_deleted_columns_are_not_chosen(tmp_path: Path):
    project = write_project(
        tmp_path,
        columns=[
            column("x", "UMAP 1", "double"),
            column("y", "UMAP 2", "double"),
            column("old_type", "Cell Type", deleted=True),
            column("bucket", "Cell Group"),
            column("diagnosis", "Diagnosis"),
            column("site", "Tissue"),
        ],
    )
    with pytest.raises(ValueError, match="cell_type"):
        add_overview_view(project, title=None, abstract=None)


def test_override_picks_a_column_by_field_or_display_name(tmp_path: Path):
    project = write_project(
        tmp_path,
        columns=[
            column("x", "UMAP 1", "double"),
            column("y", "UMAP 2", "double"),
            column("leiden", "Leiden Clusters"),
            column("bucket", "Cell Group"),
            column("diagnosis", "Diagnosis"),
            column("site", "Tissue"),
        ],
    )
    summary = add_overview_view(project, title=None, abstract=None, overrides={"cell_type": "Leiden Clusters"})
    assert summary["columns"]["cell_type"]["field"] == "leiden"


def test_cli_reports_errors_without_writing(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    project = write_project(tmp_path)
    assert main(["--project", str(project), "--write"]) == 1
    assert "--title" in capsys.readouterr().err
    assert VIEW_NAME not in read(project, "views")
