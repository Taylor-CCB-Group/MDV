"""Legacy dot plot groups merged into compound dot plots, in question views and in views.json."""

from __future__ import annotations

import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import pytest

from mdvtools.template_views import create_question_views as myeloid
from mdvtools.template_views import create_synovial_question_views as synovial

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "convert_to_compound_dot_plots.py"
_spec = importlib.util.spec_from_file_location("convert_to_compound_dot_plots", _SCRIPT)
if _spec is None or _spec.loader is None:
    raise RuntimeError(f"could not load {_SCRIPT}")
convert = importlib.util.module_from_spec(_spec)
sys.modules["convert_to_compound_dot_plots"] = convert
_spec.loader.exec_module(convert)

GENES = {"type": "RowsAsColsQuery", "linkedDsName": "genes", "maxItems": 5, "subgroupName": "gs"}


class FieldIds(defaultdict[str, str]):
    """Feature field ids for the builders: each gene's id is its name."""

    def __missing__(self, key: str) -> str:
        return key


def charts_of(view: dict[str, Any]) -> list[dict[str, Any]]:
    return [chart for charts in view["initialCharts"].values() for chart in charts]


def filtered_dot_plots(view: dict[str, Any]) -> list[dict[str, Any]]:
    return [chart for chart in charts_of(view) if convert.split_column(chart) is not None]


def compound_plots(view: dict[str, Any]) -> list[dict[str, Any]]:
    return [chart for chart in charts_of(view) if chart["type"] == "compound_dot_plot"]


MYELOID_FIELDS = myeloid.ProjectFields(features=FieldIds(), protein_batches=["batch1"])


@pytest.mark.parametrize(
    ("view_name", "split_by", "split_values"),
    [
        (myeloid.Q6_VIEW_NAME, myeloid.DIAGNOSIS, myeloid.Q6_DIAGNOSES),
        (myeloid.Q7_VIEW_NAME, myeloid.DIAGNOSIS, myeloid.Q7_DIAGNOSES),
        (myeloid.Q8_VIEW_NAME, myeloid.RESPONSE_STAGE, myeloid.RESPONSE_STAGE_GROUPS),
        (myeloid.Q9_VIEW_NAME, myeloid.TREATMENT_STAGE, list(myeloid.DEG_GROUPS)),
    ],
)
def test_myeloid_comparison_views_use_one_compound_dot_plot(view_name, split_by, split_values):
    view = myeloid.VIEW_BUILDERS[view_name](MYELOID_FIELDS)
    [plot] = compound_plots(view)
    assert (plot["split_by"], plot["split_values"]) == (split_by, split_values)
    assert filtered_dot_plots(view) == []
    assert "Background Filter" not in charts_of(view)[0]["text"]


@pytest.mark.parametrize(
    ("view_name", "split_by", "split_values"),
    [
        (synovial.Q6_VIEW_NAME, synovial.DIAGNOSIS, synovial.Q6_DIAGNOSES),
        (synovial.Q7_VIEW_NAME, synovial.DIAGNOSIS, synovial.Q6_DIAGNOSES),
        (synovial.Q8_VIEW_NAME, synovial.TNF_RESPONSE, synovial.Q8_GROUPS),
        (synovial.Q9_VIEW_NAME, synovial.TNF_RESPONSE, synovial.Q8_GROUPS),
    ],
)
def test_synovial_comparison_views_use_one_compound_dot_plot(view_name, split_by, split_values):
    view = synovial.VIEW_BUILDERS[view_name](FieldIds())
    [plot] = compound_plots(view)
    assert (plot["split_by"], plot["split_values"]) == (split_by, split_values)
    assert filtered_dot_plots(view) == []


def test_patient_plot_keeps_its_background_filter():
    view = myeloid.VIEW_BUILDERS[myeloid.Q2_VIEW_NAME](MYELOID_FIELDS)
    [patients] = filtered_dot_plots(view)
    assert patients["background_filter"]["categories"] == myeloid.Q2_PATIENTS
    assert compound_plots(view) == []


# --- converter -------------------------------------------------------------------------------


def legacy_dot_plot(title: str, column: str | None, values: list[str], position: list[int], size: list[int]) -> dict[str, Any]:
    chart = myeloid.dot_plot(title, "cell_state", GENES, only=(column, values) if column else None)
    return {**chart, "position": position, "size": size, "id": title[:6]}


OLD_Q6_TEXT = (
    "**Q6.** The two **dot plots** show the selected genes in each cell state for RA and PsA: similar "
    "dots mean the same.\n\n**To change a dot plot's disease:** click its settings icon (⚙) and choose "
    "a diagnosis under **Background Filter Categories**."
)


PATIENT_PLOT = legacy_dot_plot("Selected genes by Patient", "patient", ["P1", "P2"], [0, 0], [500, 500])


def write_project(root: Path) -> Path:
    project = root / "project"
    project.mkdir()
    views = {
        "Q6": {
            "initialCharts": {
                "cells": [
                    {"type": "text_box_chart", "title": "About this page", "param": [], "text": OLD_Q6_TEXT, "id": "text"},
                    # listed right then left; the panels should follow the layout
                    legacy_dot_plot("Selected genes in PsA", "diagnosis", ["PsA"], [400, 300], [395, 355]),
                    legacy_dot_plot("Selected genes in RA", "diagnosis", ["RA"], [0, 300], [395, 355]),
                    legacy_dot_plot("Unfiltered", None, [], [0, 700], [400, 300]),
                ]
            }
        },
        "Q2": {
            "initialCharts": {
                "cells": [PATIENT_PLOT]
            }
        },
    }
    datasources = [{"name": "cells", "columns": [{"field": "diagnosis", "name": "Diagnosis"}]}]
    (project / "views.json").write_text(json.dumps(views), encoding="utf-8")
    (project / "datasources.json").write_text(json.dumps(datasources), encoding="utf-8")
    return project


def read_views(project: Path) -> dict[str, Any]:
    return json.loads((project / "views.json").read_text(encoding="utf-8"))


def test_converter_merges_a_group_and_updates_its_text(tmp_path: Path):
    project = write_project(tmp_path)

    summary = convert.convert_project(project, write=True)

    views = read_views(project)
    charts = charts_of(views["Q6"])
    assert [chart["type"] for chart in charts] == ["text_box_chart", "compound_dot_plot", "dot_plot"]
    plot = charts[1]
    assert plot["title"] == "Selected genes in RA and PsA"
    assert (plot["split_by"], plot["split_values"]) == ("diagnosis", ["RA", "PsA"])
    assert plot["param"] == ["cell_state", GENES]
    assert (plot["position"], plot["size"]) == ([0, 300], [795, 355])
    assert "background_filter" not in plot
    assert charts[0]["text"] == (
        "**Q6.** The **dot plot** shows the selected genes in each cell state, with a panel each for RA "
        "and PsA: similar dots mean the same.\n\n**To compare other diseases:** click the dot plot's "
        "settings icon (⚙) and choose diagnoses under **Split values**."
    )
    assert summary["changes"][0]["stale_text"] == []
    assert views["Q2"]["initialCharts"]["cells"] == [PATIENT_PLOT]
    assert Path(summary["backup"]).is_file()


def test_converter_dry_run_writes_nothing_and_second_run_has_nothing_to_do(tmp_path: Path):
    project = write_project(tmp_path)
    before = (project / "views.json").read_bytes()

    dry = convert.convert_project(project)
    assert dry["changes"] and not dry["written"]
    assert (project / "views.json").read_bytes() == before

    convert.convert_project(project, write=True)
    again = convert.convert_project(project, write=True)
    assert again["changes"] == [] and again["backup"] is None


def test_converter_keeps_plots_with_different_genes_apart(tmp_path: Path):
    project = write_project(tmp_path)
    views = read_views(project)
    views["Q6"]["initialCharts"]["cells"][2]["param"] = ["cell_state", "CD14"]
    (project / "views.json").write_text(json.dumps(views), encoding="utf-8")

    assert convert.convert_project(project)["changes"] == []


def test_title_falls_back_to_the_split_column_label():
    assert convert.merged_title(["R baseline", "NR baseline"], "Response Stage") == "Selected genes by Response Stage"
    assert convert.merged_title(["At a", "At b", "At c"], "x") == "At a, b and c"
