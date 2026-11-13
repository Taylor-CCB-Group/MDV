"""Question views follow the selected genes (active links) and keep legends off the plot."""

from __future__ import annotations

from typing import Any

from mdvtools.template_views import create_question_views as myeloid
from mdvtools.template_views import create_synovial_question_views as synovial
from mdvtools.template_views.create_question_views import place_legends

MYELOID_FIELDS = myeloid.ProjectFields(features={}, protein_batches=["batch1"])


def cells_charts(view: dict[str, Any]) -> list[dict[str, Any]]:
    return view["initialCharts"]["cells"]


def by_title(view: dict[str, Any], title: str, kind: str | None = None) -> dict[str, Any]:
    return next(
        chart
        for chart in cells_charts(view)
        if chart.get("title") == title and kind in (None, chart.get("type"))
    )


def is_active(value: Any, datasource: str, max_items: int) -> bool:
    return (
        isinstance(value, dict)
        and value.get("type") == "RowsAsColsQuery"
        and value.get("linkedDsName") == datasource
        and value.get("maxItems") == max_items
    )


def test_myeloid_gene_charts_follow_the_selected_genes():
    q1 = myeloid.build_q1_view(MYELOID_FIELDS)
    assert is_active(by_title(q1, "Selected Gene Expression")["color_by"], "genes", 1)
    assert is_active(by_title(q1, "Cell State", "heat_map")["param"][1], "genes", myeloid.HEAT_MAP_FEATURES)
    assert q1["initialCharts"]["genes"][0]["filters"]["name"]["category"] == [myeloid.Q1_START_GENE]

    q3 = myeloid.build_q3_view(MYELOID_FIELDS)
    assert is_active(by_title(q3, "Selected Protein")["color_by"], "proteins", 1)
    assert q3["initialCharts"]["proteins"][0]["filters"]["name"]["category"] == [myeloid.Q3_START_PAIR[1]]

    q6 = myeloid.build_q6_view(MYELOID_FIELDS)
    assert is_active(by_title(q6, "Gene Expression by Diagnosis")["param"][1], "genes", 1)


def test_synovial_q1_follows_the_selected_genes():
    q1 = synovial.build_q1_view({})
    assert is_active(by_title(q1, "Selected Gene Expression")["color_by"], "rna", 1)
    assert is_active(by_title(q1, "Cell Group", "heat_map")["param"][1], "rna", myeloid.HEAT_MAP_FEATURES)


def test_legends_go_bottom_right_and_keep_a_set_position():
    view = {
        "initialCharts": {
            "cells": [
                {"type": "wgl_scatter_plot_dev", "color_by": "group", "size": [600, 700], "color_legend": {"display": True}},
                {"type": "heat_map", "param": ["group"], "size": [700, 500], "color_legend": {"display": True}},
                {"type": "stacked_row_chart", "param": ["a", "group"], "size": [900, 550], "color_legend": {"display": True, "pos": [1, 2]}},
            ]
        }
    }
    place_legends(view, {"cells": {"group": 3}})
    scatter, heat, bar = cells_charts(view)
    width, height = scatter["color_legend"]["size"]
    assert scatter["color_legend"]["pos"] == [600 - width - 10, 700 - height - 75]
    assert heat["color_legend"]["pos"] == [8, 380]
    assert bar["color_legend"] == {"display": True, "pos": [1, 2]}
