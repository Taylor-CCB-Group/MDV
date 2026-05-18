from __future__ import annotations

import re

import pytest

from mdvtools.llm.verification import (
    build_verification_summary,
    format_charts_from_saved_view,
)


class FakeProject:
    def __init__(self, view: dict, column_names: dict[str, str], size: int = 100):
        self._view = view
        self._column_names = column_names
        self.datasources = [{"name": "cells"}]
        self._size = size

    def get_view(self, view_name: str):
        assert view_name == "v1"
        return self._view

    def get_column_metadata(self, datasource_name: str, field: str):
        assert datasource_name == "cells"
        # Simulate datasource column metadata with separate display names.
        return {"field": field, "name": self._column_names.get(field, field)}

    def get_datasource_metadata(self, name: str):
        assert name == "cells"
        return {"name": name, "size": self._size}


def test_format_charts_scatter_includes_roles_and_color_by():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "S1",
                    "type": "wgl_scatter_plot",
                    "param": ["x_field", "y_field"],
                    "color_by": "color_field",
                }
            ]
        }
    }
    proj = FakeProject(
        view=view,
        column_names={
            "x_field": "X Display",
            "y_field": "Y Display",
            "color_field": "Color Display",
        },
    )

    md, has_sel = format_charts_from_saved_view(proj, proj.get_view("v1"))
    assert has_sel is False
    assert "### Charts (from saved view)" in md
    assert "**S1** (`wgl_scatter_plot`)" in md
    assert "- X: `x_field` (X Display)" in md
    assert "- Y: `y_field` (Y Display)" in md
    assert "- Color: `color_field` (Color Display)" in md


def test_format_charts_density_scatter_roles():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "D1",
                    "type": "density_scatter_plot",
                    "param": ["x", "y", "cat"],
                }
            ]
        }
    }
    proj = FakeProject(
        view=view,
        column_names={"x": "X", "y": "Y", "cat": "Category Display"},
    )
    md, _ = format_charts_from_saved_view(proj, proj.get_view("v1"))
    assert "- X: `x` (X)" in md
    assert "- Y: `y` (Y)" in md
    assert "- Category: `cat` (Category Display)" in md


def test_format_charts_box_plot_roles():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "B1",
                    "type": "box_plot",
                    "param": ["cat_col", "value_col"],
                }
            ]
        }
    }
    proj = FakeProject(
        view=view,
        column_names={"cat_col": "Cell Type", "value_col": "Score"},
    )
    md, _ = format_charts_from_saved_view(proj, proj.get_view("v1"))
    assert "- Category: `cat_col` (Cell Type)" in md
    assert "- Value: `value_col` (Score)" in md


def test_format_charts_selection_dialog_columns_to_filter():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "F1",
                    "type": "selection_dialog",
                    "param": ["a", "b"],
                }
            ]
        }
    }
    proj = FakeProject(view=view, column_names={"a": "A", "b": "B"})
    md, has_sel = format_charts_from_saved_view(proj, proj.get_view("v1"))
    assert has_sel is True
    assert "**F1** (`selection_dialog`)" in md
    # Comma-separated list in the single line.
    assert "Columns To filter: `a` (A), `b` (B)" in md


def test_format_charts_unknown_type_falls_back_to_param_index():
    view = {
        "initialCharts": {
            "cells": [
                {"title": "U1", "type": "some_new_chart", "param": ["p0", "p1"]}
            ]
        }
    }
    proj = FakeProject(view=view, column_names={"p0": "P0", "p1": "P1"})
    md, has_sel = format_charts_from_saved_view(proj, proj.get_view("v1"))
    assert has_sel is False
    assert "- Param 0: `p0` (P0)" in md
    assert "- Param 1: `p1` (P1)" in md


def test_format_charts_skips_text_box_chart():
    view = {
        "initialCharts": {
            "cells": [
                {"title": "What you can verify", "type": "text_box_chart", "text": "X"}
            ]
        }
    }
    proj = FakeProject(view=view, column_names={})
    md, has_sel = format_charts_from_saved_view(proj, proj.get_view("v1"))
    # Because we skip all `text_box_chart` types, there should be no charts output.
    assert md == ""
    assert has_sel is False


def test_format_charts_gene_wrapper_renders_gene_label():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "G1",
                    "type": "wgl_scatter_plot",
                    "param": ["rna_expr|TP53(rna_expr)|12", "y_field"],
                }
            ]
        }
    }
    proj = FakeProject(view=view, column_names={"y_field": "Y"})
    md, _ = format_charts_from_saved_view(proj, proj.get_view("v1"))
    assert "- X: Expression feature: `TP53` (subgroup `rna_expr`)" in md
    assert "- Y: `y_field` (Y)" in md


def test_build_verification_summary_includes_selection_dialog_note():
    view = {
        "initialCharts": {
            "cells": [
                {"title": "F1", "type": "selection_dialog", "param": ["a"]}
            ]
        }
    }
    proj = FakeProject(view=view, column_names={"a": "A"}, size=123)
    md = build_verification_summary(proj, final_code="# no hints", view_name="v1")
    assert "### Charts (from saved view)" in md
    assert "- Interactive filtering enabled via selection dialog." in md
    assert "- **Total rows in datasource `cells` (project metadata):** 123" in md
    assert "Rows after code-level filtering" in md

