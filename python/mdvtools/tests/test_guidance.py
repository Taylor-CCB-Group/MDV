from __future__ import annotations

import json

from mdvtools.llm.guidance import (
    ANALYSIS_SUMMARY_TITLE,
    CHATMDV_SUMMARY_TEXTBOX_ID,
    append_guidance_textbox_to_view,
    build_response_guidance,
)


class FakeProject:
    def __init__(self):
        self.datasources = [
            {"name": "qc_field_uniformity"},
            {"name": "qc_runs"},
            {"name": "qc_channels"},
        ]
        self._views: dict[str, dict] = {}
        self._columns = {
            "channel_name": {"field": "channel_name", "name": "channel_name"},
            "cv_pct": {"field": "cv_pct", "name": "cv_pct"},
            "X_pca_1": {"field": "X_pca_1", "name": "X_pca_1"},
            "X_pca_2": {"field": "X_pca_2", "name": "X_pca_2"},
            "leiden": {"field": "leiden", "name": "leiden"},
        }

    def get_view(self, view_name: str):
        return self._views[view_name]

    def set_view(self, view_name: str, view: dict):
        self._views[view_name] = view

    def get_column_metadata(self, datasource_name: str, field: str):
        return self._columns.get(field, {"field": field, "name": field})

    def get_datasource_metadata(self, name: str):
        return {"name": name, "size": 11, "columns": list(self._columns.values())}


def test_build_response_guidance_includes_sections():
    proj = FakeProject()
    proj._views["v1"] = {
        "initialCharts": {
            "qc_field_uniformity": [
                {
                    "title": "cv_pct by channel",
                    "type": "box_plot",
                    "param": ["channel_name", "cv_pct"],
                }
            ]
        }
    }
    code = "datasource_name = 'qc_field_uniformity'\n"
    md = build_response_guidance(
        proj,
        "What is the distribution of cv_pct across channel_name in qc_field_uniformity?",
        code,
        "v1",
        'fields "channel_name", "cv_pct"\ncharts "Box plot"',
        "| channel | cv_pct |\n|---|---|\n| C405 | 12 |",
    )
    assert "## 1. Why this chart is the best way to answer the question" in md
    assert "## 2. Biological insights that can be gained" in md
    assert "## 3. Subsequent analysis tasks" in md
    assert "### In summary" in md
    assert "The user asked:" in md
    assert "**Box Plot**" in md
    assert "qc_runs" in md or "Follow-up" in md


def test_build_response_guidance_subclustering_style_explanation():
    proj = FakeProject()
    proj.datasources = [{"name": "cells"}]
    proj._columns["ARVCF"] = {"field": "ARVCF", "name": "ARVCF"}
    gene_wrapper = "rna_expr|ARVCF(rna_expr)|42"
    proj._views["v1"] = {
        "initialCharts": {
            "cells": [
                {
                    "title": "PCA colored by ARVCF",
                    "type": "wgl_scatter_plot",
                    "param": ["X_pca_1", "X_pca_2"],
                    "color_by": gene_wrapper,
                },
                {
                    "title": "Density",
                    "type": "density_scatter_plot",
                    "param": ["X_pca_1", "X_pca_2", "leiden"],
                },
                {
                    "title": "Cells table",
                    "type": "table_chart",
                    "param": ["leiden", "X_pca_1", "X_pca_2", gene_wrapper],
                },
                {
                    "title": "Filter cells",
                    "type": "selection_dialog",
                    "param": ["leiden", "X_pca_1", "X_pca_2", gene_wrapper],
                },
            ]
        }
    }
    question = (
        "Could you perform subclustering within cluster 1 based on the "
        "expression levels of the gene ARVCF?"
    )
    code = (
        "datasource_name = 'cells'\n"
        "subset = adata.obs['leiden'] == '1'\n"
    )
    md = build_response_guidance(proj, question, code, "v1", "", None)
    assert "ARVCF" in md
    assert "**Scatter Plot (2D)**" in md
    assert "**Density Scatter Plot**" in md
    assert "**Table Plot**" in md
    assert "**Selection Dialog Plot**" in md
    assert "Subpopulation discovery" in md
    assert "Marker gene analysis" in md or "Automated subclustering" in md
    assert "cluster 1" in md.lower()


def test_append_guidance_textbox_prepends_and_replaces_summary():
    proj = FakeProject()
    proj._views["v1"] = {
        "initialCharts": {
            "qc_field_uniformity": [
                {
                    "title": ANALYSIS_SUMMARY_TITLE,
                    "type": "text_box_chart",
                    "id": CHATMDV_SUMMARY_TEXTBOX_ID,
                    "text": "old summary",
                },
                {
                    "title": "cv_pct by channel",
                    "type": "box_plot",
                    "param": ["channel_name", "cv_pct"],
                },
            ]
        }
    }
    guidance = "## 1. Why this chart is the best way to answer the question\n\nTest summary.\n"
    ok = append_guidance_textbox_to_view(
        proj, "v1", "qc_field_uniformity", guidance
    )
    assert ok is True
    charts = proj._views["v1"]["initialCharts"]["qc_field_uniformity"]
    assert len(charts) == 2
    assert charts[0]["type"] == "text_box_chart"
    assert charts[0]["id"] == CHATMDV_SUMMARY_TEXTBOX_ID
    assert charts[0]["text"] == guidance
    assert charts[1]["type"] == "box_plot"


def test_append_guidance_skips_when_no_plot_charts():
    proj = FakeProject()
    proj._views["empty"] = {"initialCharts": {"qc_field_uniformity": []}}
    ok = append_guidance_textbox_to_view(
        proj, "empty", "qc_field_uniformity", "summary"
    )
    assert ok is False
