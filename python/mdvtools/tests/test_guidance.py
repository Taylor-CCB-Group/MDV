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
    assert "## Why this visualization" in md
    assert "## What to look for" in md
    assert "## Suggested next steps" in md
    assert "qc_runs" in md or "follow-up" in md.lower()


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
    guidance = "## Why this visualization\n\nTest summary.\n"
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
