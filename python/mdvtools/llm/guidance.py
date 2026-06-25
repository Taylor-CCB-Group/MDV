"""
Post-execution analysis summary for ChatMDV.

Builds user-facing guidance (why / what to look for / next steps), persists it as a
view TextBox, and echoes the same markdown to the chat client.
"""
from __future__ import annotations

import json
from typing import Any, Optional

from mdvtools.charts.text_box_plot import TextBox
from mdvtools.llm.verification import (
    format_charts_from_saved_view,
    parse_datasource_name,
)

ANALYSIS_SUMMARY_TITLE = "Analysis summary"
CHATMDV_SUMMARY_TEXTBOX_ID = "chatmdv_analysis_summary"


def _is_summary_textbox(chart: dict[str, Any]) -> bool:
    if not isinstance(chart, dict):
        return False
    if chart.get("type") != "text_box_chart":
        return False
    if chart.get("id") == CHATMDV_SUMMARY_TEXTBOX_ID:
        return True
    title = chart.get("title")
    return isinstance(title, str) and title.strip() == ANALYSIS_SUMMARY_TITLE


def _plot_to_json(plot: TextBox) -> dict[str, Any]:
    return json.loads(json.dumps(plot.plot_data))


def _view_has_non_summary_charts(view: dict[str, Any], datasource_name: str) -> bool:
    initial = view.get("initialCharts")
    if not isinstance(initial, dict):
        return False
    charts = initial.get(datasource_name)
    if not isinstance(charts, list):
        return False
    return any(isinstance(c, dict) and not _is_summary_textbox(c) for c in charts)


def _suggest_next_steps(
    project: Any,
    question: str,
    datasource_name: str | None,
) -> list[str]:
    names: list[str] = []
    for ds in project.datasources or []:
        if isinstance(ds, dict) and ds.get("name"):
            names.append(str(ds["name"]))
    q_lower = question.lower()
    steps: list[str] = []
    for name in names:
        if name == datasource_name:
            continue
        if name.lower() in q_lower:
            continue
        steps.append(f"Compare with related table `{name}` (e.g. other QC metrics or run metadata).")
        if len(steps) >= 2:
            break
    steps.append(
        "Ask a follow-up to filter by session, channel, or assay—or to plot a related numeric column."
    )
    return steps[:3]


def build_response_guidance(
    project: Any,
    question: str,
    final_code: str,
    view_name: Optional[str],
    agent_plan: str,
    data_preview: str | None,
) -> str:
    """
    Markdown guidance echoed to chat and stored in the view TextBox.
    """
    ds_from_code = parse_datasource_name(final_code)
    lines: list[str] = []

    lines.append("## Why this visualization")
    view = None
    if view_name:
        try:
            view = project.get_view(view_name)
        except Exception:
            view = None
    charts_md = ""
    if isinstance(view, dict) and view.get("initialCharts"):
        try:
            charts_md, _ = format_charts_from_saved_view(project, view)
        except Exception:
            charts_md = ""
    if charts_md.strip():
        for part in charts_md.strip().splitlines():
            if part.startswith("### Charts"):
                lines.append(
                    "The saved view uses the chart(s) below to answer your question "
                    f"using datasource `{ds_from_code or 'see code'}`."
                )
            elif part.strip():
                lines.append(part)
    else:
        lines.append(
            f"This analysis addresses: **{question.strip()}**"
        )
        if agent_plan.strip():
            lines.append(f"- Planned fields/charts: {agent_plan.strip()[:400]}")

    lines.append("")
    lines.append("## What to look for")
    preview = (data_preview or "").strip()
    if preview:
        lines.append("- Review the **data preview** in chat for computed aggregates or row samples.")
        snippet = preview[:1200]
        if len(preview) > 1200:
            snippet += "\n\n_(preview truncated)_"
        lines.append("")
        lines.append(snippet)
    else:
        lines.append("- Compare groups or categories on the chart axes.")
        lines.append("- Note outliers, skew, or channels that differ strongly from others.")

    lines.append("")
    lines.append("## Suggested next steps")
    for step in _suggest_next_steps(project, question, ds_from_code):
        lines.append(f"- {step}")

    return "\n".join(lines).strip()


def append_guidance_textbox_to_view(
    project: Any,
    view_name: str,
    datasource_name: str,
    guidance_md: str,
) -> bool:
    """
    Prepend (or replace) the ChatMDV analysis summary TextBox on a saved view.

    Returns True when the view was updated.
    """
    if not guidance_md.strip():
        return False
    try:
        view = project.get_view(view_name)
    except Exception:
        return False
    if not isinstance(view, dict):
        return False
    if not _view_has_non_summary_charts(view, datasource_name):
        return False

    initial = view.get("initialCharts")
    if not isinstance(initial, dict):
        initial = {}
        view = {**view, "initialCharts": initial}

    charts = initial.get(datasource_name)
    if not isinstance(charts, list):
        charts = []

    kept = [c for c in charts if isinstance(c, dict) and not _is_summary_textbox(c)]

    summary = TextBox(
        title=ANALYSIS_SUMMARY_TITLE,
        param=[],
        size=[900, 320],
        position=[10, 10],
        id=CHATMDV_SUMMARY_TEXTBOX_ID,
    )
    summary.set_text(guidance_md)
    summary_json = _plot_to_json(summary)

    initial[datasource_name] = [summary_json, *kept]
    project.set_view(view_name, view)
    return True
