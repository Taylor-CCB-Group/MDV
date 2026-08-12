"""
Post-execution analysis summary for ChatMDV.

Builds user-facing guidance (why / biological insights / next steps), persists it as a
view TextBox, and echoes the same markdown to the chat client.
"""
from __future__ import annotations

import json
import re
from typing import Any, Optional

from mdvtools.charts.text_box_plot import TextBox
from mdvtools.llm.verification import (
    CHART_PARAM_ROLE_MAP,
    SavedViewChart,
    _WRAPPER_RE,
    collect_saved_view_charts,
    format_chart_param_display,
    parse_datasource_name,
)

ANALYSIS_SUMMARY_TITLE = "Analysis summary"
CHATMDV_SUMMARY_TEXTBOX_ID = "chatmdv_analysis_summary"

_CHART_TYPE_LABELS: dict[str, str] = {
    "wgl_scatter_plot": "Scatter Plot (2D)",
    "density_scatter_plot": "Density Scatter Plot",
    "box_plot": "Box Plot",
    "violin_plot": "Violin Plot",
    "selection_dialog": "Selection Dialog Plot",
    "table_chart": "Table Plot",
    "single_heat_map": "Heatmap",
    "dot_plot": "Dot Plot",
    "histogram_chart": "Histogram",
    "multi_line_chart": "Multiline Chart",
    "pie_chart": "Pie Chart",
    "ring_chart": "Ring Chart",
    "stacked_row_chart": "Stacked Row Chart",
    "row_chart": "Row Chart",
    "sankey_chart": "Sankey Plot",
    "wordcloud_chart": "Wordcloud",
    "abundance_box_plot": "Abundance Box Plot",
}


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


def _short_field_label(project: Any, datasource_name: str, token: str) -> str:
    """Compact label for narrative text (gene name, display name, or field id)."""
    m = _WRAPPER_RE.match(token.strip())
    if m:
        return m.group(2).strip()
    try:
        col = project.get_column_metadata(datasource_name, token)
        display = col.get("name") or token
        return str(display)
    except Exception:
        return token


def _chart_type_label(chart_type: str) -> str:
    return _CHART_TYPE_LABELS.get(
        chart_type, chart_type.replace("_", " ").title()
    )


def _role_params(
    project: Any, chart: SavedViewChart
) -> list[tuple[str, str, str]]:
    """Return (role, raw_token, display) for each chart param."""
    roles = CHART_PARAM_ROLE_MAP.get(chart.chart_type)
    out: list[tuple[str, str, str]] = []
    for idx, token in enumerate(chart.param_tokens):
        role = roles[idx] if roles and idx < len(roles) else f"Param {idx}"
        display = format_chart_param_display(project, chart.datasource_name, token)
        out.append((role, token, display))
    return out


def _detect_subset_context(question: str, code: str) -> str | None:
    """Best-effort description of a row subset (e.g. a cluster) from question/code."""
    for text in (question, code):
        m = re.search(
            r"(?:within|inside|in)\s+(?:cluster|group)\s+['\"]?(\w+)['\"]?",
            text,
            re.IGNORECASE,
        )
        if m:
            return f"cluster {m.group(1)}"
        m = re.search(r"cluster\s+['\"]?(\w+)['\"]?", text, re.IGNORECASE)
        if m:
            return f"cluster {m.group(1)}"
    if ".loc[" in code or ".query(" in code.lower() or "adata[" in code:
        return "the filtered cell subset"
    return None


def _gene_tokens_from_charts(
    project: Any, charts: list[SavedViewChart]
) -> list[str]:
    genes: list[str] = []
    seen: set[str] = set()
    for chart in charts:
        for token in (*chart.param_tokens, chart.color_by or ""):
            m = _WRAPPER_RE.match(token.strip()) if token else None
            if m:
                gene = m.group(2).strip()
                if gene not in seen:
                    seen.add(gene)
                    genes.append(gene)
            elif token:
                label = _short_field_label(project, chart.datasource_name, token)
                if re.search(r"gene|expr", token, re.IGNORECASE) and label not in seen:
                    seen.add(label)
                    genes.append(label)
    return genes


def _embedding_axes(charts: list[SavedViewChart]) -> list[str]:
    axes: list[str] = []
    seen: set[str] = set()
    for chart in charts:
        for token in chart.param_tokens[:2]:
            low = token.lower()
            if "umap" in low or "pca" in low or "tsne" in low:
                if token not in seen:
                    seen.add(token)
                    axes.append(token)
    return axes


def _detect_intents(
    question: str, code: str, charts: list[SavedViewChart]
) -> set[str]:
    q = question.lower()
    intents: set[str] = set()
    if any(k in q for k in ("subcluster", "sub-cluster", "within cluster", "inside cluster")):
        intents.add("subcluster")
    if any(k in q for k in ("expression", "express", "gene", "marker", "deg", "de ")):
        intents.add("expression")
    if any(k in q for k in ("marker", "rank_genes", "differential", "deg", "de gene")):
        intents.add("markers")
    if any(k in q for k in ("proportion", "abundance", "composition", "frequency", "percent")):
        intents.add("proportion")
    if any(k in q for k in ("qc", "quality", "uniformity", "cv", "channel")):
        intents.add("qc")
    if any(k in q for k in ("correlat", "relationship", "association", "compare", "versus", " vs ")):
        intents.add("relationship")
    if any(k in q for k in ("trajectory", "pseudotime", "lineage", "development")):
        intents.add("trajectory")
    if any(c.chart_type == "selection_dialog" for c in charts):
        intents.add("interactive")
    if _embedding_axes(charts):
        intents.add("embedding")
    gene_list = []
    for c in charts:
        for t in c.param_tokens:
            m = _WRAPPER_RE.match(t.strip())
            if m:
                gene_list.append(m.group(2))
    if gene_list:
        intents.add("expression")
    if not intents:
        intents.add("explore")
    return intents


def _describe_chart_rationale(
    project: Any,
    chart: SavedViewChart,
    subset: str | None,
) -> str:
    label = _chart_type_label(chart.chart_type)
    roles = _role_params(project, chart)
    subset_phrase = f" for cells in {subset}" if subset else ""

    if chart.chart_type == "wgl_scatter_plot":
        x = roles[0][2] if roles else "X"
        y = roles[1][2] if len(roles) > 1 else "Y"
        color = (
            format_chart_param_display(project, chart.datasource_name, chart.color_by)
            if chart.color_by
            else None
        )
        lines = [
            f"Plots {x} vs {y}{subset_phrase}"
            + (f", colored by {color}." if color else ".")
        ]
        if color and "Expression feature" in color:
            lines.append(
                "This allows visual identification of subclusters or gradients "
                "that may correspond to different expression levels."
            )
        elif color:
            lines.append(
                "This reveals how cells group in embedding or feature space and "
                "whether categories or expression levels form distinct regions."
            )
        else:
            lines.append(
                "This shows spatial structure in the chosen dimensions for the selected cells."
            )
        return "\n".join(lines)

    if chart.chart_type == "density_scatter_plot":
        x = roles[0][2] if roles else "X"
        y = roles[1][2] if len(roles) > 1 else "Y"
        cat = roles[2][2] if len(roles) > 2 else None
        base = f"Adds density information to a scatter of {x} vs {y}{subset_phrase}"
        if cat:
            base += f", grouped or colored by {cat}"
        base += ", highlighting regions of high cell density."
        return (
            base
            + " Dense regions may indicate stable subpopulations; sparse tails may be transitional or rare states."
        )

    if chart.chart_type == "selection_dialog":
        filters = ", ".join(display for _, _, display in roles) or "(see chart)"
        return (
            f"Allows interactive subsetting and filtering by {filters}, supporting "
            "further subclustering, focused views, or downstream analysis on selected cells."
        )

    if chart.chart_type == "table_chart":
        cols = ", ".join(display for _, _, display in roles) or "relevant columns"
        return (
            f"Provides a tabular view of {cols}{subset_phrase}, enabling detailed "
            "inspection, sorting, and row-level review of the underlying data."
        )

    if chart.chart_type in ("box_plot", "violin_plot"):
        cat = roles[0][2] if roles else "category"
        val = roles[1][2] if len(roles) > 1 else "value"
        kind = "distribution spread" if chart.chart_type == "violin_plot" else "median and spread"
        return (
            f"Compares {val} across {cat}{subset_phrase}, showing {kind} per group "
            "to identify shifts, outliers, or group-specific differences."
        )

    if chart.chart_type == "single_heat_map":
        return (
            "Summarizes expression or scores across categories and features in a matrix view, "
            "making patterns of high and low values easy to compare at a glance."
        )

    if chart.chart_type == "dot_plot":
        return (
            "Shows per-category expression or scores for multiple features simultaneously, "
            "useful for comparing marker patterns across groups."
        )

    if chart.chart_type in ("pie_chart", "ring_chart", "stacked_row_chart"):
        return (
            "Shows relative proportions or composition across categories, "
            "helping assess whether group sizes or abundances differ meaningfully."
        )

    # Generic fallback
    param_summary = ", ".join(display for _, _, display in roles)
    if param_summary:
        return (
            f"Visualizes {param_summary}{subset_phrase} using a {label.lower()} "
            "to support inspection of the patterns relevant to your question."
        )
    return f"Provides a {label.lower()} view{subset_phrase} to help answer your question."


def _synthesis_sentence(
    question: str,
    charts: list[SavedViewChart],
    genes: list[str],
    subset: str | None,
) -> str:
    parts: list[str] = []
    if subset and genes:
        parts.append(
            f"This combination directly addresses your question by visualizing heterogeneity "
            f"within {subset} as a function of {', '.join(genes)} expression"
        )
    elif genes:
        parts.append(
            f"Together, these charts link {', '.join(genes)} expression to cell grouping "
            "and support both overview and detailed inspection"
        )
    elif subset:
        parts.append(
            f"These views focus on {subset} and combine overview plots with tools for "
            "interactive exploration where appropriate"
        )
    else:
        parts.append(
            "Together, these charts answer your question with complementary views—overview "
            "plots plus tabular or interactive inspection where saved"
        )
    if any(c.chart_type == "selection_dialog" for c in charts):
        parts.append("and enables hypothesis-driven downstream analysis on selected cells")
    else:
        parts.append("and helps validate patterns seen in the plots")
    return parts[0] + ", " + parts[1] + "."


def _biological_insights(
    intents: set[str],
    genes: list[str],
    subset: str | None,
    data_preview: str | None,
) -> list[str]:
    bullets: list[str] = []
    gene_phrase = ", ".join(genes) if genes else "the chosen feature"

    if "subcluster" in intents or ("expression" in intents and subset):
        bullets.append(
            "**Subpopulation discovery:** If expression is heterogeneous"
            + (f" within {subset}" if subset else "")
            + ", scatter and density plots may reveal spatially distinct subclusters, "
            "suggesting functional or developmental diversity."
        )
    if "expression" in intents or genes:
        bullets.append(
            f"**Gene expression patterns:** The visualization may show gradients or discrete "
            f"groups of {gene_phrase}-high and {gene_phrase}-low cells, which could correspond "
            "to different cell states or subtypes."
        )
    if "markers" in intents:
        bullets.append(
            "**Marker and DE context:** Patterns in expression or heatmap/dot views can highlight "
            "candidate markers and whether differences are cluster-wide or subset-specific."
        )
    if "proportion" in intents:
        bullets.append(
            "**Composition shifts:** Proportion charts can reveal whether certain groups are "
            "over- or under-represented, which may reflect biology or technical bias."
        )
    if "qc" in intents:
        bullets.append(
            "**Data quality:** Distribution and comparison plots help spot outlier channels, "
            "runs, or assays that may need filtering before biological interpretation."
        )
    if "relationship" in intents:
        bullets.append(
            "**Association structure:** Scatter or density views can reveal correlations, "
            "anti-correlations, or non-linear relationships between measured variables."
        )
    if "trajectory" in intents or ("expression" in intents and genes):
        bullets.append(
            "**Cluster purity and annotation:** If expression aligns with known markers or "
            "functional states, this can inform annotation or refinement of cluster labels."
        )
    if (data_preview or "").strip():
        bullets.append(
            "**Computed summaries:** Review the data preview in chat for printed aggregates, "
            "counts, or table excerpts that complement what you see in the charts."
        )
    if not bullets:
        bullets.append(
            "**Pattern discovery:** Look for separated groups, gradients, outliers, or "
            "category-specific shifts that motivate a more targeted follow-up analysis."
        )
    return bullets


def _suggest_next_steps(
    project: Any,
    question: str,
    datasource_name: str | None,
    intents: set[str],
    charts: list[SavedViewChart],
    genes: list[str],
    subset: str | None,
) -> list[str]:
    steps: list[str] = []
    gene_phrase = ", ".join(genes) if genes else "the feature of interest"
    has_selection = any(c.chart_type == "selection_dialog" for c in charts)

    if has_selection or "subcluster" in intents:
        steps.append(
            "**Automated subclustering:** Use the selection dialog to select expression-high "
            "or expression-low subpopulations"
            + (f" within {subset}" if subset else "")
            + ", then run further clustering or differential expression on the subset."
        )
    if genes or "expression" in intents or "markers" in intents:
        steps.append(
            f"**Marker gene analysis:** Identify other genes that co-vary with {gene_phrase}"
            + (f" within {subset}" if subset else "")
            + " to discover new markers or pathways."
        )
    if "trajectory" in intents or ("expression" in intents and genes):
        steps.append(
            f"**Trajectory / lineage inference:** If {gene_phrase} expression forms a gradient, "
            "investigate whether this reflects a developmental trajectory or continuous state change."
        )
    if genes and ("subcluster" in intents or has_selection):
        steps.append(
            f"**Functional enrichment:** Run GO or pathway enrichment on {gene_phrase}-high vs "
            f"{gene_phrase}-low subclusters to infer biological functions."
        )
    if "proportion" in intents:
        steps.append(
            "**Stratified comparisons:** Break down proportions by sample, batch, or condition "
            "to see whether composition differences are consistent across cohorts."
        )
    if "qc" in intents:
        steps.append(
            "**QC follow-up:** Filter or flag problematic channels/runs, then re-plot key metrics "
            "to confirm improvements before downstream analysis."
        )

    q_lower = question.lower()
    names: list[str] = []
    for ds in project.datasources or []:
        if isinstance(ds, dict) and ds.get("name"):
            names.append(str(ds["name"]))
    for name in names:
        if name == datasource_name or name.lower() in q_lower:
            continue
        steps.append(
            f"**Cross-datasource check:** Compare with related table `{name}` "
            "(e.g. run metadata, alternate QC tables, or companion assays)."
        )
        break

    if len(steps) < 3:
        steps.append(
            "**Follow-up question:** Ask ChatMDV to filter by a specific group, add another "
            "gene or numeric column, or export a table of cells matching a visual pattern."
        )
    return steps[:5]


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

    view = None
    charts: list[SavedViewChart] = []
    if view_name:
        try:
            view = project.get_view(view_name)
            if isinstance(view, dict):
                charts = collect_saved_view_charts(project, view)
        except Exception:
            view = None
            charts = []

    subset = _detect_subset_context(question, final_code)
    genes = _gene_tokens_from_charts(project, charts)
    intents = _detect_intents(question, final_code, charts)

    lines.append("## 1. Why this chart is the best way to answer the question")
    lines.append("")
    lines.append("The user asked:")
    lines.append(f'> "{question.strip()}"')
    lines.append("")

    if charts:
        for chart in charts:
            label = _chart_type_label(chart.chart_type)
            title = chart.title
            lines.append(f"**{label}**" + (f" ({title})" if title and title != label else "") + ":")
            lines.append(_describe_chart_rationale(project, chart, subset))
            lines.append("")
        lines.append(_synthesis_sentence(question, charts, genes, subset))
    else:
        lines.append(
            f"This analysis addresses your question using datasource "
            f"`{ds_from_code or 'see generated code'}`."
        )
        if agent_plan.strip():
            lines.append("")
            lines.append("Planned fields and charts:")
            lines.append(f"- {agent_plan.strip()[:600]}")

    lines.append("")
    lines.append("## 2. Biological insights that can be gained")
    lines.append("")
    for bullet in _biological_insights(intents, genes, subset, data_preview):
        lines.append(f"- {bullet}")

    preview = (data_preview or "").strip()
    if preview:
        lines.append("")
        snippet = preview[:1200]
        if len(preview) > 1200:
            snippet += "\n\n_(preview truncated)_"
        lines.append(snippet)

    lines.append("")
    lines.append("## 3. Subsequent analysis tasks")
    lines.append("")
    for step in _suggest_next_steps(
        project, question, ds_from_code, intents, charts, genes, subset
    ):
        lines.append(f"- {step}")

    lines.append("")
    lines.append("### In summary")
    if charts:
        lines.append(
            "This visualization suite enables visual exploration of the patterns in your question "
            "and provides tools for interactive, hypothesis-driven downstream analysis."
        )
    else:
        lines.append(
            "Review the computed outputs and ask a follow-up to add charts or refine the analysis."
        )

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
