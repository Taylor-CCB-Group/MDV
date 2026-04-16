"""
Build a verification summary for Chat MDV from project metadata and generated code.
The summary is returned to the chat client (not embedded as a view TextBox chart).
"""
from __future__ import annotations

import ast
import re
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    from mdvtools.mdvproject import MDVProject

# Wrapper-based expression columns use this FieldName format:
#   "<subgroup>|<feature>(<subgroup>)|<index>"
# where subgroup is a rows-as-columns subgroup key.
_WRAPPER_RE = re.compile(r"([^|]+)\|([^|(]+)\(\1\)\|\s*(\d+)")
# Chart classes commonly instantiated in chat-generated scripts
_CHART_CLASS_RE = re.compile(
    r"\b(DotPlot|ScatterPlot|HeatmapPlot|HistogramPlot|BoxPlot|ViolinPlot|"
    r"DensityScatterPlot|ScatterPlot3D|RowChart|StackedRowChart|PieChart|RingChart|"
    r"AbundanceBoxPlot|MultiLinePlot|TablePlot|WordcloudPlot|SankeyPlot|"
    r"SelectionDialogPlot|RowSummaryBox|TextBox)\s*\("
)


def parse_datasource_name(code: str) -> Optional[str]:
    """Extract datasource_name = \"...\" from generated code."""
    try:
        tree = ast.parse(code)
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                if (
                    len(node.targets) == 1
                    and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id == "datasource_name"
                    and isinstance(node.value, ast.Constant)
                    and isinstance(node.value.value, str)
                ):
                    return node.value.value
    except SyntaxError:
        pass
    return None


def _project_datasource_names(project: Any) -> list[str]:
    return [ds["name"] for ds in project.datasources]


def _has_genes_table(project: Any) -> bool:
    for name in _project_datasource_names(project):
        if name == "genes" or "gene" in name.lower():
            return True
    return False


def _ensembl_style(token: str) -> bool:
    t = token.strip()
    if re.match(r"^ENS[A-Z]*G\d+", t):
        return True
    if re.match(r"^ENSMUSG\d+", t):
        return True
    if t.startswith("ENS") and len(t) > 6:
        return True
    return False


def _extract_chart_types(code: str) -> list[str]:
    seen: set[str] = set()
    order: list[str] = []
    for m in _CHART_CLASS_RE.finditer(code):
        name = m.group(1)
        if name not in seen and name != "TextBox":
            seen.add(name)
            order.append(name)
    return order


def _filter_hints(code: str) -> list[str]:
    """Heuristic: lines suggesting row subsetting in generated code."""
    hints: list[str] = []
    lowered = code.lower()
    if ".loc[" in code or ".iloc[" in code:
        hints.append("Uses pandas .loc / .iloc indexing (possible row subset).")
    if ".query(" in lowered:
        hints.append("Uses DataFrame.query(...) (possible row subset).")
    if re.search(r"adata\s*\[\s*[^,\]]+\s*,", code):
        hints.append("Uses AnnData subset indexing on adata[...].")
    if "boolean" in lowered and "[" in code and "adata" in lowered:
        pass
    return hints


def _to_param_list(p: Any) -> list[Any]:
    """Normalize `chart['param']` into a list for display."""
    if p is None:
        return []
    if isinstance(p, list):
        return p
    if isinstance(p, str):
        return [p]
    return [p]


def _format_param_value(project: Any, datasource_name: str, token: Any) -> str:
    """
    Format a saved-view param token for display.

    - Plain tokens are treated as datasource `field` ids, resolved through
      `project.get_column_metadata(datasource_name, field)`.
    - Wrapper tokens in the form `<subgroup>|<feature>(<subgroup>)|<index>` are linked-expression features.
    """
    if not isinstance(token, str):
        return str(token)

    s = token.strip()
    m = _WRAPPER_RE.match(s)
    if m:
        subgroup = m.group(1).strip()
        feature = m.group(2).strip()
        return f"Expression feature: `{feature}` (subgroup `{subgroup}`)"

    try:
        col = project.get_column_metadata(datasource_name, s)
        display = col.get("name") or s
        if str(display) != s:
            return f"`{s}` ({display})"
        return f"`{s}`"
    except Exception:
        # Token could be a virtual/linked column not present in datasource metadata.
        return f"`{s}`"


def format_charts_from_saved_view(
    project: Any, view: dict[str, Any]
) -> tuple[str, bool]:
    """
    Build a markdown subsection listing charts from the saved view JSON.

    Returns (markdown, has_selection_dialog).
    """
    initial = view.get("initialCharts")
    if not isinstance(initial, dict) or not initial:
        return "", False

    # Stable simplified role labels based on persisted `chart['type']`.
    role_map: dict[str, list[str]] = {
        "wgl_scatter_plot": ["X", "Y"],
        "density_scatter_plot": ["X", "Y", "Category"],
        "box_plot": ["Category", "Value"],
        "violin_plot": ["Category", "Value"],
    }

    has_selection_dialog = False
    out: list[str] = ["### Charts (from saved view)"]
    added_any = False

    for ds_name, charts in initial.items():
        if not isinstance(charts, list):
            continue
        for chart in charts:
            if not isinstance(chart, dict):
                continue
            ctype = chart.get("type")
            if not isinstance(ctype, str):
                continue

            # Skip any text-box type charts to avoid recursion/duplication.
            if ctype == "text_box_chart":
                continue

            # For robustness: some charts might not have a title.
            title = (chart.get("title") or "").strip() or ctype
            param_list = _to_param_list(chart.get("param"))

            if ctype == "selection_dialog":
                has_selection_dialog = True
                formatted = ", ".join(
                    _format_param_value(project, ds_name, p) for p in param_list
                )
                out.append(f"**{title}** (`{ctype}`)")
                out.append(
                    f"- Columns To filter: {formatted if formatted else '(none)'}"
                )
                # Also list color_by if a chart JSON includes it.
                if chart.get("color_by"):
                    out.append(
                        f"- Color: {_format_param_value(project, ds_name, chart.get('color_by'))}"
                    )
                out.append("")
                added_any = True
                continue

            out.append(f"**{title}** (`{ctype}`)")
            added_any = True
            roles = role_map.get(ctype)
            if roles:
                for idx, role in enumerate(roles):
                    if idx < len(param_list):
                        out.append(
                            f"- {role}: {_format_param_value(project, ds_name, param_list[idx])}"
                        )
                # Any remaining params beyond expected roles.
                for idx in range(len(roles), len(param_list)):
                    out.append(
                        f"- Param {idx}: {_format_param_value(project, ds_name, param_list[idx])}"
                    )
            else:
                for idx, p in enumerate(param_list):
                    out.append(f"- Param {idx}: {_format_param_value(project, ds_name, p)}")

            if chart.get("color_by"):
                out.append(
                    f"- Color: {_format_param_value(project, ds_name, chart.get('color_by'))}"
                )

            out.append("")

    if not added_any:
        return "", False

    return "\n".join(out).rstrip() + "\n", has_selection_dialog


def build_verification_summary(
    project: Any, final_code: str, view_name: Optional[str] = None
) -> str:
    """
    Plain-text (markdown-friendly) summary of what the user can verify from code and project data.
    """
    lines: list[str] = []
    lines.append("## What you can verify")
    lines.append("")

    ds_names = _project_datasource_names(project)
    lines.append(f"- **Datasources in project:** {', '.join(ds_names) if ds_names else '(none)'}")

    ds_from_code = parse_datasource_name(final_code)
    if ds_from_code:
        lines.append(f"- **Primary datasource used in generated code:** `{ds_from_code}`")

    view = None
    if view_name:
        try:
            view = project.get_view(view_name)
        except Exception:
            view = None

    charts_md = ""
    has_selection_dialog = False
    if isinstance(view, dict) and view.get("initialCharts"):
        try:
            charts_md, has_selection_dialog = format_charts_from_saved_view(
                project, view
            )
        except Exception:
            charts_md, has_selection_dialog = "", False

    if charts_md:
        lines.append("")
        lines.append(charts_md.rstrip())

    w_matches = list(_WRAPPER_RE.finditer(final_code))
    gene_tokens = [m.group(2).strip() for m in w_matches]

    if w_matches:
        lines.append("")
        lines.append("### Wrapper expression parameters")
        for tok in gene_tokens:
            kind = (
                "Ensembl-style ID (heuristic)"
                if _ensembl_style(tok)
                else "gene label as in `genes` / `var` `name` column (heuristic: not ENS* )"
            )
            lines.append(f"- `{tok}` — {kind}")
    else:
        lines.append("")
        lines.append(
            "- **Wrapper expression:** not used in this view (no `<subgroup>|…(<subgroup>)|…` parameters detected)."
        )
        if _has_genes_table(project) and not gene_tokens:
            lines.append(
                "- **Note:** this project includes a genes table; this view does not plot gene expression columns."
            )

    fh = _filter_hints(final_code)

    lines.append("")
    lines.append("### Rows / filtering")
    cell_ds = ds_from_code if ds_from_code and ds_from_code in ds_names else None
    if cell_ds is None:
        for candidate in ("cells", "obs", "cell"):
            for n in ds_names:
                if candidate in n.lower():
                    cell_ds = n
                    break
            if cell_ds:
                break
    if cell_ds is None and ds_names:
        cell_ds = ds_names[0]

    if cell_ds:
        try:
            md = project.get_datasource_metadata(cell_ds)
            n = md.get("size", "?")
            lines.append(
                f"- **Total rows in datasource `{cell_ds}` (project metadata):** {n}"
            )
            if fh:
                lines.append(
                    f"- **Rows after code-level filtering:** {n} (code may subset rows; this textbox uses project metadata after execution)"
                )
            else:
                lines.append(
                    f"- **Rows after code-level filtering:** {n} (no obvious code-level row subsetting detected)"
                )
        except Exception:
            lines.append(f"- Could not read size for datasource `{cell_ds}`.")
    else:
        lines.append("- Could not determine a cell-level datasource for row count.")

    if fh:
        for h in fh:
            lines.append(f"- Code hint: {h}")
    if has_selection_dialog:
        lines.append("- Interactive filtering enabled via selection dialog.")

    # Fallback to code-only details if we couldn't format saved-view charts.
    if not charts_md:
        lines.append("")
        lines.append("### Columns / parameters (from code)")
        if gene_tokens:
            lines.append(
                f"- Gene expression column(s): {', '.join(f'`{t}`' for t in gene_tokens)}"
            )

        param_lists = re.findall(r"params\s*=\s*\[([^\]]+)\]", final_code, re.DOTALL)
        cell_like: set[str] = set()
        for block in param_lists[:8]:
            for m in re.finditer(r'"([^"]+)"', block):
                s = m.group(1)
                if not s.startswith("gs|") and "|" not in s:
                    cell_like.add(s)
        if cell_like:
            sample = sorted(cell_like)[:40]
            extra = (
                f" (+{len(cell_like) - len(sample)} more)"
                if len(cell_like) > 40
                else ""
            )
            lines.append(
                f"- Other chart parameters (strings): {', '.join(f'`{c}`' for c in sample)}{extra}"
            )
        elif not gene_tokens:
            lines.append(
                "- (No `params=[...]` string list parsed; inspect the code block for column names.)"
            )

        lines.append("")
        lines.append("### Chart types (from code)")
        ctypes = _extract_chart_types(final_code)
        if ctypes:
            lines.append(f"- {', '.join(ctypes)}")
        else:
            lines.append("- (Could not infer chart class names from code.)")

    lines.append("")
    lines.append(
        "_Summary is generated from executed project metadata and saved view JSON; it may not reflect interactive UI-only filters._"
    )

    return "\n".join(lines)


