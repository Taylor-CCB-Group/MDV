from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Optional


@dataclass(frozen=True)
class RowsAsColumnsExpression:
    """
    A wrapper-capable expression datasource reachable via a rows-as-columns link.

    `subgroup_key` is the key used in the FieldName wrapper string:
        f"{subgroup_key}|{feature}({subgroup_key})|{row_index}"
    """

    datasource_name: str
    name_column: str
    subgroup_key: str
    subgroup_label: str


@dataclass(frozen=True)
class InferredDatasourceRoles:
    obs_datasource: str
    expressions: list[RowsAsColumnsExpression]

    def preferred_expression(self) -> Optional[RowsAsColumnsExpression]:
        if not self.expressions:
            return None
        # Prefer common modality names when present.
        for preferred in ("rna", "genes", "gene", "expression"):
            for e in self.expressions:
                if e.datasource_name.lower() == preferred:
                    return e
        return self.expressions[0]


def _pick_obs_datasource(project: Any) -> str:
    names = []
    try:
        names = list(project.get_datasource_names())
    except Exception:
        # Fallback: MDVProject.datasources property
        try:
            names = [ds.get("name") for ds in (project.datasources or []) if ds.get("name")]
        except Exception:
            names = []

    if "cells" in names:
        return "cells"
    # last-resort: first datasource
    if names:
        return str(names[0])
    raise ValueError("Project has no datasources; cannot infer roles")


def _pick_default_subgroup(subgroups: dict[str, Any]) -> tuple[str, str]:
    """
    Pick a reasonable default subgroup from rows-as-columns metadata.
    Prefer keys that look like primary expression (e.g. *_expr).
    """
    keys = list(subgroups.keys())
    if not keys:
        raise ValueError("rows_as_columns link has no subgroups")
    for k in keys:
        if k.lower().endswith("_expr") or k.lower() == "expr":
            sg = subgroups[k] or {}
            return k, str(sg.get("label") or sg.get("name") or k)
    k0 = keys[0]
    sg0 = subgroups[k0] or {}
    return k0, str(sg0.get("label") or sg0.get("name") or k0)


def infer_datasource_roles(project: Any) -> InferredDatasourceRoles:
    """
    Infer obs vs wrapper-expression datasources from MDV project metadata.

    This uses rows-as-columns links on the obs datasource to discover expression feature tables
    (e.g. `rna`, `protein`) and their available subgroups (e.g. `rna_expr`).
    """
    obs = _pick_obs_datasource(project)
    expressions: list[RowsAsColumnsExpression] = []

    # Discover rows-as-columns links from the obs datasource to other datasources.
    try:
        links = project.get_links(obs, "rows_as_columns")
    except Exception:
        links = []

    for item in links or []:
        try:
            target_ds = str(item.get("datasource"))
            link = item.get("link", {}).get("rows_as_columns", {})
            name_column = str(link.get("name_column") or "name")
            subgroups = link.get("subgroups") or {}
            if not isinstance(subgroups, dict) or len(subgroups) == 0:
                continue
            subgroup_key, subgroup_label = _pick_default_subgroup(subgroups)
            expressions.append(
                RowsAsColumnsExpression(
                    datasource_name=target_ds,
                    name_column=name_column,
                    subgroup_key=subgroup_key,
                    subgroup_label=subgroup_label,
                )
            )
        except Exception:
            # Best-effort inference; ignore malformed link entries.
            continue

    return InferredDatasourceRoles(obs_datasource=obs, expressions=expressions)


# MDV datasource column dict keys (see markdown_utils.create_column_markdown).
CATEGORICAL_DATATYPES = frozenset({"text", "text16", "multitext", "unique"})
NUMERIC_DATATYPES = frozenset({"integer", "double", "int32"})

CHATMDV_CATEGORICAL_FIELD_IDS_CAP = 30


def column_field_id(col: dict[str, Any]) -> str:
    """Return the chart param field id for a metadata column entry."""
    field = col.get("field")
    if field is not None and str(field).strip():
        return str(field)
    name = col.get("name")
    if name is not None and str(name).strip():
        return str(name)
    return ""


def categorical_field_ids_from_metadata(ds_meta: dict[str, Any]) -> list[str]:
    """Field ids for categorical columns from ``project.get_datasource_metadata(...)``."""
    columns = ds_meta.get("columns") or []
    if not isinstance(columns, list):
        return []
    out: list[str] = []
    for col in columns:
        if not isinstance(col, dict):
            continue
        if col.get("datatype") not in CATEGORICAL_DATATYPES:
            continue
        fid = column_field_id(col)
        if fid:
            out.append(fid)
    return out


def numeric_field_ids_from_metadata(ds_meta: dict[str, Any]) -> list[str]:
    """Field ids for numeric columns from datasource metadata."""
    columns = ds_meta.get("columns") or []
    if not isinstance(columns, list):
        return []
    out: list[str] = []
    for col in columns:
        if not isinstance(col, dict):
            continue
        if col.get("datatype") not in NUMERIC_DATATYPES:
            continue
        fid = column_field_id(col)
        if fid:
            out.append(fid)
    return out


def collect_wrapper_subgroup_keys_for_project(project: Any) -> set[str]:
    """All rows-as-columns subgroup keys declared on the observation datasource."""
    roles = infer_datasource_roles(project)
    keys: set[str] = {e.subgroup_key for e in roles.expressions}
    try:
        from mdvtools.llm.column_field_resolve import rows_as_columns_subgroup_keys_from_metadata

        md = project.get_datasource_metadata(roles.obs_datasource)
        keys |= rows_as_columns_subgroup_keys_from_metadata(md)
    except Exception:
        pass
    return keys


def build_chatmdv_roles_constants_block(project: Any) -> str:
    """
    Python source injected into ChatMDV-generated scripts so the LLM does not
    rediscover rows-as-columns metadata at runtime.
    """
    roles = infer_datasource_roles(project)
    expr = roles.preferred_expression()
    lines = [
        "# --- ChatMDV project roles (derived from project metadata; do not edit) ---",
        f"CHATMDV_OBS_DATASOURCE = {json.dumps(roles.obs_datasource)}",
    ]
    if expr is None:
        lines.extend(
            [
                "CHATMDV_EXPR_DATASOURCE = None",
                "CHATMDV_EXPR_NAME_COLUMN = None",
                "CHATMDV_EXPR_SUBGROUP_KEY = None",
            ]
        )
    else:
        lines.extend(
            [
                f"CHATMDV_EXPR_DATASOURCE = {json.dumps(expr.datasource_name)}",
                f"CHATMDV_EXPR_NAME_COLUMN = {json.dumps(expr.name_column)}",
                f"CHATMDV_EXPR_SUBGROUP_KEY = {json.dumps(expr.subgroup_key)}",
            ]
        )
    expr_entries = [
        {
            "datasource": e.datasource_name,
            "name_column": e.name_column,
            "subgroup_key": e.subgroup_key,
        }
        for e in roles.expressions
    ]
    lines.append(f"CHATMDV_EXPRESSIONS = {json.dumps(expr_entries)}")
    try:
        obs_meta = project.get_datasource_metadata(roles.obs_datasource)
        cat_ids = categorical_field_ids_from_metadata(obs_meta)
        if len(cat_ids) > CHATMDV_CATEGORICAL_FIELD_IDS_CAP:
            cat_ids = cat_ids[:CHATMDV_CATEGORICAL_FIELD_IDS_CAP]
        lines.append(f"CHATMDV_CATEGORICAL_FIELD_IDS = {json.dumps(cat_ids)}")
    except Exception:
        lines.append("CHATMDV_CATEGORICAL_FIELD_IDS = []")
    lines.append("# --- end ChatMDV roles ---")
    return "\n".join(lines)


def format_no_hallucination_chart_policy() -> str:
    """
    Prompt text: every chart type must reference real datasources and columns from project metadata only.
    """
    return (
        "- **Metadata-first chart params (all chart types):** Before building *any* chart (DotPlot, Heatmap, Scatter, "
        "BoxPlot, ViolinPlot, TablePlot, SelectionDialogPlot, RowChart, Sankey, etc.), **discover** what exists: "
        "`[d['name'] for d in project.datasources]` and `project.get_datasource_metadata(<ds>)` for each datasource "
        "you reference. **Never** invent a datasource name for `initialCharts` keys or `datasource_name=...`.\n"
        "- **Column field ids:** Every string in `params`, and column slots like `color_by` / `tooltip.column` / "
        "`background_filter.column` when used, must be a **Field ID** from that datasource’s metadata (Project Data "
        "Context), or a valid expression **wrapper** whose subgroup prefix appears under "
        "`links.*.rows_as_columns.subgroups` for the observation datasource. Do not guess `leiden`, `gs`, or `gene_ids` "
        "unless they appear in that metadata.\n"
        "- **Pre-flight check:** Before `project.set_view(...)`, verify columns resolve, e.g. "
        "`project.get_datasource_as_dataframe(obs_ds, columns=[...])` with the same field/wrapper strings as the chart "
        "(optional but recommended when using wrappers). If resolution fails, use bounded `print(...)` and omit broken "
        "charts.\n"
        "- **Column metadata schema:** Each entry in `get_datasource_metadata(ds)['columns']` uses **`datatype`** "
        "(e.g. `text`, `integer`, `double`), **`field`** (chart param id), and **`name`** (display). "
        "**Never** use `col['dtype']` or pandas/AnnData-style keys. Use `CHATMDV_CATEGORICAL_FIELD_IDS` or "
        "`categorical_field_ids_from_metadata(project.get_datasource_metadata(CHATMDV_OBS_DATASOURCE))`.\n"
    )


def format_metadata_column_schema_policy() -> str:
    """Prompt text: MDV column dict schema and multi-gene heatmap param pattern."""
    return (
        "- **Datasource column dicts:** Keys are `field`, `name`, `datatype` (and optionally `values`) — not `dtype`. "
        "Do not build categoricals with `col['dtype'] == 'text'`.\n"
        "- **Categorical discovery:** Prefer injected `CHATMDV_CATEGORICAL_FIELD_IDS` or "
        "`categorical_field_ids_from_metadata(project.get_datasource_metadata(CHATMDV_OBS_DATASOURCE))` "
        "from `mdvtools.llm.datasource_roles`.\n"
        "- **Multi-feature heatmap (e.g. several genes by cell type):** On the observation datasource, "
        "`HeatmapPlot(params=[<categorical_field_id>, <wrapper1>, <wrapper2>, ...])` where the first param groups "
        "rows (cell type / cluster field id from context) and remaining params are expression wrappers built with "
        "`build_expression_wrapper_token` after resolving each feature index from `CHATMDV_EXPR_DATASOURCE`.\n"
    )


def format_feature_table_field_policy(roles: InferredDatasourceRoles) -> str:
    """
    Prompt text for RAG: use actual feature-table field ids and name_column from metadata.

    Avoids assuming a pandas column literally named ``name`` (many projects use ``gene_ids`` etc.).
    """
    if not roles.expressions:
        return (
            "- Feature-table gene labels: inspect `df2.columns` and the Project Data Context **Field ID** "
            "for each feature datasource. Do not assume a column named `name` unless it appears in `df2` "
            "and in the context table."
        )
    lines: list[str] = []
    for e in roles.expressions:
        lines.append(
            f"- For feature datasource `{e.datasource_name}`, the rows-as-columns link declares "
            f"`name_column`=`{e.name_column}`. Use that column for feature labels in code and the **same Field ID** "
            f"in chart ``params`` (e.g. RowSummaryBox on `genes`). Do not substitute `name` unless it matches."
        )
        lines.append(
            f"- **Wrapper `param` tokens on `{roles.obs_datasource}`:** expression columns use "
            f'`{e.subgroup_key}|<gene>({e.subgroup_key})|<index>` — the first segment **must** be `{e.subgroup_key}` '
            "(the subgroup key from metadata). **Never** use example names like `rna_expr` or `protein_expr` unless "
            "that exact key appears in Project Data Context / Datasource roles for this project; wrong keys cause "
            "frontend errors (`subgroup` undefined)."
        )
    return "\n".join(lines)


# Conventional datasource name for ChatMDV-persisted Scanpy marker tables (long format).
CHAT_RANK_GENES_DATASOURCE_NAME = "chat_rank_genes_result"


def format_marker_ranking_viz_policy() -> str:
    """
    Prompt text: do not use DotPlot/Heatmap on ``cells`` as a stand-in for ``rank_genes_groups``; optional scratch DS.

    Wrapper charts on the observation datasource read the MDV matrix—they cannot reproduce the full per-cluster
    statistics table from Scanpy and often show wrong gene or cluster counts.
    """
    ds = CHAT_RANK_GENES_DATASOURCE_NAME
    # One complete call block so the model does not emit a dangling `'chat_rank_genes_result',` line (IndentationError).
    example_block = (
        "\n  ```python\n"
        "  project.add_datasource(\n"
        f'      "{ds}",\n'
        "      marker_df,\n"
        "      replace_data=True,\n"
        "      add_to_view=view_name,\n"
        "  )\n"
        "  ```\n"
    )
    return (
        "- **`rank_genes_groups` vs DotPlot / Heatmap on `cells`:** If the primary answer is a **Scanpy marker table** "
        "(`sc.tl.rank_genes_groups`, `adata.uns['rank_genes_groups']`, or the same printed long-format `DataFrame`), "
        "**do not** add **DotPlot**, **Heatmap**, or other **wrapper-`param` expression charts** on **`cells`** to "
        "represent that table—those charts use the **MDV** expression matrix and a **subset** of genes; they are **not** "
        "the ranked statistics (scores, p-values, all top-N per cluster) and often show **wrong** gene or cluster counts.\n"
        "- **Still allowed:** DotPlot / Heatmap on `cells` for **user-named genes** or exploratory expression when you "
        "are **not** claiming the chart equals the full `rank_genes_groups` output.\n"
        "- **Interpretation / cell-type questions:** Prefer **print + markdown** without `add_datasource` unless the user "
        "clearly needs the marker table **in the saved project view** (see section 7 \"Precedence vs marker persistence\").\n"
        "- **In-view table (optional):** persist the printed long-format `DataFrame` with **one** complete call—either "
        f"the single line `project.add_datasource('{ds}', marker_df, replace_data=True, add_to_view=view_name)` "
        "or the multi-line form below. Column names become field ids; a default table is attached. **Do not** then call "
        "`set_view` with a **new** `initialCharts` dict that **drops** this datasource—merge with "
        "`project.get_view(view_name)` if you add charts, or rely on `add_datasource` alone."
        + example_block
        + "- **Syntax guard:** **Never** output a **fragment** such as a lone line `'"
        + ds
        + "',` or any indented line that is only a string and comma—**that causes `IndentationError`**. Omit "
        "`add_datasource` entirely if you are not using it.\n"
        "- **Print-only view:** If you do not add that datasource, **do not** add a wrapper DotPlot/Heatmap on `cells` for "
        "the marker table; you may call `set_view` with an **empty** `initialCharts` mapping or omit `set_view` after "
        "printing.\n"
    )


def format_obs_table_chart_param_policy() -> str:
    """
    Prompt text: chart ``param`` strings must match Field IDs per datasource; Scanpy marker tables default to chat stdout.

    For marker tables: do not invent ``cells`` field ids; optional in-view copy via ``add_datasource`` to
    ``CHAT_RANK_GENES_DATASOURCE_NAME`` (see body). Avoid ad-hoc ``set_column`` on ``cells`` for long-format marker output.
    """
    ds = CHAT_RANK_GENES_DATASOURCE_NAME
    return (
        "- **Table chart / TablePlot and `params`:** Each string in `params` must be a **Field ID** listed in Project "
        "Data Context for **the same datasource** as the chart (the key under `initialCharts`, e.g. `cells`). "
        "Field IDs are per datasource: a column like `gene_ids` on the **`genes`** datasource is **not** valid as a "
        "`param` on a chart bound to **`cells`** unless `cells` also lists that Field ID.\n"
        "- **Scanpy column names on `cells`:** Names from `rank_genes_groups` or `DataFrame.columns` (e.g. `gene`, `score`) "
        "do **not** automatically exist on `cells`. Do **not** put those into `table_chart` **on `cells`** unless listed "
        "in context.\n"
        f"- **Scratch marker table datasource:** After `add_datasource('{ds}', marker_df, ...)`, charts under "
        f"`initialCharts['{ds}']` use field ids from that `DataFrame`—that is the supported in-view table path.\n"
        f"- **Marker datasource param source of truth:** For `table_chart` on `{ds}`, derive `params` from actual "
        "persisted fields (e.g. `project.get_datasource_metadata(...)['columns']` and/or `marker_df.columns`) rather "
        "than hard-coded names. If your requested marker tokens (`cluster/gene/...` or Scanpy "
        "`group/names/...`) are missing after persistence, do not save a broken chart.\n"
        "- **Post-write marker check:** After writing marker table datasource, verify every chart param token exists in "
        "that datasource field set before `set_view`; if not, fallback to bounded `print(...)` output only.\n"
        "- **Otherwise:** Show the table with **bounded `print(...)`**. If you do not use `add_datasource`, do **not** add "
        "`table_chart` on `cells` with guessed column names.\n"
        "- **Rare:** `project.set_column(...)` on `cells` only for **one value per observation row** aligned with `cells`—"
        "not for long-format marker tables.\n"
    )


def format_marker_gene_scanpy_fallback_policy(path_to_data: str) -> str:
    """
    Prompt text: marker-gene requests must not assume cluster/DE columns exist on the ``genes`` table.

    When a ``.h5ad`` path is available, prefer computing markers in Scanpy and printing to chat.
    """
    p = (path_to_data or "").strip()
    has_h5ad = bool(p) and p.lower().endswith(".h5ad")

    semantic = (
        "- **Marker genes / top genes per cluster:** Cluster labels (e.g. `leiden`, `final_analysis`) live on the "
        "**observation (cells)** datasource or `adata.obs`, not on the **genes** feature table. Per-gene DE stats "
        "(e.g. `dge_mean_diff`) may exist on `genes` *only if* listed in Project Data Context for that datasource.\n"
        "- **Never** assert that `genes` contains `leiden` or other cell-level columns. Do not `raise ValueError` "
        "requiring columns that are not present in `project.get_datasource_as_dataframe('genes').columns` / context.\n"
        "- **Do not** save a `table_chart`/`TablePlot` on `cells` using Scanpy marker `DataFrame` column names unless "
        "those exact Field IDs exist on `cells` in Project Data Context. **Primary** answer: **bounded `print(...)`** in "
        "chat; use **`add_datasource('chat_rank_genes_result', ...)`** only when an **in-view** copy of that table is "
        "appropriate (see sections 3 and 7 and \"Marker ranking vs DotPlot\").\n"
    )

    ds = CHAT_RANK_GENES_DATASOURCE_NAME
    if has_h5ad:
        return (
            semantic
            + "- If the MDV `genes` table lacks columns needed for ranking (or you need full expression), and "
            "`data_path` is a `.h5ad` file: load `adata = sc.read_h5ad(data_path)`, pick the cluster column from "
            "`adata.obs.columns`, then compute markers with Scanpy (e.g. `sc.tl.rank_genes_groups` with "
            "`groupby=<cluster_key>`, or mean expression per group). Print a **bounded** table (top 5 per cluster) "
            "via `print(...)` as the **primary** answer.\n"
            + "- **In-view table (optional):** When the user clearly needs the marker table **in the saved view**, you may "
            "call "
            f"`project.add_datasource('{ds}', marker_df, replace_data=True, add_to_view=view_name)` "
            f"to persist the long-format marker `DataFrame` under the conventional name `{ds}` (see \"Marker ranking vs "
            "DotPlot\" in Parameter Handling). **Do not** add other arbitrary datasources.\n"
            + "- When using AnnData for markers, load `MDVProject(project_path)` for context; `add_datasource` is **only** "
            f"allowed for that `{ds}` marker table when needed—not for unrelated uploads.\n"
        )

    return (
        semantic
        + "- If there is **no** `.h5ad` at `data_path`, compute from MDV only: use **cells** for cluster ids and "
        "expression via row-datasource wrappers or bulk reads; do not require DE columns on `genes` unless they exist. "
        "For marker listings, prefer **bounded `print(...)`** over a `table_chart` with guessed column names.\n"
    )


def format_visualization_consistency_policy() -> str:
    """
    Prompt text: charts must reflect the same pipeline as printed tables (no Scanpy stdout + unrelated MDV wrappers).
    """
    return (
        "- **Single source of truth:** Any chart in `project.set_view(...)` must reflect the **same** quantitative "
        "pipeline as tables or summaries you `print(...)` in this script—not a second independent computation unless "
        "you state they are equivalent.\n"
        "- **Scanpy / AnnData primary:** If results come from `adata` (e.g. `sc.tl.rank_genes_groups`, means from "
        "`adata[...].X`), **do not** default to HeatmapPlot/DotPlot built only from **wrapper** `params` "
        "(`<subgroup>|<gene>(<subgroup>)|<index>`) to show “the same” figure—wrappers read the MDV project matrix, "
        "which may differ from the `.h5ad` used in Scanpy.\n"
        "- **`rank_genes_groups` tables:** Do not use DotPlot/Heatmap on `cells` with wrappers as a substitute for the "
        "full printed marker statistics table; see \"Marker ranking vs DotPlot\" under Parameter Handling.\n"
        "- **Preferred instead:** bounded `print(...)` of the Scanpy/pandas result; optional `TablePlot`/`TextBox` "
        "from the **same** `DataFrame`; or charts using **only** MDV datasources/wrappers end-to-end with no parallel "
        "Scanpy table for the same visualization.\n"
        "- **Avoid by default:** pairing Scanpy-printed marker or expression tables with a separate wrapper-based "
        "expression heatmap that is not explicitly derived from the same array you printed.\n"
    )


def format_scanpy_hybrid_routing_policy() -> str:
    """
    Prompt text: hybrid routing contract for Scanpy outputs (row-aligned -> cells; non-row-aligned -> datasource table).
    """
    ds = CHAT_RANK_GENES_DATASOURCE_NAME
    return (
        "- **Hybrid Scanpy routing contract (core-first):** Route outputs by row alignment, not by convenience.\n"
        "- **Write to `cells` only when 1-row-per-observation:** cluster labels (`leiden`/`louvain`), subcluster labels "
        "merged back by explicit id, and embedding coordinates (`X_umap_*`, `X_pca_*`) that map directly to observations.\n"
        f"- **Use datasource tables for non-row-aligned outputs:** marker/DE long-format tables (e.g. `rank_genes_groups`) "
        f"belong in `{ds}`-style datasources, not ad-hoc columns on `cells`.\n"
        "- **ID-safe merge required for `set_column(...)`:** align by explicit observation id (e.g. `cell_id`) and never "
        "by implicit DataFrame index position.\n"
        "- **Pre-chart validation gate:** before creating charts on a new field, verify field metadata exists on the target "
        "datasource and verify non-missing/category counts are sensible; if unresolved, fall back to bounded `print(...)`.\n"
        f"- **Marker table strictness (`{ds}`):** derive marker table chart params from datasource metadata after "
        "write/replace, not assumptions. Accept common Scanpy aliases only when they resolve to real persisted fields.\n"
        "- **Intent routing examples:**\n"
        "  - marker list / top-N DE -> bounded text table first; optional persisted table datasource when user asks for a "
        "saved view or text box.\n"
        "  - text box requests -> datasource-bound table/text on marker datasource, not guessed marker-stat fields on `cells`.\n"
        "  - heatmap / dot / bubble / violin for marker genes -> expression visualization on valid expression fields/wrappers; "
        "do not claim these replace full DE-stat tables.\n"
        "  - cluster distribution/counts -> aggregate from `cells.<cluster_key>` and optionally show ring/bar chart.\n"
        "  - cell-type prediction from markers -> text-first; write `predicted_cell_type` to `cells` only when mapping is "
        "explicit and one-value-per-cell can be derived.\n"
        "- **Follow-up resolution:** phrases like “those genes”, “same”, or “visualize that” should resolve to the most "
        "recent marker result; if none exists, recompute bounded top-N markers first.\n"
        "- **When user says “use scanpy”:** prefer AnnData/Scanpy compute path when available, while still enforcing this "
        "routing contract for where outputs are persisted.\n"
    )

