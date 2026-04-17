from __future__ import annotations

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
    return "\n".join(lines)


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
    )

    if has_h5ad:
        return (
            semantic
            + "- If the MDV `genes` table lacks columns needed for ranking (or you need full expression), and "
            "`data_path` is a `.h5ad` file: load `adata = sc.read_h5ad(data_path)`, pick the cluster column from "
            "`adata.obs.columns`, then compute markers with Scanpy (e.g. `sc.tl.rank_genes_groups` with "
            "`groupby=<cluster_key>`, or mean expression per group). Print a **bounded** table (top 5 per cluster) "
            "via `print(...)` as the primary answer. Prefer this path over failing on missing MDV columns.\n"
            + "- When using AnnData for markers, you may still load `MDVProject(project_path)` read-only for context; "
            "do not call `project.add_datasource` unless explicitly creating a new project.\n"
        )

    return (
        semantic
        + "- If there is **no** `.h5ad` at `data_path`, compute from MDV only: use **cells** for cluster ids and "
        "expression via row-datasource wrappers or bulk reads; do not require DE columns on `genes` unless they exist.\n"
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
        "- **Preferred instead:** bounded `print(...)` of the Scanpy/pandas result; optional `TablePlot`/`TextBox` "
        "from the **same** `DataFrame`; or charts using **only** MDV datasources/wrappers end-to-end with no parallel "
        "Scanpy table for the same visualization.\n"
        "- **Avoid by default:** pairing Scanpy-printed marker or expression tables with a separate wrapper-based "
        "expression heatmap that is not explicitly derived from the same array you printed.\n"
    )

