"""
Project scale and memory budget helpers for ChatMDV prompt routing.

Used to inject MDV-first / Scanpy-last-resort guidance into agent and RAG prompts.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any

from mdvtools.llm.datasource_roles import infer_datasource_roles

_BYTES_PER_CELL_ESTIMATE = 16


def _large_project_row_threshold() -> int:
    raw = os.environ.get("CHATMDV_LARGE_PROJECT_ROWS", "100000").strip()
    try:
        return max(1, int(raw))
    except ValueError:
        return 100_000


def _available_ram_mb() -> float:
    try:
        import psutil

        return psutil.virtual_memory().available / (1024 * 1024)
    except Exception:
        return 0.0


@dataclass(frozen=True)
class ProjectScale:
    obs_rows: int
    obs_columns: int
    estimated_obs_df_mb: float
    available_ram_mb: float
    is_large: bool
    has_h5ad: bool
    obs_datasource: str


def assess_project_scale(project: Any, path_to_data: str) -> ProjectScale:
    """Estimate observation-table scale and whether the project is 'large' for ChatMDV."""
    roles = infer_datasource_roles(project)
    obs_ds = roles.obs_datasource
    obs_rows = 0
    obs_columns = 0
    try:
        meta = project.get_datasource_metadata(obs_ds)
        obs_rows = int(meta.get("size") or 0)
        cols = meta.get("columns") or []
        obs_columns = len(cols)
    except Exception:
        pass

    estimated_obs_df_mb = (obs_rows * obs_columns * _BYTES_PER_CELL_ESTIMATE) / (1024 * 1024)
    available_ram_mb = _available_ram_mb()
    threshold = _large_project_row_threshold()
    p = (path_to_data or "").strip()
    has_h5ad = bool(p) and p.lower().endswith(".h5ad")

    return ProjectScale(
        obs_rows=obs_rows,
        obs_columns=obs_columns,
        estimated_obs_df_mb=round(estimated_obs_df_mb, 1),
        available_ram_mb=round(available_ram_mb, 1),
        is_large=obs_rows >= threshold,
        has_h5ad=has_h5ad,
        obs_datasource=obs_ds,
    )


_AGENT_OBS_COLUMN_CAP = 20
_NUMERIC_DTYPES = frozenset({"integer", "double", "number", "float"})
_CATEGORICAL_DTYPES = frozenset({"text", "string", "multitext", "category"})
_EMBEDDING_HINTS = ("umap", "pca", "tsne", "embed")


def agent_probe_columns_from_metadata(
    project: Any,
    datasource_name: str,
    *,
    cap: int = _AGENT_OBS_COLUMN_CAP,
) -> list[str]:
    """Field ids to load for the pandas agent on large projects (schema probing only)."""
    try:
        meta = project.get_datasource_metadata(datasource_name)
    except Exception:
        return []
    columns = meta.get("columns") or []
    picked: list[str] = []
    for col in columns:
        field = col.get("field")
        if not field:
            continue
        dtype = str(col.get("datatype") or "").lower()
        field_lower = field.lower()
        if dtype in _CATEGORICAL_DTYPES or dtype in _NUMERIC_DTYPES:
            picked.append(field)
        elif any(hint in field_lower for hint in _EMBEDDING_HINTS):
            picked.append(field)
        if len(picked) >= cap:
            break
    if not picked and columns:
        first = columns[0].get("field")
        if first:
            picked.append(first)
    return picked


def agent_expression_probe_columns(
    project: Any,
    expression_datasource: str,
    name_column: str,
    *,
    cap: int = 12,
) -> list[str]:
    """Feature-table columns for agent probing (name column + stat columns, not gene matrix)."""
    try:
        meta = project.get_datasource_metadata(expression_datasource)
    except Exception:
        return [name_column] if name_column else []
    columns = meta.get("columns") or []
    picked: list[str] = []
    if name_column:
        picked.append(name_column)
    for col in columns:
        field = col.get("field")
        if not field or field in picked:
            continue
        dtype = str(col.get("datatype") or "").lower()
        if dtype in _NUMERIC_DTYPES:
            picked.append(field)
        if len(picked) >= cap:
            break
    return picked or ([name_column] if name_column else [])


def load_agent_dataframes(
    project: Any,
    roles: Any,
    scale: ProjectScale,
    *,
    extra_datasource: str | None = None,
) -> list[Any]:
    """
    Return [obs_df] or [obs_df, expr_df] for the pandas agent.

    Large projects load column subsets only; small projects load full tables.
    When there is no expression datasource and the project has multiple tabular tables,
    ``extra_datasource`` (typically question-resolved) is loaded as df2 for agent probing.
    """
    obs_ds = roles.obs_datasource
    if scale.is_large:
        obs_cols = agent_probe_columns_from_metadata(project, obs_ds)
        df1 = (
            project.get_datasource_as_dataframe(obs_ds, columns=obs_cols)
            if obs_cols
            else project.get_datasource_as_dataframe(obs_ds)
        )
    else:
        df1 = project.get_datasource_as_dataframe(obs_ds)

    primary_expr = roles.preferred_expression()
    if primary_expr is None:
        names: list[str] = []
        try:
            names = list(project.get_datasource_names())
        except Exception:
            pass
        if (
            extra_datasource
            and extra_datasource != obs_ds
            and len(names) > 2
            and "cells" not in names
        ):
            try:
                df_extra = project.get_datasource_as_dataframe(extra_datasource)
                return [df1, df_extra]
            except Exception:
                pass
        return [df1]

    expr_ds = primary_expr.datasource_name
    if scale.is_large:
        expr_cols = agent_expression_probe_columns(
            project, expr_ds, primary_expr.name_column
        )
        df2 = (
            project.get_datasource_as_dataframe(expr_ds, columns=expr_cols)
            if expr_cols
            else project.get_datasource_as_dataframe(expr_ds)
        )
    else:
        df2 = project.get_datasource_as_dataframe(expr_ds)
    return [df1, df2]


def format_scale_context_block(scale: ProjectScale) -> str:
    """Markdown snippet for agent/RAG prompts."""
    ram_note = (
        f"available RAM ~{scale.available_ram_mb:.0f} MB"
        if scale.available_ram_mb > 0
        else "available RAM unknown"
    )
    size_note = (
        f"~{scale.estimated_obs_df_mb:.0f} MB if all {scale.obs_columns} metadata columns are loaded"
        if scale.obs_columns > 0
        else "metadata column count unknown"
    )
    large_note = (
        "Treat as a **large project**: prefer MDV chart APIs and column-subset loads; "
        "avoid full `get_datasource_as_dataframe` and full in-memory AnnData."
        if scale.is_large
        else "Small/medium project: column-subset loads are still preferred when only a few fields are needed."
    )
    return (
        f"- Observation datasource `{scale.obs_datasource}`: **{scale.obs_rows:,}** rows, "
        f"**{scale.obs_columns}** metadata columns ({size_note}; {ram_note}).\n"
        f"- {large_note}"
    )
