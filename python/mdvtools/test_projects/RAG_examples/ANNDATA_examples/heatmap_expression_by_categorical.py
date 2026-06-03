"""
Example: Heatmap of multiple gene expression columns grouped by a categorical obs field.

Uses MDV metadata schema (``datatype``, ``field``) — not ``col['dtype']``.
Pattern: HeatmapPlot(params=[<categorical_field_id>, <wrapper1>, <wrapper2>, ...]) on obs datasource.
"""
from __future__ import annotations

import json
import os

from mdvtools.charts.heatmap_plot import HeatmapPlot
from mdvtools.llm.column_field_resolve import build_expression_wrapper_token
from mdvtools.llm.datasource_roles import (
    categorical_field_ids_from_metadata,
    infer_datasource_roles,
)
from mdvtools.mdvproject import MDVProject


def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))


def _resolve_gene_wrapper(project, expr, gene: str) -> str:
    df_var = project.get_datasource_as_dataframe(
        expr.datasource_name, columns=[expr.name_column]
    )
    names = df_var[expr.name_column].astype(str).tolist()
    idx = names.index(gene)
    return build_expression_wrapper_token(expr.subgroup_key, gene, idx)


def main():
    project_path = os.environ.get("MDV_PROJECT", "/app/mdv/pbmc3k_chat")
    view_name = "Heatmap: gene expression by cell group (example)"
    obs_ds = "cells"
    genes = os.environ.get("GENES", "NKG7,GNLY").split(",")
    group_field = os.environ.get("GROUP_FIELD", "")

    project = MDVProject(project_path, delete_existing=False)
    roles = infer_datasource_roles(project)
    expr = roles.preferred_expression()
    if expr is None:
        raise RuntimeError("No rows-as-columns expression link; cannot build gene wrappers.")

    obs_meta = project.get_datasource_metadata(roles.obs_datasource)
    cat_fields = categorical_field_ids_from_metadata(obs_meta)
    if not group_field:
        group_field = cat_fields[0] if cat_fields else ""
    if not group_field:
        raise RuntimeError("No categorical field id available for heatmap grouping axis.")

    wrappers = [_resolve_gene_wrapper(project, expr, g.strip()) for g in genes if g.strip()]
    if len(wrappers) < 1:
        raise RuntimeError("At least one gene wrapper is required for heatmap params.")

    heatmap = HeatmapPlot(
        title="Expression by cell group",
        params=[group_field, *wrappers],
        size=[900, 500],
        position=[10, 10],
    )
    heatmap.set_color_scale({"log": False})

    project.set_view(
        view_name,
        {"initialCharts": {obs_ds: [convert_plot_to_json(heatmap)]}},
    )
    project.set_editable(True)


if __name__ == "__main__":
    main()
