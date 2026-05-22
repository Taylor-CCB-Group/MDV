"""
Example: Leiden cluster marker heatmap using Scanpy DE genes and metadata-valid wrappers.

Before `set_view`, discover the rows-as-columns subgroup key from the project (never hardcode `gs` / `rna_expr`).
Build wrappers with `mdvtools.llm.column_field_resolve.build_expression_wrapper_token`.

Requires an existing MDV project with `cells` + expression link (e.g. from `convert_scanpy_to_mdv`).

HeatmapPlot uses set_x_axis / set_y_axis / set_axis below — not set_axis_properties. For BoxPlot, ViolinPlot,
ScatterPlot, DotPlot use set_axis_properties("x", {...}) only.
"""
from __future__ import annotations

import json
import os

import scanpy as sc

from mdvtools.charts.heatmap_plot import HeatmapPlot
from mdvtools.charts.selection_dialog_plot import SelectionDialogPlot
from mdvtools.llm.column_field_resolve import build_expression_wrapper_token
from mdvtools.llm.datasource_roles import infer_datasource_roles
from mdvtools.mdvproject import MDVProject


def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))


def _pick_cluster_key(adata: sc.AnnData) -> str:
    for k in ("leiden", "louvain"):
        if k in adata.obs.columns:
            return k
    raise ValueError("No leiden/louvain in adata.obs; pick a categorical cluster column.")


def main():
    project_path = os.environ.get("MDV_PROJECT", os.path.expanduser("~/mdv/project"))
    data_path = os.environ.get("H5AD_PATH", "path_to_data.h5ad")
    view_name = "Heatmap of rank_genes markers per Leiden"
    obs_ds = "cells"

    project = MDVProject(project_path, delete_existing=False)
    roles = infer_datasource_roles(project)
    expr = roles.preferred_expression()
    if expr is None:
        raise RuntimeError("No rows-as-columns expression link; cannot build gene wrappers.")
    subgroup_key = expr.subgroup_key
    name_column = expr.name_column

    adata = sc.read_h5ad(data_path)
    cluster_key = _pick_cluster_key(adata)
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")

    # Top N markers per cluster (union), capped for a readable heatmap
    n_per = 3
    cap = 24
    deg = sc.get.rank_genes_groups_df(adata, group=None)
    seen: set[str] = set()
    ordered: list[str] = []
    for grp in deg["group"].unique():
        sub = deg[deg["group"] == grp].head(n_per)
        for g in sub["names"].astype(str):
            if g not in seen:
                seen.add(g)
                ordered.append(g)
            if len(ordered) >= cap:
                break
        if len(ordered) >= cap:
            break

    var_names = [str(x) for x in adata.var_names]
    label_to_idx: dict[str, int] = {}
    if name_column in adata.var.columns:
        labels = adata.var[name_column].astype(str)
        for i, lab in enumerate(labels):
            label_to_idx.setdefault(str(lab), i)
    for g in ordered:
        if g not in label_to_idx:
            try:
                label_to_idx[g] = var_names.index(g)
            except ValueError:
                label_to_idx[g] = -1
    ordered = [g for g in ordered if label_to_idx.get(g, -1) >= 0]

    gene_params = [
        build_expression_wrapper_token(subgroup_key, g, label_to_idx[g]) for g in ordered
    ]
    heatmap_params = [cluster_key] + gene_params

    heatmap_plot = HeatmapPlot(
        title="Marker expression (rank_genes_groups) by cluster",
        params=heatmap_params,
        size=[1200, 472],
        position=[10, 10],
    )
    heatmap_plot.set_color_scale({"log": False})
    heatmap_plot.set_x_axis(axis_labels="Cluster", axis_title=cluster_key)
    heatmap_plot.set_y_axis(axis_labels="Genes", axis_title="Markers")
    heatmap_plot.set_axis(
        xtextSize=10,
        ytextSize=10,
        xsize=110,
        ysize=110,
        xtickfont=8,
        ytickfont=8,
        xrotate_status=True,
        yrotate_status=True,
    )

    sel = SelectionDialogPlot(
        title="Filter clusters / genes",
        params=heatmap_params,
        size=[400, 300],
        position=[10, 500],
    )

    view_config = {
        "initialCharts": {
            obs_ds: [
                convert_plot_to_json(heatmap_plot),
                convert_plot_to_json(sel),
            ]
        }
    }
    project.set_view(view_name, view_config)
    project.set_editable(True)


if __name__ == "__main__":
    main()
