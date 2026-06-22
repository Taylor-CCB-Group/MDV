"""
Example: UMAP scatter colored by gene expression, plus SelectionDialogPlot for subset filtering.

Demonstrates the ChatMDV-safe pattern for "expression on embedding within a group":
- Full-dataset ScatterPlot with wrapper-based `set_color_by`
- SelectionDialogPlot with categorical obs field ids (interactive filter in the UI)
- No row-index or invented chart filter methods (`set_row_indices`, `set_background_filter`, etc.)
"""
from __future__ import annotations

import json
import os

from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.selection_dialog_plot import SelectionDialogPlot
from mdvtools.llm.column_field_resolve import build_expression_wrapper_token
from mdvtools.llm.datasource_roles import infer_datasource_roles
from mdvtools.mdvproject import MDVProject


def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))


def main():
    project_path = os.environ.get("MDV_PROJECT", "/app/mdv/pbmc3k_chat")
    view_name = "UMAP gene expression with selection filter (example)"
    obs_ds = "cells"
    gene = os.environ.get("GENE", "RALY")
    # Categorical field id for interactive subset filtering (use Project Data Context field id).
    filter_field = os.environ.get("FILTER_FIELD", "leiden")

    project = MDVProject(project_path, delete_existing=False)
    roles = infer_datasource_roles(project)
    expr = roles.preferred_expression()
    if expr is None:
        raise RuntimeError("No rows-as-columns expression link; cannot build gene wrappers.")

    df_var = project.get_datasource_as_dataframe(
        expr.datasource_name, columns=[expr.name_column]
    )
    names = df_var[expr.name_column].astype(str).tolist()
    idx = names.index(gene)
    wrapper = build_expression_wrapper_token(expr.subgroup_key, gene, idx)

    scatter = ScatterPlot(
        title=f"{gene} expression on UMAP",
        params=["X_umap_1", "X_umap_2"],
        size=[792, 472],
        position=[10, 10],
    )
    scatter.set_color_by(wrapper)
    scatter.set_default_color(wrapper)
    scatter.set_brush("poly")
    scatter.set_opacity(0.85)
    scatter.set_radius(4)

    selection = SelectionDialogPlot(
        title="Filter cells",
        params=[filter_field],
        size=[320, 472],
        position=[820, 10],
    )

    project.set_view(
        view_name,
        {
            "initialCharts": {
                obs_ds: [
                    convert_plot_to_json(scatter),
                    convert_plot_to_json(selection),
                ]
            }
        },
    )
    project.set_editable(True)


if __name__ == "__main__":
    main()
