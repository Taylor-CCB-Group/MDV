"""
Example: UMAP scatter colored by gene expression via metadata-valid wrapper.

Uses injected-style constants pattern: discover roles from the project, never hardcode gs/genes.
"""
from __future__ import annotations

import json
import os

from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.llm.column_field_resolve import build_expression_wrapper_token
from mdvtools.llm.datasource_roles import infer_datasource_roles
from mdvtools.mdvproject import MDVProject


def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))


def main():
    project_path = os.environ.get("MDV_PROJECT", "/app/mdv/pbmc3k_chat")
    view_name = "UMAP gene expression (example)"
    obs_ds = "cells"
    gene = os.environ.get("GENE", "RALY")

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

    project.set_view(
        view_name,
        {"initialCharts": {obs_ds: [convert_plot_to_json(scatter)]}},
    )
    project.set_editable(True)


if __name__ == "__main__":
    main()
