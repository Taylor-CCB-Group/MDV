"""Name clashes, the gates column and UMAP roles in project curation."""

import json
from pathlib import Path

from mdvtools.examples.curate_mdv_project import column_role, gene_panel
from mdvtools.examples.rename_softdelete_mdv_project_fields import (
    apply_plan,
    build_soft_delete_plan,
    planned_soft_deletes,
)
from mdvtools.template_views.create_question_views import compound_dot_plot, genes_panel


def write_project(tmp_path: Path, columns: list[dict]) -> Path:
    (tmp_path / "datasources.json").write_text(
        json.dumps([{"name": "cells", "columns": columns}]), encoding="utf-8"
    )
    return tmp_path


def rename(field: str, target: str) -> dict:
    return {
        "row": 2,
        "field": field,
        "display_name": field,
        "current_name": field,
        "target_name": target,
        "action": "rename",
    }


COLUMNS = [
    {"field": "spatial_region", "name": "spatial_region"},
    {"field": "region", "name": "region"},
]


def test_rename_clashing_with_a_visible_column_is_blocked(tmp_path):
    project = write_project(tmp_path, COLUMNS)
    _summary, results, _log = apply_plan(project, "cells", [rename("spatial_region", "Region")], True)
    assert results[0]["status"] == "blocked"


def test_rename_ignores_a_column_about_to_be_soft_deleted(tmp_path):
    project = write_project(tmp_path, COLUMNS)
    _summary, results, _log = apply_plan(
        project, "cells", [rename("spatial_region", "Region")], True, {"region"}
    )
    assert results[0]["status"] == "planned"


def test_rename_ignores_an_already_soft_deleted_column(tmp_path):
    columns = [COLUMNS[0], {**COLUMNS[1], "deleted": True}]
    project = write_project(tmp_path, columns)
    _summary, results, _log = apply_plan(project, "cells", [rename("spatial_region", "Region")], True)
    assert results[0]["status"] == "planned"


def test_gates_column_is_never_soft_deleted(tmp_path):
    project = write_project(
        tmp_path,
        [{"field": "__gates__", "name": "gates"}, {"field": "nCount_RNA", "name": "nCount_RNA"}],
    )
    plan = build_soft_delete_plan(project, "cells", set(), set(), {})
    assert planned_soft_deletes(plan) == {"nCount_RNA"}


def test_a_column_used_in_a_view_is_not_counted_as_hidden(tmp_path):
    project = write_project(tmp_path, COLUMNS)
    views = {"View": {"initialCharts": {"cells": [{"param": ["region"]}]}}}
    plan = build_soft_delete_plan(project, "cells", {"spatial_region"}, set(), views)
    assert planned_soft_deletes(plan) == set()


def test_umap_dimension_names_are_umap_roles():
    assert column_role({"name": "UMAP Dimension 1"}) == "umap_x"
    assert column_role({"name": "UMAP Dimension 2"}) == "umap_y"


def test_dot_plot_without_a_split_has_one_panel():
    chart = compound_dot_plot("Genes by Cell Type", "cell_type", ["gene_a"])
    assert chart["type"] == "compound_dot_plot"
    assert chart["split_by"] is None
    assert chart["split_values"] == []


def test_genes_pane_uses_the_gene_info_chart():
    charts = gene_panel("genes", ["CD14"])
    assert [chart["type"] for chart in charts] == ["selection_dialog", "geneNetwork"]
    assert charts[1]["param"] == ["name"]


def test_myeloid_genes_pane_uses_the_gene_info_chart():
    assert [chart["type"] for chart in genes_panel()] == ["selection_dialog", "geneNetwork"]
