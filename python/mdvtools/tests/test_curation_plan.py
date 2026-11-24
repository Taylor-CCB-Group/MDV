"""Name clashes, the gates column and UMAP roles in project curation."""

import json
from pathlib import Path

from mdvtools.examples.curate_mdv_project import (
    CONCORD_VIEW_NAME,
    PROTEIN_VIEW_NAME,
    Q7_VIEW_NAME,
    arrange_views,
    choose_roles,
    column_role,
    gene_panel,
    opening_names,
    plan_question_views,
    split_category_values,
    stored_name,
    table_kind,
)
from mdvtools.examples.rename_softdelete_mdv_project_fields import (
    apply_plan,
    build_plan,
    build_soft_delete_plan,
    planned_soft_deletes,
)
from mdvtools.template_views.create_question_views import (
    Q1_VIEW_NAME,
    Q3_VIEW_NAME,
    compound_dot_plot,
    genes_panel,
)


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


def test_workbook_renames_bucket_and_annotation(tmp_path):
    project = write_project(
        tmp_path,
        [
            {"field": "bucket", "name": "bucket"},
            {"field": "annotation", "name": "annotation"},
            {"field": "sub_bucket", "name": "sub_bucket"},
        ],
    )
    workbook = tmp_path / "cleaning.xlsx"
    import pandas as pd

    pd.DataFrame(
        [
            {
                "datasource": "cells",
                "original name": "bucket",
                "displayed name": "bucket",
                "new_display_name": "Cell Group",
                "Final Column Name": "Cell Group",
                "Project id": "mdv_gut",
            },
            {
                "datasource": "cells",
                "original name": "annotation",
                "displayed name": "annotation",
                "new_display_name": "Cell Type",
                "Final Column Name": "Cell Type",
                "Project id": "mdv_gut",
            },
        ]
    ).to_excel(workbook, sheet_name="Data Cleaning", index=False)
    plan, _warnings, _errors, kept, _ambiguous = build_plan(project, "cells", workbook, "mdv_gut")
    by_field = {item["field"]: item for item in plan}
    assert by_field["bucket"]["action"] == "rename"
    assert by_field["bucket"]["target_name"] == "Cell Group"
    assert by_field["annotation"]["target_name"] == "Cell Type"
    assert by_field["sub_bucket"]["action"] == "keep"
    assert "sub_bucket" in kept


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


def test_gene_prompt_sits_above_the_list():
    charts = gene_panel("genes", ["PTPRC"], "Select gene(s) to view")
    assert [chart["type"] for chart in charts] == ["text_box_chart", "selection_dialog", "geneNetwork"]
    assert charts[0]["title"] == "Select gene(s) to view"
    assert charts[0]["text"] == "Select gene(s) to view"
    assert charts[1]["filters"]["name"]["category"] == ["PTPRC"]


def test_comparison_lists_start_on_five_genes():
    names = ["TNF", "TNFRSF1A", "TNFRSF1B", "IFNG", "CXCL9", "CXCL10"]
    assert opening_names(names, 5) == names[:5]
    assert opening_names(names, 1) == ["TNF"]
    assert opening_names(["CD14"], 5) == ["CD14"]


def test_split_uses_two_to_five_category_values():
    assert split_category_values({"values": ["UC", "CD", "Healthy"]}) == ["UC", "CD", "Healthy"]
    assert split_category_values({"values": ["skin"]}) is None
    assert split_category_values({"values": [f"type-{i}" for i in range(21)]}) is None


def test_a_short_category_becomes_a_split_dot_plot():
    values = split_category_values({"values": ["UC", "CD"]})
    chart = compound_dot_plot(
        "Cell Group, split by Diagnosis",
        "bucket",
        {"type": "RowsAsColsQuery"},
        split=("Disease", values or []),
    )
    assert chart["split_by"] == "Disease"
    assert chart["split_values"] == ["UC", "CD"]
    assert "pos" not in chart["color_legend"]


def test_major_cell_type_is_preferred_over_a_long_cell_state():
    roles = choose_roles(
        [
            {"name": "Cell State", "field": "final", "values": [str(index) for index in range(40)]},
            {"name": "Major Cell Type", "field": "major", "values": ["A", "B", "C"]},
        ]
    )
    assert roles["cell_type"]["field"] == "major"


def test_derm_and_taurus_column_names_get_roles():
    assert column_role({"name": "Diagnosed Disease"}) == "diagnosis"
    assert column_role({"name": "Condition"}) == "condition"
    assert column_role({"name": "Treatment"}) == "treatment"
    assert column_role({"name": "Presence of Inflammation"}) == "inflammation"
    assert column_role({"name": "Biopsy Site"}) == "tissue"
    assert column_role({"name": "Remission Status"}) == "response"
    assert column_role({"name": "Major Cell Type"}) == "cell_type"


def test_unrenamed_cell_fields_and_prot_table_are_recognized():
    assert column_role({"name": "bucket"}) == "cell_group"
    assert column_role({"name": "sub_bucket"}) == "cell_subgroup"
    assert column_role({"name": "annotation"}) == "cell_type"
    assert table_kind("prot") == 2
    assert table_kind("protein") == 2
    assert table_kind("rna") == 0


def test_protein_names_match_hyphen_or_underscore():
    available = {"adt-CD14", "adt-CD1c", "adt-PD-L1"}
    assert stored_name(available, "adt_CD14") == "adt-CD14"
    assert stored_name(available, "adt_CD1c") == "adt-CD1c"
    assert stored_name(available, "adt_PD-L1") == "adt-PD-L1"
    assert stored_name(available, "adt_CD99") is None


def test_new_question_views_keep_protein_beside_the_gene_view():
    ordered = arrange_views(
        ["Cell atlas", "summary"],
        {
            Q1_VIEW_NAME: {},
            PROTEIN_VIEW_NAME: {},
            "Is the expression different between disease stages/patient groups/tissue types?": {},
            CONCORD_VIEW_NAME: {},
        },
    )
    assert ordered[:2] == ["Cell atlas", "summary"]
    assert ordered[ordered.index(Q1_VIEW_NAME) + 1] == PROTEIN_VIEW_NAME
    assert ordered[ordered.index(PROTEIN_VIEW_NAME) + 1] == CONCORD_VIEW_NAME


def test_protein_view_follows_the_gene_view_and_omitted_views_are_dropped():
    ordered = arrange_views(
        [
            "Cell atlas",
            Q1_VIEW_NAME,
            "Is the expression different between disease stages/patient groups/tissue types?",
            CONCORD_VIEW_NAME,
            Q7_VIEW_NAME,
            "MNP cell signatures",
        ],
        {PROTEIN_VIEW_NAME: {}, CONCORD_VIEW_NAME: {}, Q1_VIEW_NAME: {}},
    )
    assert Q7_VIEW_NAME not in ordered
    assert Q3_VIEW_NAME not in ordered
    assert ordered[ordered.index(Q1_VIEW_NAME) + 1] == PROTEIN_VIEW_NAME
    assert CONCORD_VIEW_NAME in ordered
    assert ordered.index(CONCORD_VIEW_NAME) < ordered.index("MNP cell signatures")


def test_signature_view_is_not_created(tmp_path):
    plan = plan_question_views(tmp_path, [], {})
    assert Q7_VIEW_NAME not in plan["views"]
    assert any(item["question"] == Q7_VIEW_NAME for item in plan["skipped"])


def test_myeloid_genes_pane_uses_the_gene_info_chart():
    assert [chart["type"] for chart in genes_panel()] == ["selection_dialog", "geneNetwork"]
