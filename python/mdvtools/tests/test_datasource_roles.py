from mdvtools.llm.datasource_roles import (
    format_feature_table_field_policy,
    format_marker_gene_scanpy_fallback_policy,
    format_visualization_consistency_policy,
    infer_datasource_roles,
)


class FakeProject:
    def __init__(self):
        self._names = ["cells", "rna", "protein"]

    def get_datasource_names(self):
        return self._names

    def get_links(self, datasource, filter=None):
        assert datasource == "cells"
        assert filter == "rows_as_columns"
        return [
            {
                "datasource": "rna",
                "link": {
                    "rows_as_columns": {
                        "name_column": "name",
                        "name": "rna",
                        "subgroups": {
                            "rna_expr": {"name": "rna_expr", "label": "rna_expr", "type": "sparse"},
                            "rna_logged_counts": {"name": "rna_logged_counts", "label": "rna_logged_counts", "type": "sparse"},
                        },
                    }
                },
            },
            {
                "datasource": "protein",
                "link": {
                    "rows_as_columns": {
                        "name_column": "name",
                        "name": "protein",
                        "subgroups": {
                            "protein_expr": {"name": "protein_expr", "label": "protein_expr", "type": "sparse"},
                            "protein_clr": {"name": "protein_clr", "label": "protein_clr", "type": "dense"},
                        },
                    }
                },
            },
        ]


def test_infer_datasource_roles_prefers_cells_and_expr_subgroups():
    roles = infer_datasource_roles(FakeProject())
    assert roles.obs_datasource == "cells"
    assert [e.datasource_name for e in roles.expressions] == ["rna", "protein"]
    assert roles.expressions[0].name_column == "name"
    assert roles.expressions[0].subgroup_key == "rna_expr"
    assert roles.expressions[1].subgroup_key == "protein_expr"


def test_format_feature_table_field_policy_lists_name_columns():
    roles = infer_datasource_roles(FakeProject())
    text = format_feature_table_field_policy(roles)
    assert "name_column" in text
    assert "`rna`" in text or "rna" in text
    assert "Do not substitute" in text or "gene_ids" in text or "name_column" in text


def test_format_visualization_consistency_policy():
    text = format_visualization_consistency_policy()
    assert "Single source of truth" in text
    assert "Scanpy" in text or "AnnData" in text
    assert "wrapper" in text.lower()


def test_format_marker_gene_scanpy_fallback_with_h5ad():
    text = format_marker_gene_scanpy_fallback_policy("/data/adata.h5ad")
    assert "rank_genes_groups" in text
    assert "Never" in text or "never" in text.lower()
    assert "leiden" in text.lower()


def test_format_marker_gene_scanpy_fallback_without_h5ad():
    text = format_marker_gene_scanpy_fallback_policy("")
    assert "Marker genes" in text or "marker" in text.lower()
    assert "MDV only" in text or "wrappers" in text.lower()


def test_format_feature_table_field_policy_without_expressions():
    class NoLinkProject:
        def get_datasource_names(self):
            return ["cells"]

        def get_links(self, *_a, **_k):
            return []

    roles = infer_datasource_roles(NoLinkProject())
    text = format_feature_table_field_policy(roles)
    assert "df2.columns" in text
    assert "Field ID" in text

