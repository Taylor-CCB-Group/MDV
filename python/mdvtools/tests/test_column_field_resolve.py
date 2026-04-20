from mdvtools.llm.column_field_resolve import (
    normalize_view_chart_params,
    prune_view_charts_with_invalid_params,
    rows_as_columns_subgroup_keys_from_metadata,
    build_expression_wrapper_token,
)


class FakeProjectCells:
    """cells has leiden; rows-as-columns link with subgroup ``rna_expr`` for wrapper tests."""

    datasources = [{"name": "cells"}]

    def get_datasource_metadata(self, name: str):
        assert name == "cells"
        return {
            "name": "cells",
            "columns": [
                {"name": "leiden", "field": "leiden", "datatype": "text"},
            ],
            "links": {
                "genes": {
                    "rows_as_columns": {
                        "name_column": "name",
                        "subgroups": {
                            "rna_expr": {
                                "name": "gene_scores",
                                "label": "Gene Expr",
                                "type": "dense",
                            }
                        },
                    }
                }
            },
        }


class FakeProjectMarkerAliases:
    def get_datasource_metadata(self, name: str):
        assert name == "chat_rank_genes_result"
        return {
            "name": "chat_rank_genes_result",
            "columns": [
                {"name": "cluster", "field": "cluster", "datatype": "text"},
                {"name": "gene", "field": "gene", "datatype": "text"},
                {"name": "score", "field": "score", "datatype": "double"},
                {"name": "logfoldchange", "field": "logfoldchange", "datatype": "double"},
                {"name": "pval", "field": "pval", "datatype": "double"},
                {"name": "pval_adj", "field": "pval_adj", "datatype": "double"},
            ],
        }


class FakeProjectScanpyMarkerColumnNames:
    """Metadata fields match Scanpy ``rank_genes_groups`` DataFrame columns (e.g. ``group``, ``names``)."""

    datasources = [{"name": "chat_rank_genes_result"}]

    def get_datasource_metadata(self, name: str):
        assert name == "chat_rank_genes_result"
        return {
            "name": "chat_rank_genes_result",
            "columns": [
                {"name": "group", "field": "group", "datatype": "integer"},
                {"name": "names", "field": "names", "datatype": "text"},
                {"name": "logfoldchanges", "field": "logfoldchanges", "datatype": "double"},
                {"name": "pvals_adj", "field": "pvals_adj", "datatype": "double"},
                {"name": "scores", "field": "scores", "datatype": "double"},
            ],
        }


def test_prune_drops_table_chart_with_unknown_param_tokens():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "markers",
                    "type": "table_chart",
                    "param": ["leiden", "gene_ids", "score"],
                }
            ]
        }
    }
    proj = FakeProjectCells()
    pruned, n = prune_view_charts_with_invalid_params(view, proj)
    assert n == 1
    assert pruned["initialCharts"]["cells"] == []


def test_prune_keeps_wrapper_params():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "x",
                    "type": "dot_plot",
                    "param": ["rna_expr|CD3D(rna_expr)|12", "leiden"],
                }
            ]
        }
    }
    proj = FakeProjectCells()
    pruned, n = prune_view_charts_with_invalid_params(view, proj)
    assert n == 0
    assert len(pruned["initialCharts"]["cells"]) == 1


def test_prune_drops_hallucinated_subgroup_wrapper():
    """Wrappers must use a subgroup key present in datasource metadata."""
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "bad",
                    "type": "heat_map",
                    "param": ["leiden", "made_up_subgroup|CD3D(made_up_subgroup)|0"],
                }
            ]
        }
    }
    proj = FakeProjectCells()
    pruned, n = prune_view_charts_with_invalid_params(view, proj)
    assert n == 1
    assert pruned["initialCharts"]["cells"] == []


def test_prune_drops_unknown_datasource_key():
    view = {
        "initialCharts": {
            "nonexistent_ds": [{"title": "x", "type": "table_chart", "param": ["a"]}],
        }
    }
    proj = FakeProjectCells()
    pruned, n = prune_view_charts_with_invalid_params(view, proj)
    assert n == 1
    assert "nonexistent_ds" not in pruned.get("initialCharts", {})


def test_subgroup_keys_from_metadata():
    md = FakeProjectCells().get_datasource_metadata("cells")
    assert rows_as_columns_subgroup_keys_from_metadata(md) == {"rna_expr"}


def test_build_expression_wrapper_token():
    assert build_expression_wrapper_token("gs", "MS4A1", 3) == "gs|MS4A1(gs)|3"


def test_normalize_then_prune_removes_bad_marker_chart():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "markers",
                    "type": "table_chart",
                    "param": ["leiden", "gene_ids"],
                }
            ]
        }
    }
    proj = FakeProjectCells()
    normalized = normalize_view_chart_params(view, proj)
    pruned, n = prune_view_charts_with_invalid_params(normalized, proj)
    assert n == 1
    assert pruned["initialCharts"]["cells"] == []


def test_normalize_resolves_color_by_display_name():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "s",
                    "type": "scatter_plot",
                    "param": ["n_genes", "total_counts"],
                    "color_by": "Leiden",
                }
            ]
        }
    }

    class P:
        datasources = [{"name": "cells"}]

        def get_datasource_metadata(self, name: str):
            return {
                "name": "cells",
                "columns": [
                    {"name": "Leiden", "field": "leiden", "datatype": "text"},
                    {"name": "n_genes", "field": "n_genes", "datatype": "integer"},
                    {"name": "total_counts", "field": "total_counts", "datatype": "double"},
                ],
            }

    out = normalize_view_chart_params(view, P())
    assert out["initialCharts"]["cells"][0]["color_by"] == "leiden"


def test_normalize_and_prune_keeps_scanpy_style_marker_table_params():
    """When HDF5/metadata use Scanpy column names, params stay on those field ids and prune keeps the chart."""
    view = {
        "initialCharts": {
            "chat_rank_genes_result": [
                {
                    "title": "Top DE genes",
                    "type": "table_chart",
                    "param": [
                        "group",
                        "names",
                        "logfoldchanges",
                        "pvals_adj",
                        "scores",
                    ],
                }
            ]
        }
    }
    proj = FakeProjectScanpyMarkerColumnNames()
    normalized = normalize_view_chart_params(view, proj)
    assert normalized["initialCharts"]["chat_rank_genes_result"][0]["param"] == [
        "group",
        "names",
        "logfoldchanges",
        "pvals_adj",
        "scores",
    ]
    pruned, n = prune_view_charts_with_invalid_params(normalized, proj)
    assert n == 0
    assert len(pruned["initialCharts"]["chat_rank_genes_result"]) == 1


def test_normalize_marker_alias_params_for_chat_rank_genes_result():
    view = {
        "initialCharts": {
            "chat_rank_genes_result": [
                {
                    "title": "markers",
                    "type": "table_chart",
                    "param": ["group", "names", "scores", "logfoldchanges", "pvals", "pvals_adj"],
                }
            ]
        }
    }
    proj = FakeProjectMarkerAliases()
    normalized = normalize_view_chart_params(view, proj)
    assert normalized["initialCharts"]["chat_rank_genes_result"][0]["param"] == [
        "cluster",
        "gene",
        "score",
        "logfoldchange",
        "pval",
        "pval_adj",
    ]


def test_prune_drops_unresolved_marker_param_after_alias_normalization():
    view = {
        "initialCharts": {
            "chat_rank_genes_result": [
                {
                    "title": "markers",
                    "type": "table_chart",
                    "param": ["group", "names", "unknown_token"],
                }
            ]
        }
    }
    proj = FakeProjectMarkerAliases()
    normalized = normalize_view_chart_params(view, proj)
    pruned, n = prune_view_charts_with_invalid_params(normalized, proj)
    assert n == 1
    assert pruned["initialCharts"]["chat_rank_genes_result"] == []


def test_client_needs_refresh_after_chat():
    from mdvtools.llm.chat_client_refresh import client_needs_refresh_after_chat

    assert client_needs_refresh_after_chat(
        "project.add_datasource('chat_rank_genes_result', df, replace_data=True)"
    )
    assert client_needs_refresh_after_chat("add_datasource_polars(name, lf)")
    assert client_needs_refresh_after_chat("delete_datasource('x')")
    assert not client_needs_refresh_after_chat("print('add_datasource')")
    assert not client_needs_refresh_after_chat("")
