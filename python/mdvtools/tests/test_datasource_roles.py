import inspect

from mdvtools.llm.datasource_roles import (
    CHAT_RANK_GENES_DATASOURCE_NAME,
    build_chatmdv_roles_constants_block,
    build_datasource_field_index,
    categorical_field_ids_from_metadata,
    column_field_id,
    find_datasources_for_fields,
    format_feature_table_field_policy,
    format_marker_gene_scanpy_fallback_policy,
    format_marker_ranking_viz_policy,
    format_metadata_column_schema_policy,
    format_mdv_first_data_access_policy,
    format_no_hallucination_chart_policy,
    format_obs_annadata_alignment_policy,
    format_obs_table_chart_param_policy,
    format_proportion_chart_policy,
    format_gene_signature_chart_policy,
    format_column_subset_completeness_policy,
    format_targeted_chart_policies,
    format_scanpy_hybrid_routing_policy,
    format_visualization_consistency_policy,
    infer_datasource_roles,
    is_multi_table_tabular_project,
    resolve_datasource_from_question,
)
from mdvtools.llm.dataset_scale import ProjectScale, assess_project_scale


class FakeProject:
    def __init__(self):
        self._names = ["cells", "rna", "protein"]

    def get_datasource_names(self):
        return self._names

    def get_datasource_metadata(self, name):
        if name == "cells":
            return {
                "name": "cells",
                "size": 100,
                "columns": [
                    {"field": "major", "name": "Major", "datatype": "text"},
                    {"field": "n_genes", "name": "n genes", "datatype": "integer"},
                ],
            }
        return {"name": name, "size": 0, "columns": []}

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


def test_infer_datasource_roles_skips_missing_datasource_on_link():
    class ProjectWithBadLink(FakeProject):
        def get_links(self, datasource, filter=None):
            links = super().get_links(datasource, filter)
            links = list(links)
            links.append(
                {
                    "datasource": None,
                    "link": {
                        "rows_as_columns": {
                            "name_column": "name",
                            "subgroups": {
                                "bogus_expr": {
                                    "name": "bogus_expr",
                                    "label": "bogus_expr",
                                    "type": "sparse",
                                },
                            },
                        }
                    },
                }
            )
            return links

    roles = infer_datasource_roles(ProjectWithBadLink())
    assert [e.datasource_name for e in roles.expressions] == ["rna", "protein"]
    assert "None" not in [e.datasource_name for e in roles.expressions]


def test_format_feature_table_field_policy_lists_name_columns():
    roles = infer_datasource_roles(FakeProject())
    text = format_feature_table_field_policy(roles)
    assert "name_column" in text
    assert "`rna`" in text or "rna" in text
    assert "Do not substitute" in text or "gene_ids" in text or "name_column" in text
    assert "Wrapper `param` tokens" in text
    assert "rna_expr" in text
    assert "rna_expr|" in text or "subgroup_key" in text


def test_format_no_hallucination_chart_policy():
    text = format_no_hallucination_chart_policy()
    assert "Metadata-first chart params" in text
    assert "Never** invent" in text or "invent a datasource" in text
    assert "color_by" in text
    assert "get_datasource_as_dataframe" in text
    assert "datatype" in text
    assert "dtype" in text


def test_categorical_field_ids_from_metadata():
    meta = {
        "columns": [
            {"field": "major", "datatype": "text", "name": "Major"},
            {"field": "X_umap_1", "datatype": "double", "name": "UMAP1"},
        ]
    }
    assert categorical_field_ids_from_metadata(meta) == ["major"]
    assert column_field_id(meta["columns"][0]) == "major"


def test_build_chatmdv_roles_constants_includes_categorical_field_ids():
    block = build_chatmdv_roles_constants_block(FakeProject())
    assert "CHATMDV_CATEGORICAL_FIELD_IDS" in block
    assert "major" in block


def test_format_metadata_column_schema_policy():
    text = format_metadata_column_schema_policy()
    assert "datatype" in text
    assert "dtype" in text
    assert "HeatmapPlot" in text
    assert "build_expression_wrapper_token" in text


def test_format_visualization_consistency_policy():
    text = format_visualization_consistency_policy()
    assert "Single source of truth" in text
    assert "Scanpy" in text or "AnnData" in text
    assert "wrapper" in text.lower()
    assert "rank_genes_groups" in text


def test_format_obs_table_chart_param_policy():
    doc = inspect.getdoc(format_obs_table_chart_param_policy)
    assert doc is not None
    assert "CHAT_RANK_GENES_DATASOURCE_NAME" in doc or "chat stdout" in doc.lower()
    assert "avoid persisting extra datasources" not in doc.lower()

    text = format_obs_table_chart_param_policy()
    assert "Field ID" in text
    assert "cells" in text
    assert "print" in text.lower() or "`print" in text
    assert CHAT_RANK_GENES_DATASOURCE_NAME in text
    assert "derive `params` from actual persisted fields" in text
    assert "Post-write marker check" in text


def test_format_marker_ranking_viz_policy():
    text = format_marker_ranking_viz_policy()
    assert "rank_genes_groups" in text
    assert "DotPlot" in text
    assert "Heatmap" in text
    assert "Interpretation / cell-type questions" in text
    assert CHAT_RANK_GENES_DATASOURCE_NAME in text
    assert "project.add_datasource(" in text
    assert "IndentationError" in text


def test_format_targeted_chart_policies():
    compact = format_targeted_chart_policies(compact=True)
    assert "StackedRowChart" in compact
    assert "stacked_row_plot" in compact
    assert "param=" in compact
    assert "obs_df.index" in compact or "AnnData" in compact

    full = format_targeted_chart_policies(compact=False)
    assert "Gene signature" in full or "signature" in full.lower()
    assert format_obs_annadata_alignment_policy()[:20] in full
    assert "PieChart" in format_proportion_chart_policy()


def test_format_mdv_first_includes_hybrid_bans():
    large = ProjectScale(
        obs_rows=500_000,
        obs_columns=40,
        estimated_obs_df_mb=300.0,
        available_ram_mb=2048.0,
        is_large=True,
        has_h5ad=True,
        obs_datasource="cells",
    )
    text = format_mdv_first_data_access_policy(large, "/data/x.h5ad", compact=True)
    assert "combine `read_h5ad`" in text
    assert "invent `bucket`" in text


def test_format_marker_gene_scanpy_fallback_with_h5ad():
    text = format_marker_gene_scanpy_fallback_policy("/data/adata.h5ad")
    assert "rank_genes_groups" in text
    assert "Never" in text or "never" in text.lower()
    assert "leiden" in text.lower()
    assert "table_chart" in text or "TablePlot" in text
    assert "**Primary** answer" in text or "Primary" in text
    assert CHAT_RANK_GENES_DATASOURCE_NAME in text
    assert "add_datasource" in text


def test_format_marker_gene_scanpy_fallback_large_project_requires_backed():
    large = ProjectScale(
        obs_rows=500_000,
        obs_columns=40,
        estimated_obs_df_mb=300.0,
        available_ram_mb=2048.0,
        is_large=True,
        has_h5ad=True,
        obs_datasource="cells",
    )
    text = format_marker_gene_scanpy_fallback_policy("/data/adata.h5ad", large)
    assert "backed='r'" in text
    assert ".copy()" in text or "copy" in text.lower()


def test_format_mdv_first_data_access_policy_large_project():
    large = ProjectScale(
        obs_rows=987_743,
        obs_columns=45,
        estimated_obs_df_mb=675.0,
        available_ram_mb=2700.0,
        is_large=True,
        has_h5ad=True,
        obs_datasource="cells",
    )
    text = format_mdv_first_data_access_policy(large, "/data/adata.h5ad", compact=False)
    assert "MDV-first" in text
    assert "backed='r'" in text
    assert "columns=[...]" in text
    assert "Scanpy last resort" in text
    assert "987,743" in text

    compact = format_mdv_first_data_access_policy(large, "/data/adata.h5ad", compact=True)
    assert "MDV-first data access" in compact
    assert "backed='r'" in compact
    assert "add_datasource" in compact


def test_assess_project_scale_marks_large_obs_table():
    scale = assess_project_scale(FakeProject(), "")
    assert scale.obs_datasource == "cells"
    assert scale.obs_rows == 100
    assert scale.is_large is False

    class LargeProject(FakeProject):
        def get_datasource_metadata(self, name):
            if name == "cells":
                return {
                    "name": "cells",
                    "size": 250_000,
                    "columns": [{"field": f"c{i}", "datatype": "text"} for i in range(10)],
                }
            return super().get_datasource_metadata(name)

    large = assess_project_scale(LargeProject(), "/data/x.h5ad")
    assert large.is_large is True
    assert large.has_h5ad is True


def test_format_marker_gene_scanpy_fallback_without_h5ad():
    text = format_marker_gene_scanpy_fallback_policy("")
    assert "Marker genes" in text or "marker" in text.lower()
    assert "MDV only" in text or "wrappers" in text.lower()
    assert "print" in text.lower()


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


def test_format_scanpy_hybrid_routing_policy():
    text = format_scanpy_hybrid_routing_policy()
    assert "Hybrid Scanpy routing contract" in text
    assert "1-row-per-observation" in text
    assert "cell_id" in text
    assert CHAT_RANK_GENES_DATASOURCE_NAME in text
    assert "text box" in text.lower()
    assert "dot / bubble / violin" in text.lower()
    assert "cluster distribution" in text.lower()
    assert "predicted_cell_type" in text
    assert "Marker table strictness" in text


class MicronLikeProject:
    """Minimal metadata shaped like MICRON-scanr-02 for routing tests."""

    datasources = [
        {"name": "qc_channels", "columns": [{"field": "channel_name"}, {"field": "run_id"}]},
        {
            "name": "qc_field_uniformity",
            "columns": [{"field": "channel_name"}, {"field": "cv_pct"}, {"field": "uniformity_score"}],
        },
        {"name": "qc_runs", "columns": [{"field": "assay"}, {"field": "run_id"}]},
        {"name": "qc_sessions", "columns": [{"field": "session_id"}]},
    ]

    dir = "/app/mdv/MICRON-scanr-02"

    def get_datasource_names(self):
        return [ds["name"] for ds in self.datasources]

    def get_datasource_metadata(self, name):
        for ds in self.datasources:
            if ds["name"] == name:
                return {"name": name, "size": 10, "columns": ds["columns"]}
        raise KeyError(name)

    def get_links(self, *_a, **_k):
        return []


def test_micron_resolve_datasource_from_question_by_name():
    project = MicronLikeProject()
    q = (
        "What is the distribution of cv_pct across different channel_name values "
        "in qc_field_uniformity?"
    )
    assert resolve_datasource_from_question(project, q) == "qc_field_uniformity"


def test_micron_resolve_datasource_from_question_by_fields():
    project = MicronLikeProject()
    q = "Show cv_pct and uniformity_score by channel"
    assert resolve_datasource_from_question(project, q) == "qc_field_uniformity"


def test_micron_resolve_datasource_qc_runs():
    project = MicronLikeProject()
    q = "What is the breakdown of assay types present in qc_runs?"
    assert resolve_datasource_from_question(project, q) == "qc_runs"


def test_resolve_datasource_n_genes_does_not_match_genes_substring():
    """``genes`` must not match inside the field token ``n_genes``."""

    class PbmcLikeProject:
        datasources = [{"name": "cells"}, {"name": "genes"}]

        def get_datasource_names(self):
            return ["cells", "genes"]

        def get_datasource_metadata(self, name):
            if name == "cells":
                return {
                    "name": "cells",
                    "size": 2700,
                    "columns": [
                        {"field": "leiden", "name": "leiden", "datatype": "text"},
                        {"field": "n_genes", "name": "n genes", "datatype": "integer"},
                    ],
                }
            return {
                "name": "genes",
                "size": 1838,
                "columns": [
                    {"field": "gene_ids", "name": "gene ids", "datatype": "text"},
                    {"field": "n_cells_by_counts", "name": "n cells", "datatype": "integer"},
                ],
            }

    project = PbmcLikeProject()
    q = "What is the distribution of n_genes across different leiden clusters?"
    assert resolve_datasource_from_question(project, q) == "cells"


def test_find_datasources_for_fields_maps_cv_pct():
    project = MicronLikeProject()
    index = build_datasource_field_index(project)
    owners = find_datasources_for_fields(index, ["cv_pct"])
    assert owners["cv_pct"] == ["qc_field_uniformity"]


def test_micron_infer_obs_datasource_prefers_qc_runs():
    project = MicronLikeProject()
    roles = infer_datasource_roles(project)
    assert roles.obs_datasource == "qc_runs"


def test_is_multi_table_tabular_project_for_micron():
    project = MicronLikeProject()
    roles = infer_datasource_roles(project)
    assert is_multi_table_tabular_project(project, roles) is True


def test_rag_prompt_includes_all_tables_and_resolved_datasource():
    from mdvtools.llm.templates import get_createproject_prompt_RAG

    project = MicronLikeProject()
    prompt = get_createproject_prompt_RAG(
        project,
        "",
        "qc_field_uniformity",
        'fields "channel_name", "cv_pct"\ncharts "Violin plot"',
        "What is the distribution of cv_pct in qc_field_uniformity?",
        compact=True,
    )
    assert "Primary datasource for this question: **qc_field_uniformity**" in prompt
    assert "qc_runs" in prompt
    assert "datasource_name = 'qc_field_uniformity'" in prompt

