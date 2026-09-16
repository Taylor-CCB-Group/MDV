import json
from pathlib import Path

from click.testing import CliRunner

import mdvtools.cli as cli_module
from mdvtools.template_views import (
    VIEW_NAMES,
    create_default_views,
    infer_roles,
    unique_top_markers,
)
from mdvtools.template_views.create_default_views import Grid, MarkerResult, pretty_field


def _tiny_datasources() -> list[dict]:
    return [
        {
            "name": "cells",
            "size": 10,
            "columns": [
                {
                    "field": "annotation",
                    "name": "annotation",
                    "datatype": "text",
                    "values": ["T", "B"],
                },
                {
                    "field": "tissue",
                    "name": "tissue",
                    "datatype": "text",
                    "values": ["blood", "skin"],
                },
                {
                    "field": "diagnosis",
                    "name": "diagnosis",
                    "datatype": "text",
                    "values": ["healthy", "disease"],
                },
                {
                    "field": "X_umap_1",
                    "name": "X_umap_1",
                    "datatype": "double",
                    "minMax": [0.0, 1.0],
                },
                {
                    "field": "X_umap_2",
                    "name": "X_umap_2",
                    "datatype": "double",
                    "minMax": [0.0, 1.0],
                },
                {
                    "field": "n_genes_by_counts",
                    "name": "n_genes_by_counts",
                    "datatype": "double",
                    "minMax": [100.0, 2000.0],
                },
            ],
            "links": {
                "genes": {
                    "rows_as_columns": {
                        "name_column": "name",
                        "subgroups": {"gs": {"name": "gs"}},
                    }
                }
            },
        },
        {
            "name": "genes",
            "size": 3,
            "columns": [
                {
                    "field": "name",
                    "name": "name",
                    "datatype": "text",
                    "values": ["CD14", "CD3D", "MS4A1"],
                }
            ],
        },
    ]


def _tiny_project(tmp_path: Path) -> Path:
    project = tmp_path / "project"
    project.mkdir()
    (project / "datasources.json").write_text(json.dumps(_tiny_datasources(), indent=2) + "\n")
    (project / "views.json").write_text("{}\n")
    (project / "state.json").write_text(json.dumps({"all_views": ["default"]}) + "\n")
    return project


def test_infer_roles_picks_obs_cell_type_and_embedding():
    roles = infer_roles(_tiny_datasources())
    assert roles.obs == "cells"
    assert roles.cell_type == "annotation"
    assert roles.tissue == "tissue"
    assert roles.disease == "diagnosis"
    assert roles.embedding == ["X_umap_1", "X_umap_2"]
    assert roles.rna is not None
    assert roles.rna.ds_name == "genes"
    assert roles.rna.names == ["CD14", "CD3D", "MS4A1"]
    assert "CD14" in roles.gene_wrappers


def test_pretty_field_rewrites_generated_cluster_ids():
    assert (
        pretty_field("RNA_nbclust_5a743e84.988e.47bf.ada7.33b307458379_1_clusters")
        == "Clusters"
    )
    assert pretty_field("rna:RNA_nbclust_abc_1_clusters") == "Clusters"
    assert pretty_field("spatialclust_region_assignments") == "Spatial niche"
    assert pretty_field("tissue") == "tissue"
    assert pretty_field("rna:annotation") == "annotation"


def test_grid_wraps_at_twelve_columns():
    grid = Grid()
    first, size_a = grid.place(6, 5)
    second, size_b = grid.place(6, 5)
    third, size_c = grid.place(6, 4)
    assert first == [0, 0] and size_a == [6, 5]
    assert second == [6, 0] and size_b == [6, 5]
    assert third == [0, 5] and size_c == [6, 4]


def test_unique_top_markers_first_cluster_wins_and_caps():
    rows = [
        {"contrast": "cluster vs rest", "cluster": "A", "gene": "CD14", "rank": 1},
        {"contrast": "cluster vs rest", "cluster": "A", "gene": "LYZ", "rank": 2},
        {"contrast": "cluster vs rest", "cluster": "B", "gene": "CD14", "rank": 1},
        {"contrast": "cluster vs rest", "cluster": "B", "gene": "MS4A1", "rank": 2},
        {"contrast": "treatment", "cluster": "within cell type", "gene": "TNF", "rank": 1},
    ]
    picked = unique_top_markers(MarkerResult(rows=rows), per_cluster=2, cap=2)
    assert picked == [("CD14", "A"), ("LYZ", "A")]


def test_create_default_views_writes_expected_names(tmp_path):
    project = _tiny_project(tmp_path)
    summary = create_default_views(project, skip_markers=True)

    assert "1 Study overview" in summary["added"]
    assert "2 Cell atlas" in summary["added"]
    assert "3 Marker genes and signatures" in summary["added"]
    assert "4 Tissue and disease context" in summary["added"]
    assert "5 RNA-protein concordance" not in summary["added"]

    views = json.loads((project / "views.json").read_text())
    state = json.loads((project / "state.json").read_text())
    for name in summary["added"]:
        assert name in views
        assert name in state["all_views"]
    assert "default" in state["all_views"]
    assert (project / "views.json.bak").is_file()
    assert set(VIEW_NAMES) >= set(summary["added"])


def test_create_default_views_cli_forwards_skip_markers(tmp_path, monkeypatch):
    captured = {}

    def fake_create_default_views(project, *, skip_markers=False, recompute_markers=False):
        captured["project"] = project
        captured["skip_markers"] = skip_markers
        captured["recompute_markers"] = recompute_markers

    import mdvtools.template_views as template_views_pkg

    monkeypatch.setattr(template_views_pkg, "create_default_views", fake_create_default_views)

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        ["create-default-views", str(tmp_path / "proj"), "--skip-markers"],
    )
    assert result.exit_code == 0, result.output
    assert captured["skip_markers"] is True
    assert captured["recompute_markers"] is False
