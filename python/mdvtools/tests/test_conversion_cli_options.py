from click.testing import CliRunner

import mdvtools.cli as cli_module
import mdvtools.conversions as conversions_module
import mdvtools.spatial.conversion as spatial_conversion_module
import scanpy as sc


def test_convert_scanpy_cli_forwards_x_umap_options(monkeypatch):
    captured = {}

    monkeypatch.setattr(sc, "read_h5ad", lambda path: object())

    def fake_convert_scanpy_to_mdv(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs

    monkeypatch.setattr(conversions_module, "convert_scanpy_to_mdv", fake_convert_scanpy_to_mdv)

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        [
            "convert-scanpy",
            "out",
            "input.h5ad",
            "--obs-datasource-name",
            "observations",
            "--var-datasource-name",
            "features",
            "--link-name-column",
            "display_name",
            "--compute-x-umap",
            "--leiden-resolution",
            "0.75",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["args"][5] == "observations"
    assert captured["args"][6] == "features"
    assert captured["kwargs"]["link_name_column"] == "display_name"
    assert captured["kwargs"]["compute_x_umap"] is True
    assert captured["kwargs"]["leiden_resolution"] == 0.75


def test_convert_spatial_cli_forwards_x_umap_options(monkeypatch):
    captured = {}

    def fake_convert_spatialdata_to_mdv(args):
        captured["args"] = args

    monkeypatch.setattr(spatial_conversion_module, "convert_spatialdata_to_mdv", fake_convert_spatialdata_to_mdv)

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        [
            "convert-spatial",
            "out",
            "spatial",
            "--obs-datasource-name",
            "observations",
            "--var-datasource-name",
            "features",
            "--link-name-column",
            "display_name",
            "--compute-x-umap",
            "--leiden-resolution",
            "1.25",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["args"].obs_datasource_name == "observations"
    assert captured["args"].var_datasource_name == "features"
    assert captured["args"].link_name_column == "display_name"
    assert captured["args"].compute_x_umap is True
    assert captured["args"].leiden_resolution == 1.25
