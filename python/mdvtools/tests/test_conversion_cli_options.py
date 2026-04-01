from click.testing import CliRunner

import mdvtools.cli as cli_module


def test_convert_scanpy_cli_forwards_x_umap_options(monkeypatch):
    captured = {}

    monkeypatch.setattr(cli_module.sc, "read_h5ad", lambda path: object())

    def fake_convert_scanpy_to_mdv(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs

    monkeypatch.setattr(cli_module, "convert_scanpy_to_mdv", fake_convert_scanpy_to_mdv)

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        [
            "convert-scanpy",
            "out",
            "input.h5ad",
            "--compute-x-umap",
            "--leiden-resolution",
            "0.75",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["kwargs"]["compute_x_umap"] is True
    assert captured["kwargs"]["leiden_resolution"] == 0.75


def test_convert_spatial_cli_forwards_x_umap_options(monkeypatch):
    captured = {}

    def fake_convert_spatialdata_to_mdv(args):
        captured["args"] = args

    monkeypatch.setattr(cli_module, "convert_spatialdata_to_mdv", fake_convert_spatialdata_to_mdv)

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        [
            "convert-spatial",
            "out",
            "spatial",
            "--compute-x-umap",
            "--leiden-resolution",
            "1.25",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["args"].compute_x_umap is True
    assert captured["args"].leiden_resolution == 1.25
