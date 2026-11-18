import json
from pathlib import Path

import pytest
from click.testing import CliRunner

import mdvtools.cli as cli_module
import mdvtools.spatial.conversion as spatial_conversion_module
from mdvtools.spatial.conversion import SpatialDataConversionArgs, convert_spatialdata_to_mdv


def _write_store(path: Path) -> None:
    # Built by hand rather than with spatialdata.datasets.blobs, which runs numba code
    # that fails when another test module has set NUMBA_DISABLE_JIT.
    import dask.array as da
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from shapely.geometry import Point
    from spatialdata import SpatialData
    from spatialdata.models import Image2DModel, ShapesModel, TableModel

    rng = np.random.default_rng(0)
    pixels = da.from_array(rng.integers(0, 255, size=(1, 64, 64)).astype(np.uint8))
    centres = rng.uniform(0, 64, size=(10, 2))
    circles = gpd.GeoDataFrame(
        {"radius": np.full(10, 2.0)},
        geometry=[Point(x, y) for x, y in centres],
    )
    table = AnnData(
        X=rng.random((10, 3)).astype(np.float32),
        obs=pd.DataFrame(
            {"region": pd.Categorical(["cells"] * 10), "instance_id": np.arange(10)},
            index=pd.Index([str(i) for i in range(10)]),
        ),
    )
    table.obsm["spatial"] = centres
    SpatialData(
        images={"image": Image2DModel.parse(pixels, dims=("c", "y", "x"))},
        shapes={"cells": ShapesModel.parse(circles)},
        tables={
            "table": TableModel.parse(table, region="cells", region_key="region", instance_key="instance_id")
        },
    ).write(path)


def _region_spatial_entries(project: Path) -> list[dict]:
    datasources = json.loads((project / "datasources.json").read_text())
    return [
        region["spatial"]
        for datasource in datasources
        for region in datasource.get("regions", {}).get("all_regions", {}).values()
        if "spatial" in region
    ]


def _convert(tmp_path: Path, source: Path, **kwargs) -> Path:
    project = tmp_path / "project"
    temp = tmp_path / "temp"
    temp.mkdir()
    convert_spatialdata_to_mdv(
        SpatialDataConversionArgs(
            spatialdata_path=str(source),
            output_folder=str(project),
            temp_folder=str(temp),
            **kwargs,
        )
    )
    return project


def test_single_store_regions_read_from_the_url_and_nothing_is_copied(tmp_path):
    source = tmp_path / "sample.zarr"
    _write_store(source)

    project = _convert(tmp_path, source, spatial_url="http://127.0.0.1:8091/")

    spatial = _region_spatial_entries(project)
    assert spatial
    assert all(entry["url"] == "http://127.0.0.1:8091/" for entry in spatial)
    assert all(entry["file"] == "sample.zarr" for entry in spatial)
    assert not (project / "spatial").exists()


def test_regions_have_no_url_without_the_option(tmp_path):
    source = tmp_path / "sample.zarr"
    _write_store(source)

    project = _convert(tmp_path, source)

    spatial = _region_spatial_entries(project)
    assert spatial
    assert all("url" not in entry for entry in spatial)
    assert (project / "spatial" / "sample.zarr").is_dir()


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"spatial_url": " "}, "must not be empty"),
        ({"spatial_url": "http://x/", "link": True}, "cannot be combined"),
        ({"spatial_url": "http://x/", "batch": True}, "cannot be combined with --batch"),
    ],
)
def test_conversion_rejects_invalid_spatial_url_before_reading(tmp_path, kwargs, message):
    args = SpatialDataConversionArgs(
        spatialdata_path=str(tmp_path / "missing.zarr"),
        output_folder=str(tmp_path / "project"),
        temp_folder=str(tmp_path),
        **kwargs,
    )
    with pytest.raises(ValueError, match=message):
        convert_spatialdata_to_mdv(args)


def test_cli_forwards_spatial_url(monkeypatch):
    captured = {}

    def fake_convert_spatialdata_to_mdv(args):
        captured["args"] = args

    monkeypatch.setattr(spatial_conversion_module, "convert_spatialdata_to_mdv", fake_convert_spatialdata_to_mdv)

    result = CliRunner().invoke(
        cli_module.cli,
        ["convert-spatial", "out", "spatial", "--spatial-url", "http://127.0.0.1:8091/"],
    )

    assert result.exit_code == 0, result.output
    assert captured["args"].spatial_url == "http://127.0.0.1:8091/"
    assert captured["args"].link is False


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        (["--link"], "cannot be combined"),
        (["--batch"], "cannot be combined with --batch"),
        (["--csc-tables"], "cannot be combined with --spatial-url"),
    ],
)
def test_cli_rejects_spatial_url_with_batch_link_or_rewrites(monkeypatch, extra, message):
    monkeypatch.setattr(
        spatial_conversion_module,
        "convert_spatialdata_to_mdv",
        lambda args: pytest.fail("conversion should not run"),
    )

    result = CliRunner().invoke(
        cli_module.cli,
        ["convert-spatial", "out", "spatial", "--spatial-url", "http://x/", *extra],
    )

    assert result.exit_code != 0
    assert message in result.output
