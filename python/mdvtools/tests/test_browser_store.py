import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

import mdvtools.cli as cli_module
import mdvtools.spatial.conversion as spatial_conversion_module
from mdvtools.spatial.browser_store import (
    BROWSER_PRESETS,
    BrowserStoreOptions,
    resolve_browser_store_options,
    write_spatialdata_store,
)
from mdvtools.spatial.conversion import SpatialDataConversionArgs, convert_spatialdata_to_mdv
from mdvtools.spatial.conversion_options import BROWSER_PRESET_CHOICES


def test_presets_match_cli_choices():
    assert tuple(BROWSER_PRESETS) == BROWSER_PRESET_CHOICES


def test_default_preset_rewrites_nothing():
    options = resolve_browser_store_options()
    assert options == BrowserStoreOptions()
    assert not options.enabled
    assert options.describe() == "none"


def test_standard_and_experimental_presets():
    standard = resolve_browser_store_options("standard")
    assert standard == BrowserStoreOptions(
        image_codec="htj2k",
        image_preset="balanced",
        csc_tables=True,
        morton_points=True,
        points_zarr=False,
    )
    assert resolve_browser_store_options("experimental").points_zarr is True


def test_overrides_apply_on_top_of_preset():
    options = resolve_browser_store_options(
        "experimental",
        image_codec="none",
        csc_tables=False,
        image_preset="small",
    )
    assert options.image_codec is None
    assert options.image_preset == "small"
    assert options.csc_tables is False
    assert options.morton_points is True
    assert options.points_zarr is True

    assert resolve_browser_store_options(morton_points=True).enabled


def test_invalid_values_are_rejected():
    with pytest.raises(ValueError, match="browser preset"):
        resolve_browser_store_options("fastest")
    with pytest.raises(ValueError, match="image codec"):
        BrowserStoreOptions(image_codec="png")
    with pytest.raises(ValueError, match="image preset"):
        BrowserStoreOptions(image_codec="htj2k", image_preset="tiny")


def test_cli_forwards_browser_store_options(monkeypatch):
    captured = {}

    def fake_convert_spatialdata_to_mdv(args):
        captured["args"] = args

    monkeypatch.setattr(spatial_conversion_module, "convert_spatialdata_to_mdv", fake_convert_spatialdata_to_mdv)

    result = CliRunner().invoke(
        cli_module.cli,
        ["convert-spatial", "out", "spatial", "--browser-preset", "standard", "--no-morton-points", "--points-zarr"],
    )

    assert result.exit_code == 0, result.output
    assert captured["args"].browser_store == BrowserStoreOptions(
        image_codec="htj2k",
        image_preset="balanced",
        csc_tables=True,
        morton_points=False,
        points_zarr=True,
    )


def test_cli_defaults_to_no_browser_store_rewrites(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        spatial_conversion_module,
        "convert_spatialdata_to_mdv",
        lambda args: captured.setdefault("args", args),
    )

    result = CliRunner().invoke(cli_module.cli, ["convert-spatial", "out", "spatial"])

    assert result.exit_code == 0, result.output
    assert not captured["args"].browser_store.enabled


def test_cli_rejects_browser_store_rewrites_with_link(monkeypatch):
    monkeypatch.setattr(
        spatial_conversion_module,
        "convert_spatialdata_to_mdv",
        lambda args: pytest.fail("conversion should not run"),
    )

    result = CliRunner().invoke(
        cli_module.cli,
        ["convert-spatial", "out", "spatial", "--link", "--csc-tables"],
    )

    assert result.exit_code != 0
    assert "cannot be combined with --link" in result.output


def test_conversion_rejects_browser_store_rewrites_with_link_before_reading(tmp_path):
    args = SpatialDataConversionArgs(
        spatialdata_path=str(tmp_path / "missing.zarr"),
        output_folder=str(tmp_path / "out"),
        temp_folder=str(tmp_path),
        link=True,
        browser_store=BrowserStoreOptions(morton_points=True),
    )
    with pytest.raises(ValueError, match="cannot be combined with --link"):
        convert_spatialdata_to_mdv(args)


def _zarr_json(path: Path) -> dict:
    return json.loads((path / "zarr.json").read_text())


def _first_level_codecs(raster_group: Path) -> list[str]:
    level = _zarr_json(raster_group)["attributes"]["ome"]["multiscales"][0]["datasets"][0]["path"]
    return [codec["name"] for codec in _zarr_json(raster_group / level)["codecs"]]


@pytest.fixture
def small_sdata():
    # Built by hand rather than with spatialdata.datasets.blobs, which runs numba code
    # that fails when another test module has set NUMBA_DISABLE_JIT.
    import dask.array as da
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from anndata import AnnData
    from spatialdata import SpatialData
    from spatialdata.models import Image2DModel, PointsModel, TableModel

    rng = np.random.default_rng(0)
    pixels = da.from_array(rng.integers(0, 255, size=(3, 256, 256)))
    counts = rng.random((20, 5))
    counts[counts < 0.7] = 0
    points = pd.DataFrame(
        {
            "x": rng.uniform(0, 256, 500),
            "y": rng.uniform(0, 256, 500),
            "genes": pd.Categorical(rng.choice(["a", "b", "c"], 500)),
        }
    )
    table = AnnData(
        X=sp.csr_matrix(counts.astype(np.float32)),
        obs=pd.DataFrame(index=pd.Index([str(i) for i in range(20)])),
    )
    return SpatialData(
        images={
            "blobs_uint8": Image2DModel.parse(pixels.astype(np.uint8), dims=("c", "y", "x")),
            # Browser codecs can't take float pixels.
            "blobs_image": Image2DModel.parse(pixels.astype(np.float64), dims=("c", "y", "x")),
        },
        points={"blobs_points": PointsModel.parse(points, feature_key="genes")},
        tables={"table": TableModel.parse(table)},
    )


def test_write_without_rewrites_matches_plain_write(small_sdata, tmp_path):
    import spatialdata as sd

    dest = tmp_path / "plain.zarr"
    write_spatialdata_store(small_sdata, dest, BrowserStoreOptions(), progress=lambda message: None)

    assert "experimental.openjph_htj2k" not in _first_level_codecs(dest / "images" / "blobs_uint8")
    assert not (dest / "points" / "blobs_points" / "experimental.points-zarr").exists()
    assert set(sd.read_zarr(dest).images) == set(small_sdata.images)


def test_write_applies_experimental_rewrites(small_sdata, tmp_path):
    import spatialdata as sd
    from spatialdata_js_util import htj2k_available

    if not htj2k_available():
        pytest.skip("no HTJ2K encoder backend available")

    dest = tmp_path / "browser.zarr"
    messages: list[str] = []
    write_spatialdata_store(small_sdata, dest, BROWSER_PRESETS["experimental"], progress=messages.append)

    assert "experimental.openjph_htj2k" in _first_level_codecs(dest / "images" / "blobs_uint8")
    # Float images are left as spatialdata wrote them rather than failing the write.
    assert "experimental.openjph_htj2k" not in _first_level_codecs(dest / "images" / "blobs_image")
    assert any("Not recompressing images/blobs_image" in message for message in messages)
    assert _zarr_json(dest / "tables" / "table" / "X")["attributes"]["encoding-type"] == "csc_matrix"

    points_element = dest / "points" / "blobs_points"
    assert (points_element / "experimental.points-zarr" / "zarr.json").is_file()
    assert "morton_code_2d" in pq.read_schema(points_element / "points.parquet").names

    # The rewritten store is still a valid SpatialData store.
    sdata = sd.read_zarr(dest)
    assert set(sdata.images) == set(small_sdata.images)
    assert set(sdata.points) == set(small_sdata.points)
    assert sdata.tables["table"].n_obs == small_sdata.tables["table"].n_obs
    assert any("Morton-sorting" in message for message in messages)


def test_write_skips_points_zarr_for_points_without_feature_key(small_sdata, tmp_path):
    from spatialdata.models import PointsModel

    points = small_sdata.points["blobs_points"].compute()
    small_sdata.points["blobs_points"] = PointsModel.parse(points[["x", "y"]])
    dest = tmp_path / "no_feature.zarr"
    messages: list[str] = []
    write_spatialdata_store(small_sdata, dest, BrowserStoreOptions(points_zarr=True), progress=messages.append)

    assert "No points-zarr copy of points/blobs_points: it declares no feature_key" in messages
    assert not (dest / "points" / "blobs_points" / "experimental.points-zarr").exists()


def test_write_keeps_points_that_are_already_morton_indexed(small_sdata, tmp_path):
    first = tmp_path / "first.zarr"
    write_spatialdata_store(small_sdata, first, BrowserStoreOptions(morton_points=True), progress=lambda message: None)
    parquet = first / "points" / "blobs_points" / "points.parquet"
    sorted_once = pq.read_table(parquet)

    import spatialdata as sd

    second = tmp_path / "second.zarr"
    messages: list[str] = []
    write_spatialdata_store(
        sd.read_zarr(first),
        second,
        BrowserStoreOptions(morton_points=True),
        progress=messages.append,
    )

    assert "Keeping points/blobs_points: already Morton-indexed" in messages
    resorted = pq.read_table(second / "points" / "blobs_points" / "points.parquet")
    assert resorted.column("morton_code_2d").equals(sorted_once.column("morton_code_2d"))
