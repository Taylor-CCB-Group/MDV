"""Tests for the chunk shapes convert-spatial writes for images and labels.

Viv fetches one chunk of one channel per tile, so chunk size sets how many bytes each tile costs.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from mdvtools.spatial.conversion import SpatialDataConversionArgs, convert_spatialdata_to_mdv


def _write_source(path: Path, chunk_edge: int) -> None:
    import dask.array as da
    import spatialdata as sd
    from spatialdata.models import Image2DModel, Labels2DModel

    image = Image2DModel.parse(
        da.zeros((2, 64, 64), dtype=np.uint16),
        dims=("c", "y", "x"),
        chunks=(2, chunk_edge, chunk_edge),
        scale_factors=[2],
    )
    labels = Labels2DModel.parse(
        da.zeros((64, 64), dtype=np.uint32),
        dims=("y", "x"),
        chunks=(chunk_edge, chunk_edge),
    )
    sd.SpatialData(images={"img": image}, labels={"lab": labels}).write(str(path))


def _chunk_shape(array_dir: Path) -> list[int]:
    meta = json.loads((array_dir / "zarr.json").read_text())
    return meta["chunk_grid"]["configuration"]["chunk_shape"]


def _convert(tmp_path: Path, chunk_edge: int, **kwargs) -> Path:
    source = tmp_path / "source.zarr"
    _write_source(source, chunk_edge)
    temp_folder = tmp_path / "temp"
    temp_folder.mkdir()
    output = tmp_path / "project"
    convert_spatialdata_to_mdv(
        SpatialDataConversionArgs(
            spatialdata_path=str(source),
            output_folder=str(output),
            temp_folder=str(temp_folder),
            verbose=False,
            **kwargs,
        )
    )
    return output / "spatial" / "source.zarr"


@pytest.mark.filterwarnings("ignore")
def test_copy_caps_image_and_label_chunks(tmp_path):
    store = _convert(tmp_path, chunk_edge=64, image_chunk_size=16)

    assert _chunk_shape(store / "images" / "img" / "0") == [1, 16, 16]
    assert _chunk_shape(store / "images" / "img" / "1") == [1, 16, 16]
    assert _chunk_shape(store / "labels" / "lab" / "0") == [16, 16]


@pytest.mark.filterwarnings("ignore")
def test_copy_keeps_chunks_already_smaller_than_the_cap(tmp_path):
    store = _convert(tmp_path, chunk_edge=8, image_chunk_size=16)

    assert _chunk_shape(store / "images" / "img" / "0") == [1, 8, 8]
    assert _chunk_shape(store / "labels" / "lab" / "0") == [8, 8]


@pytest.mark.filterwarnings("ignore")
def test_zero_chunk_size_keeps_source_chunks(tmp_path):
    store = _convert(tmp_path, chunk_edge=64, image_chunk_size=0)

    assert _chunk_shape(store / "images" / "img" / "0") == [2, 64, 64]
    assert _chunk_shape(store / "labels" / "lab" / "0") == [64, 64]


def test_default_chunk_size_is_1024():
    args = SpatialDataConversionArgs(spatialdata_path="in", output_folder="out", temp_folder="tmp")

    assert args.image_chunk_size == 1024


@pytest.mark.filterwarnings("ignore")
def test_link_warns_about_large_chunks_and_leaves_source_unchanged(tmp_path, capsys):
    store = _convert(tmp_path, chunk_edge=64, image_chunk_size=16, link=True)

    assert "WARNING: --link keeps the source chunks" in capsys.readouterr().out
    assert _chunk_shape(store / "images" / "img" / "0") == [2, 64, 64]


@pytest.mark.filterwarnings("ignore")
def test_negative_chunk_size_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="0 or more"):
        _convert(tmp_path, chunk_edge=64, image_chunk_size=-1)
