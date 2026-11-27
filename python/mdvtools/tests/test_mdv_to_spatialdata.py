"""MDV project → SpatialData stores, and the converted project with a markers datasource."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile
from click.testing import CliRunner

import mdvtools.cli as cli_module
from mdvtools.mdvproject import MDVProject
from mdvtools.spatial.mdv_to_spatialdata import (
    _read_segmentation,
    export_converted_project,
    export_mdv_to_spatialdata,
    expression_field,
)


def _write_tiff(path: Path, shape: tuple[int, int, int], channels: list[str], physical: float | None) -> None:
    metadata: dict = {"axes": "CYX", "Channel": {"Name": channels}}
    if physical is not None:
        metadata["PhysicalSizeX"] = physical
        metadata["PhysicalSizeY"] = physical
        metadata["PhysicalSizeXUnit"] = "µm"
        metadata["PhysicalSizeYUnit"] = "µm"
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, np.zeros(shape, dtype=np.uint16), ome=True, metadata=metadata)


def _write_geojson(path: Path, cells: list[tuple[int, float, float]]) -> None:
    features = []
    for cell_id, x, y in cells:
        ring = [[x, y], [x + 1, y], [x + 1, y + 1], [x, y + 1], [x, y]]
        features.append(
            {
                "type": "Feature",
                "properties": {"DN": cell_id},
                "geometry": {"type": "Polygon", "coordinates": [ring]},
            }
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))


def _cells() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Object Id": 1,
                "sample_id": "wide",
                "leiden": "a",
                "area_pixels": 5.0,
                "x": 1.0,
                "y": 2.0,
                "spatial_X_coord": 1.0,
                "spatial_Y_coord": 2.0,
                "X_pca_1": 0.1,
                "X_pca_2": 0.2,
                "log2_DAPI": 1.0,
                "log2_CD3": 2.0,
                "zscore_DAPI": 3.0,
                "zscore_CD3": 4.0,
            },
            {
                "Object Id": 2,
                "sample_id": "wide",
                "leiden": "b",
                "area_pixels": 6.0,
                "x": 3.0,
                "y": 4.0,
                "spatial_X_coord": 3.0,
                "spatial_Y_coord": 4.0,
                "X_pca_1": 0.3,
                "X_pca_2": 0.4,
                "log2_DAPI": 5.0,
                "log2_CD3": 6.0,
                "zscore_DAPI": 7.0,
                "zscore_CD3": 8.0,
            },
            {
                "Object Id": 5,
                "sample_id": "small",
                "leiden": "a",
                "area_pixels": 7.0,
                "x": 8.0,
                "y": 9.0,
                "spatial_X_coord": 8.0,
                "spatial_Y_coord": 9.0,
                "X_pca_1": 0.5,
                "X_pca_2": 0.6,
                "log2_DAPI": 9.0,
                "log2_CD3": 10.0,
                "zscore_DAPI": 11.0,
                "zscore_CD3": 12.0,
            },
        ]
    )


def _build_project(root: Path, *, geojson_ids: dict[str, list[int]] | None = None, marker_mismatch: bool = False) -> Path:
    project_dir = root / "project"
    project = MDVProject(str(project_dir), delete_existing=True)
    frame = _cells()
    if marker_mismatch:
        frame = frame.rename(columns={"zscore_CD3": "zscore_CD4"})
    project.add_datasource("cells", frame, add_to_view=None)
    project.set_column_group("cells", "X_pca", ["X_pca_2", "X_pca_1"])
    project.set_column_group("cells", "spatial", ["spatial_X_coord", "spatial_Y_coord"])
    project.set_column_group("cells", "log2", ["log2_DAPI", "log2_CD3"])
    if marker_mismatch:
        project.set_column_group("cells", "zscore", ["zscore_CD4", "zscore_DAPI"])
    else:
        project.set_column_group("cells", "zscore", ["zscore_CD3", "zscore_DAPI"])

    avivator = project_dir / "images" / "avivator"
    _write_tiff(avivator / "wide.ome.tiff", (2, 4, 1200), ["DAPI_FINAL", "CD3"], 0.325)
    _write_tiff(avivator / "small.ome.tiff", (2, 6, 8), ["DAPI_FINAL", "CD3"], None)
    ids = geojson_ids or {"wide": [1, 2], "small": [5]}
    _write_geojson(project_dir / "json" / "wide.json", [(cell_id, 1.0, 2.0) for cell_id in ids["wide"]])
    if "small" in ids:
        _write_geojson(project_dir / "json" / "small.json", [(cell_id, 8.0, 9.0) for cell_id in ids["small"]])

    metadata = project.get_datasource_metadata("cells")
    metadata["regions"] = {
        "position_fields": ["x", "y"],
        "region_field": "sample_id",
        "default_color": "leiden",
        "scale_unit": "µm",
        "scale": 0.5,
        "avivator": {"base_url": "images/avivator/", "default_channels": [{"name": "DAPI"}]},
        "all_regions": {
            "wide": {
                "roi": {"min_x": 0, "min_y": 0, "max_x": 1200, "max_y": 4},
                "images": {},
                "json": "json/wide.json",
                "viv_image": {"file": "wide.ome.tiff"},
            },
            "small": {
                "roi": {"min_x": 0, "min_y": 0, "max_x": 8, "max_y": 6},
                "images": {},
                "json": "json/small.json",
                "viv_image": {"file": "small.ome.tiff"},
            },
        },
    }
    project.set_datasource_metadata(metadata)
    project.set_view(
        "qc",
        {
            "initialCharts": {
                "cells": [
                    {
                        "type": "wgl_scatter_plot",
                        "param": ["log2_DAPI", "leiden"],
                        "color_by": "zscore_CD3",
                    }
                ]
            }
        },
        make_default=True,
    )
    return project_dir


def _level_shapes(image) -> list[tuple[int, int]]:
    from xarray import DataTree

    if isinstance(image, DataTree):
        return [
            (int(image[name]["image"].sizes["y"]), int(image[name]["image"].sizes["x"]))
            for name in image
        ]
    return [(int(image.sizes["y"]), int(image.sizes["x"]))]


def _channel_names(image) -> list[str]:
    from xarray import DataTree

    data = image["scale0"]["image"] if isinstance(image, DataTree) else image
    return [str(value) for value in data.coords["c"].values]


def _scale_xy(element) -> np.ndarray:
    from spatialdata.transformations import Scale, get_transformation

    transform = get_transformation(element, to_coordinate_system="global")
    assert isinstance(transform, Scale)
    return np.asarray(transform.to_scale_vector(("x", "y")), dtype=np.float64)


def test_export_writes_one_store_per_region_with_shared_names(tmp_path):
    project_dir = _build_project(tmp_path)
    output = tmp_path / "stores"
    written = export_mdv_to_spatialdata(
        project_dir,
        output,
        channel_renames={"DAPI_FINAL": "DAPI"},
    )

    assert [path.name for path in written] == ["wide.zarr", "small.zarr"]
    from spatialdata import read_zarr

    wide = read_zarr(output / "wide.zarr")
    small = read_zarr(output / "small.zarr")
    for store in (wide, small):
        assert list(store.images) == ["morphology"]
        assert list(store.shapes) == ["cell_boundaries"]
        assert list(store.tables) == ["table"]

    table = wide.tables["table"]
    assert list(table.var_names) == ["DAPI", "CD3"]
    np.testing.assert_allclose(np.asarray(table.X), np.asarray(table.layers["log2"]))
    np.testing.assert_allclose(np.asarray(table.layers["log2"]), [[1, 2], [5, 6]])
    np.testing.assert_allclose(np.asarray(table.layers["zscore"]), [[3, 4], [7, 8]])
    np.testing.assert_allclose(np.asarray(table.obsm["X_pca"]), [[0.1, 0.2], [0.3, 0.4]])
    np.testing.assert_allclose(np.asarray(table.obsm["spatial"]), [[1, 2], [3, 4]])
    assert list(table.obs["cell_id"].astype(str)) == ["1", "2"]
    assert list(table.obs["segmentation_parts"]) == [1, 1]
    assert set(table.obs["region"].astype(str)) == {"cell_boundaries"}
    assert list(table.obs["leiden"].astype(str)) == ["a", "b"]
    assert "log2_DAPI" not in table.obs.columns
    assert list(wide.shapes["cell_boundaries"].index.astype(str)) == ["1", "2"]
    assert _channel_names(wide.images["morphology"]) == ["DAPI", "CD3"]
    np.testing.assert_allclose(_scale_xy(wide.images["morphology"]), [0.325, 0.325])
    np.testing.assert_allclose(_scale_xy(wide.shapes["cell_boundaries"]), [0.325, 0.325])
    np.testing.assert_allclose(_scale_xy(small.images["morphology"]), [0.5, 0.5])

    wide_levels = _level_shapes(wide.images["morphology"])
    assert wide_levels[0] == (4, 1200)
    assert len(wide_levels) > 1
    assert max(wide_levels[-1]) <= 1024
    for previous, current in zip(wide_levels, wide_levels[1:]):
        assert current == (previous[0] // 2, previous[1] // 2)
    assert _level_shapes(small.images["morphology"]) == [(6, 8)]


def test_converted_project_uses_marker_link_fields_in_views(tmp_path):
    project_dir = _build_project(tmp_path)
    output = tmp_path / "converted"
    result = CliRunner().invoke(
        cli_module.cli,
        [
            "export-spatial",
            str(project_dir),
            str(output),
            "--as-project",
            "--rename-channel",
            "DAPI_FINAL=DAPI",
        ],
    )
    assert result.exit_code == 0, result.output

    project = MDVProject(str(output))
    cells = project.get_datasource_metadata("cells")
    fields = [column["field"] for column in cells["columns"]]
    assert "log2_DAPI" not in fields
    assert "leiden" in fields
    assert "x" in fields
    group_names = [group["name"] for group in cells["columnGroups"]]
    assert group_names == ["X_pca", "spatial"]
    assert cells["regions"]["all_regions"]["wide"]["spatial"] == {
        "coordinate_system": "global",
        "file": "wide.zarr",
    }
    assert (output / "spatial" / "wide.zarr").is_dir()
    assert (output / "images" / "avivator" / "wide.ome.tiff").is_file()
    assert (output / "json" / "wide.json").is_file()

    markers = project.get_datasource_as_dataframe("markers")
    assert markers["name"].tolist() == ["DAPI", "CD3"]
    links = project.get_links("cells", "rows_as_columns")
    assert len(links) == 1
    assert links[0]["datasource"] == "markers"
    subgroups = links[0]["link"]["rows_as_columns"]["subgroups"]
    assert set(subgroups) == {"log2", "zscore"}

    chart = project.get_view("qc")["initialCharts"]["cells"][0]
    assert chart["param"] == [expression_field("log2", "DAPI", 0), "leiden"]
    assert chart["color_by"] == expression_field("zscore", "CD3", 1)
    assert expression_field("log2", "DAPI", 0) == "log2|DAPI (log2)|0"


def test_marker_groups_must_share_a_marker_set(tmp_path):
    project_dir = _build_project(tmp_path, marker_mismatch=True)
    with pytest.raises(ValueError, match="do not share a marker set"):
        export_mdv_to_spatialdata(project_dir, tmp_path / "stores")


def test_repeated_dn_is_one_multipolygon(tmp_path):
    path = tmp_path / "cells.json"
    _write_geojson(path, [(1, 1.0, 2.0), (1, 4.0, 5.0), (2, 3.0, 4.0)])
    polygons = _read_segmentation(path)
    first = polygons[1]
    second = polygons[2]
    assert first.parts == 2
    assert getattr(first.geometry, "geom_type", None) == "MultiPolygon"
    assert len(list(getattr(first.geometry, "geoms"))) == 2
    assert second.parts == 1
    assert getattr(second.geometry, "geom_type", None) == "Polygon"


def test_segmentation_ids_must_match_object_ids(tmp_path):
    project_dir = _build_project(tmp_path, geojson_ids={"wide": [1, 99], "small": [5]})
    with pytest.raises(ValueError, match="GeoJSON DN and Object Id differ"):
        export_mdv_to_spatialdata(project_dir, tmp_path / "stores")


def test_missing_segmentation_path_is_not_replaced(tmp_path):
    project_dir = _build_project(tmp_path)
    metadata = MDVProject(str(project_dir)).get_datasource_metadata("cells")
    metadata["regions"]["all_regions"]["wide"]["json"] = "json/missing.json"
    MDVProject(str(project_dir)).set_datasource_metadata(metadata)
    with pytest.raises(FileNotFoundError, match="missing.json"):
        export_mdv_to_spatialdata(project_dir, tmp_path / "stores")


def test_export_converted_project_refuses_a_nonempty_directory(tmp_path):
    project_dir = _build_project(tmp_path)
    output = tmp_path / "converted"
    output.mkdir()
    (output / "existing").write_text("keep")
    with pytest.raises(FileExistsError):
        export_converted_project(project_dir, output)
