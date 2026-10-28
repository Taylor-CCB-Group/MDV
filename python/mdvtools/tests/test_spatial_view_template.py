"""Tests for the default view that convert-spatial writes from spatial_view_template.json."""

import os

from mdvtools.spatial.project_helpers import (
    build_spatial_regions_metadata,
    create_empty_spatial_mdv_project,
    set_default_spatial_image_view,
)

import mdvtools.spatial.conversion as spatial_conversion

TEMPLATE_PATH = os.path.join(os.path.dirname(spatial_conversion.__file__), "spatial_view_template.json")


def _project_with_region(tmp_path, region: dict):
    mdv = create_empty_spatial_mdv_project(str(tmp_path), True, "spatial_region", "obs", "var")
    obs_md = mdv.get_datasource_metadata("obs")
    obs_md["regions"] = build_spatial_regions_metadata({"sample_global": region}, "spatial_region")
    mdv.set_datasource_metadata(obs_md)
    return mdv


def test_default_spatial_view_uses_spatialdata_chart(tmp_path):
    mdv = _project_with_region(
        tmp_path,
        {
            "spatial": {"coordinate_system": "global", "file": "sample.zarr"},
            "viv_image": {"file": "sample.zarr/images/morphology_focus"},
        },
    )

    set_default_spatial_image_view(mdv, TEMPLATE_PATH, False, lambda *a, **k: None, "obs", "var")

    view = mdv.get_view("default")
    assert view is not None
    [chart] = view["initialCharts"]["obs"]
    assert chart["type"] == "SpatialDataMdvRegionReact"
    assert chart["region"] == "sample_global"
    assert chart["background_filter"]["category"] == "sample_global"
    assert chart["color_by"]["linkedDsName"] == "var"
    assert "viv" not in chart
    assert "renderStack" not in chart
    assert set(view["initialCharts"]) == {"obs", "var"}
