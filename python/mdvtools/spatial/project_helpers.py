from __future__ import annotations

import json
import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from mdvtools.mdvproject import MDVProject


def build_spatial_regions_metadata(
    all_regions: dict[str, dict],
    region_field: str,
) -> dict:
    return {
        "position_fields": ["x", "y"],
        "region_field": region_field,
        "default_color": "x",
        "scale_unit": "µm",
        "scale": 1.0,
        "avivator": {
            "default_channels": [],
            "base_url": "spatial/",
        },
        "all_regions": all_regions,
    }


def set_default_spatial_image_view(
    mdv: "MDVProject",
    template_path: str,
    density: bool,
    emit,
) -> None:
    for region_name, region_data in mdv.get_datasource_metadata("cells")["regions"]["all_regions"].items():
        if "viv_image" in region_data:
            break
    else:
        raise ValueError("No region with a viv_image found")

    emit(f"Using region '{region_name}' for default view", verbose_only=True)
    with open(template_path, "r") as f:
        view_str = f.read().replace("<SPATIAL_REGION_NAME>", region_name)
        density_block = """
        {
            "linkedDsName": "genes", "maxItems": 15, "type": "RowsAsColsQuery"
        }
        """
        view_str = view_str.replace('"<DENSITY_FIELDS>"', density_block if density else "")
        mdv.set_view("default", json.loads(view_str), True)


def create_empty_spatial_mdv_project(
    output_folder: str,
    delete_existing: bool,
    region_field: str,
) -> "MDVProject":
    import pandas
    from mdvtools.mdvproject import MDVProject

    mdv = MDVProject(output_folder, delete_existing=delete_existing)

    cells_df = pandas.DataFrame(
        {
            "x": pandas.Series(dtype="float64"),
            "y": pandas.Series(dtype="float64"),
            region_field: pandas.Series(dtype="object"),
        }
    )
    mdv.add_datasource(
        "cells",
        cells_df,
        columns=[
            {"name": "x", "field": "x", "datatype": "double"},
            {"name": "y", "field": "y", "datatype": "double"},
            {"name": region_field, "field": region_field, "datatype": "text"},
        ],
        supplied_columns_only=True,
    )

    genes_df = pandas.DataFrame(
        {
            "name": pandas.Series(dtype="object"),
            "object_source": pandas.Series(dtype="object"),
        }
    )
    mdv.add_datasource(
        "genes",
        genes_df,
        columns=[
            {"name": "name", "field": "name", "datatype": "text"},
            {"name": "object_source", "field": "object_source", "datatype": "multitext", "delimiter": "|"},
        ],
        supplied_columns_only=True,
    )

    return mdv
