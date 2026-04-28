#!/usr/bin/env python3
"""
Generate synthetic SpatialData-backed MDV projects for local testing.

The generated MDV project path is intentionally a direct child of ~/mdv by
default so it can be discovered by the existing project rescan flow.
"""

from __future__ import annotations

import argparse
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from mdvtools.spatial.conversion import (
    SpatialDataConversionArgs,
    convert_spatialdata_to_mdv,
)


PROFILE_CHOICES = ("scatter-table", "spatial-overview")


def _size_label(n_cells: int) -> str:
    if n_cells >= 1_000_000 and n_cells % 1_000_000 == 0:
        return f"{n_cells // 1_000_000}m"
    if n_cells >= 1_000 and n_cells % 1_000 == 0:
        return f"{n_cells // 1_000}k"
    return str(n_cells)


def _default_output(profile: str, n_cells: int) -> Path:
    return Path.home() / "mdv" / f"synth-spatial--{profile}--{_size_label(n_cells)}"


def _configure_cache_dirs() -> None:
    temp_root = Path(tempfile.gettempdir())
    mpl_dir = temp_root / "mdv-matplotlib"
    cache_dir = temp_root / "mdv-cache"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))


def _create_anndata(n_cells: int, n_genes: int, image_size: int, seed: int):
    try:
        from dummy_spatialdata import generate_anndata
    except ImportError as error:
        raise SystemExit(
            "dummy-spatialdata is required. Install the Python dev dependencies "
            "or run `../venv/bin/pip install dummy-spatialdata==0.1.9`."
        ) from error

    rng = np.random.default_rng(seed)
    np.random.seed(seed)
    adata = generate_anndata(n_obs=n_cells, n_vars=n_genes)
    adata.X = np.asarray(adata.X, dtype=np.float32)
    adata.obsm["spatial"] = np.column_stack(
        [
            rng.uniform(0, image_size, n_cells),
            rng.uniform(0, image_size, n_cells),
        ]
    ).astype(np.float32, copy=False)
    adata.obs["cell_type"] = [f"type_{i % 6}" for i in range(n_cells)]
    adata.obs["sample_id"] = [f"sample_{i % 3}" for i in range(n_cells)]
    adata.obs["quality_score"] = rng.normal(0, 1, n_cells)
    adata.obs["total_counts"] = rng.poisson(2000, n_cells)
    return adata


def _create_spatialdata(
    *,
    n_cells: int,
    n_genes: int,
    image_size: int,
    seed: int,
):
    try:
        import dummy_spatialdata as ds
    except ImportError as error:
        raise SystemExit(
            "dummy-spatialdata is required. Install the Python dev dependencies "
            "or run `../venv/bin/pip install dummy-spatialdata==0.1.9`."
        ) from error

    adata = _create_anndata(n_cells, n_genes, image_size, seed)
    image_shape = {"x": image_size, "y": image_size}
    coordinate_systems = {
        "global": {"transformations": ["identity"], "shape": image_shape}
    }

    if n_cells <= 100_000:
        return ds.generate_dataset(
            images=[
                {
                    "type": "rgb",
                    "shape": image_shape,
                    "coordinate_system": ["global"],
                }
            ],
            shapes=[
                {
                    "n": n_cells,
                    "type": "circle",
                    "shape": image_shape,
                    "coordinate_system": ["global"],
                }
            ],
            tables=[{"table": adata, "element": "shape", "element_index": 0}],
            coordinate_systems=coordinate_systems,
            SEED=seed,
        )

    try:
        import spatialdata as sd
        from spatialdata.models import TableModel
    except ImportError as error:
        raise SystemExit(
            "spatialdata is required by dummy-spatialdata. Install the Python "
            "dev dependencies before generating synthetic SpatialData projects."
        ) from error

    image = ds.generate_imagemodel(
        {
            "type": "rgb",
            "shape": image_shape,
            "coordinate_system": ["global"],
        },
        key="image_0",
        coordinate_systems=coordinate_systems,
    )
    adata.obs["instance_id"] = adata.obs.index
    adata.obs["region"] = "image_0"
    table = TableModel.parse(
        adata,
        region="image_0",
        region_key="region",
        instance_key="instance_id",
    )
    images: Any = {"image_0": image}
    tables: Any = {"table_0": table}
    return sd.SpatialData(images=images, tables=tables)


def _table_columns(mdv, datasource: str) -> list[str]:
    metadata = mdv.get_datasource_metadata(datasource)
    return [
        column["field"]
        for column in metadata["columns"]
        if not column["field"].startswith("__")
    ]


def _scatter_config(
    *,
    chart_id: str,
    chart_type: str,
    title: str,
    gsposition: list[int],
    image_size: int,
) -> dict:
    if chart_type == "wgl_scatter_plot_dev":
        return {
            "course_radius": 1,
            "radius": 4,
            "opacity": 0.8,
            "tooltip": {"show": False},
            "category_filters": [],
            "zoom_on_filter": False,
            "point_shape": "circle",
            "category1": [],
            "category2": [],
            "contour_fill": False,
            "contour_bandwidth": 0.1,
            "contour_intensity": 1,
            "contour_opacity": 0.5,
            "contour_fillThreshold": 2,
            "densityFields": [],
            "dimension": "2d",
            "on_filter": "grey",
            "selectionFeatureCollection": {
                "type": "FeatureCollection",
                "features": [],
            },
            "field_legend": {"display": True},
            "axis": {
                "x": {"rotate_labels": False, "size": 20, "tickfont": 10},
                "y": {"rotate_labels": False, "size": 40, "tickfont": 10},
            },
            "viewState": {
                "target": [image_size / 2, image_size / 2, 0],
                "zoom": 0,
            },
            "title": title,
            "legend": "",
            "type": chart_type,
            "param": ["x", "y"],
            "id": chart_id,
            "size": [600, 360],
            "color_by": "cell_type",
            "gssize": [4, 3],
            "gsposition": gsposition,
        }
    return {
        "title": title,
        "legend": "",
        "type": chart_type,
        "param": ["x", "y"],
        "axis": {
            "x": {
                "size": 30,
                "label": "x",
                "textsize": 13,
                "textSize": 13,
                "tickfont": 10,
            },
            "y": {
                "size": 45,
                "label": "y",
                "textsize": 13,
                "textSize": 13,
                "tickfont": 10,
            },
        },
        "id": chart_id,
        "size": [600, 360],
        "tooltip": {"show": False},
        "default_color": "#377eb8",
        "brush": "poly",
        "on_filter": "hide",
        "radius": 2,
        "opacity": 0.8,
        "gssize": [4, 3],
        "gsposition": gsposition,
    }


def _add_scatter_table_comparison_view(
    *,
    mdv,
    image_size: int,
    obs_datasource_name: str,
    var_datasource_name: str,
) -> None:
    table_columns = _table_columns(mdv, obs_datasource_name)
    view = {
        "dataSources": {
            obs_datasource_name: {"layout": "gridstack", "panelWidth": 100},
            var_datasource_name: {"layout": "gridstack", "panelWidth": 0},
        },
        "initialCharts": {
            obs_datasource_name: [
                {
                    "title": "cell_type",
                    "legend": "",
                    "type": "row_chart",
                    "param": ["cell_type"],
                    "id": "synth_row_cell_type",
                    "size": [400, 300],
                    "axis": {
                        "x": {
                            "textSize": 13,
                            "label": "",
                            "size": 25,
                            "tickfont": 10,
                        }
                    },
                    "gssize": [4, 2],
                    "gsposition": [0, 0],
                },
                _scatter_config(
                    chart_id="synth_scatter_react",
                    chart_type="wgl_scatter_plot_dev",
                    title="React scatter x/y",
                    gsposition=[4, 0],
                    image_size=image_size,
                ),
                _scatter_config(
                    chart_id="synth_scatter_classic",
                    chart_type="wgl_scatter_plot",
                    title="Classic scatter x/y",
                    gsposition=[8, 0],
                    image_size=image_size,
                ),
                {
                    "title": "React table all columns",
                    "legend": "",
                    "type": "table_chart_react",
                    "param": table_columns,
                    "column_widths": {},
                    "order": {},
                    "include_index": False,
                    "sort": None,
                    "id": "synth_table_react",
                    "size": [900, 500],
                    "gssize": [6, 4],
                    "gsposition": [0, 2],
                },
                {
                    "title": "Classic table all columns",
                    "legend": "",
                    "type": "table_chart",
                    "param": table_columns,
                    "id": "synth_table_classic",
                    "size": [900, 500],
                    "gssize": [6, 4],
                    "gsposition": [6, 2],
                    "column_widths": {},
                },
            ],
            var_datasource_name: [],
        },
    }
    mdv.set_view("Scatter/table comparison", view)


def generate_project(
    *,
    output: Path,
    profile: str,
    n_cells: int,
    n_genes: int,
    image_size: int,
    seed: int,
    force: bool,
) -> None:
    if output.exists():
        if not force:
            raise SystemExit(f"{output} already exists. Pass --force to replace it.")
        shutil.rmtree(output)

    output.parent.mkdir(parents=True, exist_ok=True)
    _configure_cache_dirs()

    with tempfile.TemporaryDirectory(prefix="mdv-synth-spatial-") as temp_dir:
        temp_path = Path(temp_dir)
        source_path = temp_path / "source.zarr"
        conversion_temp = temp_path / "conversion"
        conversion_temp.mkdir()

        sdata = _create_spatialdata(
            n_cells=n_cells,
            n_genes=n_genes,
            image_size=image_size,
            seed=seed,
        )
        sdata.write(str(source_path))

        args = SpatialDataConversionArgs(
            spatialdata_path=str(source_path),
            output_folder=str(output),
            temp_folder=str(conversion_temp),
            preserve_existing=False,
            output_geojson=True,
            link=False,
            density=False,
            obs_datasource_name="cells",
            var_datasource_name="genes",
            verbose=False,
        )
        mdv = convert_spatialdata_to_mdv(args)
        if profile == "scatter-table":
            _add_scatter_table_comparison_view(
                mdv=mdv,
                image_size=image_size,
                obs_datasource_name=args.obs_datasource_name,
                var_datasource_name=args.var_datasource_name,
            )
        state = mdv.state
        state["provenance"] = {
            "created_by": "mdvtools synthetic spatialdata generator",
            "generator": "dummy-spatialdata",
            "profile": profile,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "image_size": image_size,
            "seed": seed,
            "cleanup_group": "synth-spatial",
        }
        mdv.state = state
        print(f"Created synthetic SpatialData MDV project: {output}")
        print("Run /rescan_projects or restart the DB-backed app to register it.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a synthetic SpatialData-backed MDV project.",
    )
    parser.add_argument(
        "--profile",
        choices=PROFILE_CHOICES,
        default="scatter-table",
        help="Named project profile to generate.",
    )
    parser.add_argument("--n-cells", type=int, default=1000)
    parser.add_argument("--n-genes", type=int, default=50)
    parser.add_argument("--image-size", type=int, default=512)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output MDV project path. Defaults to a flat ~/mdv/synth-spatial--... path.",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output or _default_output(args.profile, args.n_cells)
    generate_project(
        output=output.expanduser(),
        profile=args.profile,
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        image_size=args.image_size,
        seed=args.seed,
        force=args.force,
    )


if __name__ == "__main__":
    main()
