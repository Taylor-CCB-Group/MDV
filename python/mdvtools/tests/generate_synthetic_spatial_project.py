#!/usr/bin/env python3
"""
Generate synthetic SpatialData-backed MDV projects for local testing.

The generated MDV project path is intentionally a direct child of ~/mdv by
default so it can be discovered by the existing project rescan flow.
"""

from __future__ import annotations

import argparse
import copy
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

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


def _parse_size_label(value: str) -> int:
    normalized = value.strip().lower().replace("_", "")
    if not normalized:
        raise argparse.ArgumentTypeError("cell count cannot be empty")
    multiplier = 1
    if normalized[-1] == "k":
        multiplier = 1_000
        normalized = normalized[:-1]
    elif normalized[-1] == "m":
        multiplier = 1_000_000
        normalized = normalized[:-1]
    try:
        count = int(float(normalized) * multiplier)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"invalid cell count: {value}") from error
    if count <= 0:
        raise argparse.ArgumentTypeError("cell counts must be positive")
    return count


def _parse_cell_count_list(value: str) -> list[int]:
    counts = [_parse_size_label(part) for part in value.split(",")]
    if not counts:
        raise argparse.ArgumentTypeError("at least one coordinate-system count is required")
    return counts


def _default_coordinate_system_cell_counts(n_cells: int, n_coordinate_systems: int) -> list[int]:
    if n_cells < n_coordinate_systems:
        raise SystemExit("n-cells must be at least n-coordinate-systems")
    if n_coordinate_systems == 1:
        return [n_cells]

    base = [
        1_000,
        5_000,
        10_000,
        25_000,
        50_000,
        100_000,
        250_000,
        500_000,
        1_000_000,
        2_000_000,
    ]
    if n_coordinate_systems <= len(base):
        counts = base[:n_coordinate_systems]
    else:
        counts = base + [2_000_000] * (n_coordinate_systems - len(base))

    base_total = sum(counts)
    if base_total == n_cells:
        return counts
    if base_total > n_cells:
        weights = np.asarray(counts, dtype=np.float64)
        scaled = np.maximum(1, np.floor(weights / weights.sum() * n_cells).astype(np.int64))
        diff = n_cells - int(scaled.sum())
        for index in range(abs(diff)):
            target = index % n_coordinate_systems
            if diff > 0:
                scaled[target] += 1
            elif scaled[target] > 1:
                scaled[target] -= 1
        return scaled.astype(int).tolist()

    max_per_coordinate_system = max(2_000_000, int(np.ceil(n_cells / n_coordinate_systems)))
    counts = counts.copy()
    remaining = n_cells - base_total
    index = len(counts) - 1
    while remaining > 0:
        capacity = max_per_coordinate_system - counts[index]
        if capacity > 0:
            added = min(capacity, remaining)
            counts[index] += added
            remaining -= added
        index = (index - 1) % len(counts)
    return counts


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


def _create_anndata(
    n_cells: int,
    n_genes: int,
    image_size: int,
    seed: int,
    coordinate_system_cell_counts: list[int],
):
    try:
        from dummy_spatialdata import generate_anndata
    except ImportError as error:
        raise SystemExit(
            "dummy-spatialdata is required. Install the Python dev dependencies "
            "or run `../venv/bin/pip install dummy-spatialdata==0.1.9`."
        ) from error

    rng = np.random.default_rng(seed)
    np.random.seed(seed)
    if n_cells > 1_000_000:
        try:
            import anndata as ad
        except ImportError as error:
            raise SystemExit(
                "anndata is required by dummy-spatialdata. Install the Python "
                "dev dependencies before generating large synthetic projects."
            ) from error
        adata = ad.AnnData(
            X=rng.random((n_cells, n_genes), dtype=np.float32),
            obs=pd.DataFrame(index=pd.RangeIndex(n_cells)),
            var=pd.DataFrame(index=pd.Index([f"var_{i}" for i in range(n_genes)])),
        )
    else:
        adata = generate_anndata(n_obs=n_cells, n_vars=n_genes)
        adata.X = np.asarray(adata.X, dtype=np.float32)
    adata.obsm["spatial"] = np.column_stack(
        [
            rng.uniform(0, image_size, n_cells),
            rng.uniform(0, image_size, n_cells),
        ]
    ).astype(np.float32, copy=False)
    cell_type_codes = np.arange(n_cells, dtype=np.int32) % 6
    sample_codes = np.arange(n_cells, dtype=np.int32) % 3
    adata.obs["cell_type"] = pd.Categorical.from_codes(
        cell_type_codes,
        categories=[f"type_{i}" for i in range(6)],
    )
    adata.obs["sample_id"] = pd.Categorical.from_codes(
        sample_codes,
        categories=[f"sample_{i}" for i in range(3)],
    )
    if len(coordinate_system_cell_counts) > 1:
        region_codes = np.repeat(
            np.arange(len(coordinate_system_cell_counts), dtype=np.int32),
            coordinate_system_cell_counts,
        )
        adata.obs["region"] = pd.Categorical.from_codes(
            region_codes,
            categories=[f"image_{i}" for i in range(len(coordinate_system_cell_counts))],
        )
    adata.obs["quality_score"] = rng.normal(0, 1, n_cells).astype(np.float32)
    adata.obs["total_counts"] = rng.poisson(2000, n_cells).astype(np.int32)
    return adata


def _create_spatialdata(
    *,
    n_cells: int,
    n_genes: int,
    image_size: int,
    seed: int,
    coordinate_system_cell_counts: list[int],
):
    try:
        import dummy_spatialdata as ds
    except ImportError as error:
        raise SystemExit(
            "dummy-spatialdata is required. Install the Python dev dependencies "
            "or run `../venv/bin/pip install dummy-spatialdata==0.1.9`."
        ) from error

    n_coordinate_systems = len(coordinate_system_cell_counts)
    adata = _create_anndata(
        n_cells,
        n_genes,
        image_size,
        seed,
        coordinate_system_cell_counts,
    )
    image_shape = {"x": image_size, "y": image_size}
    coordinate_system_names = [
        "global" if n_coordinate_systems == 1 else f"coordinate_system_{i}"
        for i in range(n_coordinate_systems)
    ]
    coordinate_systems = {
        name: {"transformations": ["identity"], "shape": image_shape}
        for name in coordinate_system_names
    }

    if n_cells <= 100_000 and n_coordinate_systems == 1:
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

    images = {
        f"image_{index}": ds.generate_imagemodel(
            {
                "type": "rgb",
                "shape": image_shape,
                "coordinate_system": [coordinate_system_names[index]],
            },
            key=f"image_{index}",
            coordinate_systems=coordinate_systems,
        )
        for index in range(n_coordinate_systems)
    }
    adata.obs["instance_id"] = adata.obs.index
    if n_coordinate_systems == 1:
        adata.obs["region"] = "image_0"
    regions = [f"image_{index}" for index in range(n_coordinate_systems)]
    table = TableModel.parse(
        adata,
        region=regions if n_coordinate_systems > 1 else "image_0",
        region_key="region",
        instance_key="instance_id",
    )
    images_for_sdata: Any = images
    tables: Any = {"table_0": table}
    return sd.SpatialData(images=images_for_sdata, tables=tables)


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
    row_chart = {
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
    }
    react_scatter = _scatter_config(
        chart_id="synth_scatter_react",
        chart_type="wgl_scatter_plot_dev",
        title="React scatter x/y",
        gsposition=[4, 0],
        image_size=image_size,
    )
    classic_scatter = _scatter_config(
        chart_id="synth_scatter_classic",
        chart_type="wgl_scatter_plot",
        title="Classic scatter x/y",
        gsposition=[8, 0],
        image_size=image_size,
    )
    react_table = {
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
    }
    classic_table = {
        "title": "Classic table all columns",
        "legend": "",
        "type": "table_chart",
        "param": table_columns,
        "id": "synth_table_classic",
        "size": [900, 500],
        "gssize": [6, 4],
        "gsposition": [6, 2],
        "column_widths": {},
    }

    def view_for(charts: list[dict]) -> dict:
        return {
            "dataSources": {
                obs_datasource_name: {"layout": "gridstack", "panelWidth": 100},
                var_datasource_name: {"layout": "gridstack", "panelWidth": 0},
            },
            "initialCharts": {
                obs_datasource_name: charts,
                var_datasource_name: [],
            },
        }

    comparison_view = view_for(
        [
            copy.deepcopy(row_chart),
            copy.deepcopy(react_scatter),
            copy.deepcopy(classic_scatter),
            copy.deepcopy(react_table),
            copy.deepcopy(classic_table),
        ]
    )
    react_view = view_for(
        [
            copy.deepcopy(row_chart),
            copy.deepcopy(react_table),
            copy.deepcopy(react_scatter),
        ]
    )
    classic_view = view_for(
        [
            copy.deepcopy(classic_scatter),
            copy.deepcopy(row_chart),
            copy.deepcopy(classic_table),
        ]
    )
    react_view["initialCharts"][obs_datasource_name][1]["gsposition"] = [4, 0]
    react_view["initialCharts"][obs_datasource_name][2]["gsposition"] = [0, 2]
    classic_view["initialCharts"][obs_datasource_name][0]["gsposition"] = [0, 2]
    classic_view["initialCharts"][obs_datasource_name][2]["gsposition"] = [4, 0]

    mdv.set_view("Scatter/table comparison", comparison_view)
    mdv.set_view("react view", react_view)
    mdv.set_view("classic view", classic_view)


def _make_safe_large_default(mdv) -> None:
    views = mdv.views
    state = mdv.state
    original_default = views.get("default")
    data_summary = views.get("Data summary")
    if original_default is not None:
        views["Spatial overview (loads all rows)"] = original_default
    if data_summary is not None:
        views["default"] = copy.deepcopy(data_summary)
    mdv.views = views

    all_views = state.get("all_views", [])
    if "Spatial overview (loads all rows)" not in all_views and original_default is not None:
        all_views.append("Spatial overview (loads all rows)")
    state["all_views"] = all_views
    state["initial_view"] = "default"
    mdv.state = state


def generate_project(
    *,
    output: Path,
    profile: str,
    n_cells: int,
    n_genes: int,
    image_size: int,
    seed: int,
    n_coordinate_systems: int,
    coordinate_system_cell_counts: list[int],
    force: bool,
) -> None:
    if n_cells != sum(coordinate_system_cell_counts):
        raise SystemExit("n-cells must match the sum of coordinate-system cell counts")
    if n_coordinate_systems != len(coordinate_system_cell_counts):
        raise SystemExit("n-coordinate-systems must match coordinate-system cell counts")
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
            coordinate_system_cell_counts=coordinate_system_cell_counts,
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
            if n_cells >= 5_000_000 and n_coordinate_systems == 1:
                _make_safe_large_default(mdv)
        state = mdv.state
        state["provenance"] = {
            "created_by": "mdvtools synthetic spatialdata generator",
            "generator": "dummy-spatialdata",
            "profile": profile,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "image_size": image_size,
            "seed": seed,
            "n_coordinate_systems": n_coordinate_systems,
            "coordinate_system_cell_counts": coordinate_system_cell_counts,
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
    parser.add_argument("--n-coordinate-systems", type=int, default=1)
    parser.add_argument(
        "--coordinate-system-cell-counts",
        type=_parse_cell_count_list,
        default=None,
        help=(
            "Comma-separated per-coordinate-system row counts, e.g. "
            "1k,5k,10k,25k,50k,100k,250k,500k,1m,2m. "
            "When supplied, --n-cells is ignored and inferred from the sum."
        ),
    )
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
    coordinate_system_cell_counts = args.coordinate_system_cell_counts
    n_cells = args.n_cells
    n_coordinate_systems = args.n_coordinate_systems
    if coordinate_system_cell_counts is not None:
        n_cells = sum(coordinate_system_cell_counts)
        n_coordinate_systems = len(coordinate_system_cell_counts)
    else:
        coordinate_system_cell_counts = _default_coordinate_system_cell_counts(
            n_cells,
            n_coordinate_systems,
        )

    output = args.output or _default_output(args.profile, n_cells)
    generate_project(
        output=output.expanduser(),
        profile=args.profile,
        n_cells=n_cells,
        n_genes=args.n_genes,
        image_size=args.image_size,
        seed=args.seed,
        n_coordinate_systems=n_coordinate_systems,
        coordinate_system_cell_counts=coordinate_system_cell_counts,
        force=args.force,
    )


if __name__ == "__main__":
    main()
