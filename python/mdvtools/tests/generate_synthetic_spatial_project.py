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
    adata = generate_anndata(n_obs=n_cells, n_vars=n_genes)
    adata.obsm["spatial"] = np.column_stack(
        [
            rng.uniform(0, image_size, n_cells),
            rng.uniform(0, image_size, n_cells),
        ]
    )
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
        coordinate_systems={
            "global": {"transformations": ["identity"], "shape": image_shape}
        },
        SEED=seed,
    )


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
