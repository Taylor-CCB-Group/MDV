#!/usr/bin/env python3
"""
Generate synthetic AnnData-backed MDV projects for local testing.

Unlike ``generate_synthetic_spatial_project`` (SpatialData / dummy-spatialdata),
this entrypoint uses ``mdvtools.tests.mock_anndata`` plus ``convert_scanpy_to_mdv``.
Typical uses: chart/UI testing and imports without SpatialData dependencies.

Output defaults to a direct child of ``~/mdv`` so ``GET /rescan_projects`` can
register it (same flat-folder convention as other synthetic generators).
"""

from __future__ import annotations

import argparse
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import scanpy as sc

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.tests.mock_anndata import (
    MockAnnDataFactory,
    create_minimal_anndata,
)
from mdvtools.tests.mock_expression_layers import add_synth_expression_layers

PROFILE_CHOICES = ("minimal", "realistic", "large", "memory-efficient")


def _size_label(n_cells: int) -> str:
    if n_cells >= 1_000_000 and n_cells % 1_000_000 == 0:
        return f"{n_cells // 1_000_000}m"
    if n_cells >= 1_000 and n_cells % 1_000 == 0:
        return f"{n_cells // 1_000}k"
    return str(n_cells)


def _default_output(profile: str, n_cells: int, n_genes: int) -> Path:
    return Path.home() / "mdv" / f"synth-anndata--{profile}--{_size_label(n_cells)}--g{n_genes}"


def _configure_cache_dirs() -> None:
    temp_root = Path(tempfile.gettempdir())
    mpl_dir = temp_root / "mdv-matplotlib"
    cache_dir = temp_root / "mdv-cache"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))


def _build_anndata(
    *,
    profile: str,
    n_cells: int,
    n_genes: int,
    seed: int,
    add_missing: bool,
    sparse_density: float,
) -> sc.AnnData:
    np.random.seed(seed)

    if profile == "minimal":
        return create_minimal_anndata(n_cells, n_genes, add_missing=add_missing)

    factory = MockAnnDataFactory(random_seed=seed)
    if profile == "realistic":
        return factory.create_realistic(n_cells, n_genes, add_missing=add_missing)
    if profile == "large":
        return factory.create_large(
            n_cells, n_genes, add_missing=add_missing, density=sparse_density
        )
    if profile == "memory-efficient":
        return factory.create_memory_efficient_large(
            n_cells, n_genes, add_missing=add_missing, density=sparse_density
        )
    raise SystemExit(f"unknown profile: {profile}")


def generate_project(
    *,
    output: Path,
    profile: str,
    n_cells: int,
    n_genes: int,
    seed: int,
    add_missing: bool,
    sparse_density: float,
    chunk_data: bool,
    compute_x_umap: bool,
    extra_expression_layers: bool,
    force: bool,
) -> None:
    if output.exists():
        if not force:
            raise SystemExit(f"{output} already exists. Pass --force to replace it.")
        shutil.rmtree(output)

    output.parent.mkdir(parents=True, exist_ok=True)
    _configure_cache_dirs()

    adata = _build_anndata(
        profile=profile,
        n_cells=n_cells,
        n_genes=n_genes,
        seed=seed,
        add_missing=add_missing,
        sparse_density=sparse_density,
    )
    if extra_expression_layers:
        add_synth_expression_layers(adata)
    mdv = convert_scanpy_to_mdv(
        str(output),
        adata,
        delete_existing=True,
        chunk_data=chunk_data,
        compute_x_umap=compute_x_umap,
    )

    state = mdv.state
    state["provenance"] = {
        "created_by": "mdvtools synthetic anndata generator",
        "generator": "mock_anndata",
        "profile": profile,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "seed": seed,
        "add_missing": add_missing,
        "sparse_density": sparse_density,
        "chunk_data": chunk_data,
        "compute_x_umap": compute_x_umap,
        "extra_expression_layers": extra_expression_layers,
        "cleanup_group": "synth-anndata",
    }
    mdv.state = state

    print(f"Created synthetic AnnData MDV project: {output}")
    print("Run /rescan_projects or restart the DB-backed app to register it.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an MDV project from synthetic AnnData (no SpatialData).",
    )
    parser.add_argument(
        "--profile",
        choices=PROFILE_CHOICES,
        default="minimal",
        help=(
            "minimal: small categorical/metadata patterns matching legacy zip mocks; "
            "realistic: reductions/layers/uns like typical SC tooling; "
            "large / memory-efficient: stress-oriented sparse matrices"
        ),
    )
    parser.add_argument(
        "--n-cells",
        type=int,
        default=1000,
        help="Number of observation rows.",
    )
    parser.add_argument(
        "--n-genes",
        type=int,
        default=50,
        help="Number of variable columns.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--add-missing",
        action="store_true",
        help=(
            "Inject missing values where supported (minimal: fixed positions in obs/var; "
            "realistic/large/memory-efficient: mock_anndata add_missing)."
        ),
    )
    parser.add_argument(
        "--sparse-density",
        type=float,
        default=0.1,
        help="Sparsity hint for large / memory-efficient profiles (passed through mock_anndata).",
    )
    parser.add_argument(
        "--chunk-data",
        action="store_true",
        help="Use chunked transpose during conversion (slower, lower peak RAM).",
    )
    parser.add_argument(
        "--compute-x-umap",
        action="store_true",
        help="Run neighbors/UMAP/Leiden from X inside convert_scanpy_to_mdv.",
    )
    parser.add_argument(
        "--extra-expression-layers",
        action="store_true",
        help=(
            "Add AnnData.layers synth_layer_a and synth_layer_b so the project has "
            "multiple rows_as_columns subgroups after conversion."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output MDV project directory. Defaults to ~/mdv/synth-anndata--...",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    for flag, value in (
        ("--n-cells", args.n_cells),
        ("--n-genes", args.n_genes),
        ("--sparse-density", args.sparse_density),
    ):
        if value <= 0:
            parser.error(f"{flag} must be greater than zero")
        if flag == "--sparse-density" and value > 1:
            parser.error("--sparse-density must be between 0 and 1")

    return args


def main() -> None:
    args = parse_args()
    output = args.output or _default_output(args.profile, args.n_cells, args.n_genes)
    generate_project(
        output=output.expanduser(),
        profile=args.profile,
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        seed=args.seed,
        add_missing=args.add_missing,
        sparse_density=args.sparse_density,
        chunk_data=args.chunk_data,
        compute_x_umap=args.compute_x_umap,
        extra_expression_layers=args.extra_expression_layers,
        force=args.force,
    )


if __name__ == "__main__":
    main()
