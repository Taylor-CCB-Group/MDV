#!/usr/bin/env python3
"""
Create the PBMC3K MDV project used by ``prompts_pbmc3k.csv`` ChatMDV CSV benchmarks.

Default output: ``/app/mdv/pbmc3k_chat`` (override with ``--output`` or env
``MDV_PBMC3K_CHAT_PROJECT``).

The project is built from Scanpy's ``pbmc3k_processed`` dataset with Leiden
clusters computed on the existing neighbor graph when ``leiden`` is not already
present. Conversion uses ``convert_scanpy_to_mdv``, which exports:

- ``cells``: UMAP/PCA embeddings, ``leiden``, QC fields (``n_genes``, ``n_counts``,
  ``percent_mito``), and a rows-as-columns gene expression link.
- ``genes``: gene metadata plus expression matrix via subgroup ``gs``.

Example (from repository ``python/`` directory):

  python mdvtools/test_projects/chatmdv_csv_eval/setup_pbmc3k_chat.py

  python mdvtools/test_projects/chatmdv_csv_eval/setup_pbmc3k_chat.py --force

Then run benchmarks:

  python mdvtools/test_projects/chatmdv_csv_eval/run_prompt_csv.py \\
    --csv mdvtools/test_projects/chatmdv_csv_eval/prompts_pbmc3k.csv --limit 1
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import scanpy as sc

_PY_ROOT = Path(__file__).resolve().parents[3]
if str(_PY_ROOT) not in sys.path:
    sys.path.insert(0, str(_PY_ROOT))

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject

DEFAULT_OUTPUT = Path(os.environ.get("MDV_PBMC3K_CHAT_PROJECT", "/app/mdv/pbmc3k_chat"))


def _ensure_leiden(adata: sc.AnnData, *, resolution: float) -> None:
    if "leiden" in adata.obs.columns:
        adata.obs["leiden"] = pd.Categorical(adata.obs["leiden"].astype(str))
        return
    if "neighbors" not in adata.uns:
        raise RuntimeError(
            "AnnData has no 'leiden' column and no 'neighbors' graph in .uns; "
            "cannot compute Leiden clusters for the chat benchmark project."
        )
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    adata.obs["leiden"] = pd.Categorical(adata.obs["leiden"].astype(str))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"MDV project directory (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing project directory.",
    )
    parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=0.9,
        help="Leiden resolution when computing clusters (default: 0.9).",
    )
    return parser.parse_args()


def setup_pbmc3k_chat(
    output_path: Path,
    *,
    force: bool = False,
    leiden_resolution: float = 0.9,
) -> MDVProject:
    output_path = output_path.expanduser().resolve()
    if output_path.exists() and not force:
        raise FileExistsError(
            f"{output_path} already exists. Pass --force to recreate the project."
        )

    print("Loading scanpy dataset pbmc3k_processed...")
    adata = sc.datasets.pbmc3k_processed()
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    _ensure_leiden(adata, resolution=leiden_resolution)
    print(f"  leiden clusters: {adata.obs['leiden'].nunique()}")

    print(f"Creating MDV project at {output_path}...")
    convert_scanpy_to_mdv(str(output_path), adata, delete_existing=force or not output_path.exists())

    project = MDVProject(str(output_path))
    state = project.state
    state["provenance"] = {
        "created_by": "setup_pbmc3k_chat.py",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": "scanpy",
        "parameters": {
            "dataset": "pbmc3k_processed",
            "leiden_resolution": leiden_resolution,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        },
    }
    project.state = state
    project.set_editable(True)

    names = project.get_datasource_names()
    required_obs = ("leiden", "X_umap_1", "X_umap_2", "X_pca_1", "n_genes", "percent_mito")
    cells_meta = project.get_datasource_metadata("cells")
    fields = {col["field"] for col in cells_meta["columns"]}
    missing = [c for c in required_obs if c not in fields]
    if missing:
        raise RuntimeError(f"Project cells datasource missing expected columns: {missing}")
    if "genes" not in names:
        raise RuntimeError("Project is missing 'genes' datasource.")

    print(f"Created MDV project at {output_path}")
    print(f"  datasources: {', '.join(names)}")
    return project


def main() -> int:
    args = _parse_args()
    try:
        setup_pbmc3k_chat(
            args.output,
            force=args.force,
            leiden_resolution=args.leiden_resolution,
        )
    except FileExistsError as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
