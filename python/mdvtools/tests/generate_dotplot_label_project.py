#!/usr/bin/env python3
"""
Generate a synthetic MDV project that reproduces dot plot x-axis label issues.

Builds AnnData with gene names of varied length, converts it with
``convert_scanpy_to_mdv`` and saves a ``dotplot labels`` view with one dot plot
per scenario:

- few_short: 4 short genes in a wide chart; labels fit flat but are rotated
- long_clipped: long gene names; labels cut off by the fixed bottom margin
- crowded: 60 genes in a narrow chart; labels overlap
- long_first: long first label; spills left under the y-axis labels
- mixed_subgroups: fields from two subgroups; the suffix must stay

Gene expression is added as an ``rna_expr`` layer so x labels read like
``CD4(rna_expr)``; X becomes the ``gs`` subgroup used by ``mixed_subgroups``.

Charts carry no ``axis`` config so the DotPlot constructor applies its defaults.

Output defaults to a direct child of ``~/mdv`` so ``GET /rescan_projects`` can
register it. Inside the dev container, pass ``--output /app/mdv/<name>``.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import scanpy as sc

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.llm.column_field_resolve import build_expression_wrapper_token
from mdvtools.tests.generate_synthetic_anndata_project import _configure_cache_dirs
from mdvtools.tests.mock_anndata import MockAnnDataFactory

DEFAULT_NAME = "synth-dotplot-labels"
VIEW_NAME = "dotplot labels"
CATEGORY_COLUMN = "cell_type"
EXPR_SUBGROUP = "rna_expr"
X_SUBGROUP = "gs"

SHORT_GENES = ["CD4", "CD8A", "GZMB", "NKG7", "MS4A1", "LYZ", "CD14", "IL7R"]
LONG_GENES = [
    "ENSG00000228253-AS1",
    "LINC01409_novel_transcript",
    "HLA-DRB1_alt_haplotype",
    "RP11-206L10.9",
    "MTRNR2L12_pseudogene",
    "ENSG00000254876.2",
    "SNHG29_antisense_RNA",
    "LOC105369174_uncharacterized",
    "AC011462.4_lincRNA",
    "TRBV20-1_variable_segment",
]
VERY_LONG_GENE = "ENSG00000283907_readthrough_long_noncoding"
N_CROWDED = 60


def _gene_names() -> list[str]:
    crowded = [f"GENE{i:02d}-{'x' * (i % 9)}" for i in range(N_CROWDED)]
    return SHORT_GENES + LONG_GENES + [VERY_LONG_GENE] + crowded


def _build_anndata(n_cells: int, seed: int) -> sc.AnnData:
    np.random.seed(seed)
    genes = _gene_names()
    factory = MockAnnDataFactory(random_seed=seed)
    adata = factory.create_with_specific_features(n_cells=n_cells, n_genes=len(genes))
    adata.var_names = genes
    adata.layers[EXPR_SUBGROUP] = np.asarray(adata.X, dtype=np.float32).copy()
    return adata


def _scenarios(gene_index: dict[str, int]) -> list[dict[str, Any]]:
    def fields(genes: list[str], subgroup: str = EXPR_SUBGROUP) -> list[str]:
        return [build_expression_wrapper_token(subgroup, g, gene_index[g]) for g in genes]

    crowded = [g for g in gene_index if g.startswith("GENE")]
    return [
        {
            "id": "few_short",
            "title": "few_short: 4 short labels (fit flat, still rotated)",
            "fields": fields(SHORT_GENES[:4]),
            "size": [700, 400],
            "position": [10, 10],
        },
        {
            "id": "long_clipped",
            "title": "long_clipped: long labels (cut off at bottom)",
            "fields": fields(LONG_GENES),
            "size": [700, 400],
            "position": [720, 10],
        },
        {
            "id": "crowded",
            "title": "crowded: 60 labels in narrow chart (overlap)",
            "fields": fields(crowded),
            "size": [450, 400],
            "position": [10, 420],
        },
        {
            "id": "long_first",
            "title": "long_first: long first label (spills under y-axis)",
            "fields": fields([VERY_LONG_GENE] + SHORT_GENES[4:8]),
            "size": [500, 400],
            "position": [470, 420],
        },
        {
            "id": "mixed_subgroups",
            "title": "mixed_subgroups: two subgroups (keep suffix)",
            "fields": fields(SHORT_GENES[:3]) + fields(SHORT_GENES[:3], X_SUBGROUP),
            "size": [700, 400],
            "position": [10, 830],
        },
    ]


def _dot_plot(scenario: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": scenario["id"],
        "title": scenario["title"],
        "legend": "",
        "type": "dot_plot",
        "param": [CATEGORY_COLUMN, *scenario["fields"]],
        "size": scenario["size"],
        "position": scenario["position"],
    }


def generate_project(*, output: Path, n_cells: int, seed: int, force: bool) -> None:
    if output.exists():
        if not force:
            raise SystemExit(f"{output} already exists. Pass --force to replace it.")
        shutil.rmtree(output)

    output.parent.mkdir(parents=True, exist_ok=True)
    _configure_cache_dirs()

    adata = _build_anndata(n_cells, seed)
    gene_index = {g: i for i, g in enumerate(adata.var_names)}
    mdv = convert_scanpy_to_mdv(str(output), adata, delete_existing=True)

    charts = [_dot_plot(s) for s in _scenarios(gene_index)]
    mdv.set_view(VIEW_NAME, {"initialCharts": {"cells": charts, "genes": []}}, True)

    state = mdv.state
    state["provenance"] = {
        "created_by": "mdvtools dotplot label generator",
        "n_cells": n_cells,
        "seed": seed,
        "scenarios": [c["id"] for c in charts],
        "cleanup_group": "synth-dotplot-labels",
    }
    mdv.state = state

    print(f"Created dot plot label MDV project: {output}")
    print(f"Open view '{VIEW_NAME}'. Run /rescan_projects to register it.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an MDV project whose dot plots reproduce x-axis label issues.",
    )
    parser.add_argument("--n-cells", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=f"Output MDV project directory. Defaults to ~/mdv/{DEFAULT_NAME}",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.n_cells <= 0:
        parser.error("--n-cells must be greater than zero")
    return args


def main() -> None:
    args = parse_args()
    output = args.output or Path.home() / "mdv" / DEFAULT_NAME
    generate_project(
        output=output.expanduser(), n_cells=args.n_cells, seed=args.seed, force=args.force
    )


if __name__ == "__main__":
    main()
