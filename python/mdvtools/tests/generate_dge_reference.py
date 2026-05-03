"""Generate Scanpy DGE reference data for mdv-dge validation.

Runs rank_genes_groups on the PBMC3k dataset for all leiden clusters vs rest,
exports results as JSON fixtures for use in vitest browser-side DGE tests.

Usage:
    cd /app/python && poetry run python -m mdvtools.tests.generate_dge_reference
"""

import json
import sys
import numpy as np
import scanpy as sc


H5AD_PATH = "/app/mdv/16/scanpy-pbmc3k.h5ad"
FIXTURE_DIR = "/app/src/tests/fixtures"
TOP_N_SMALL = 50
VALIDATION_GENES = 100


def extract_results(adata, method: str) -> dict:
    """Extract rank_genes_groups results into a serializable dict."""
    result = sc.get.rank_genes_groups_df(adata, group=None)

    groups = {}
    for group_name in result["group"].unique():
        gdf = result[result["group"] == group_name].copy()
        gdf = gdf.sort_values("pvals")
        groups[str(group_name)] = {
            "genes": gdf["names"].tolist(),
            "scores": [float(x) for x in gdf["scores"]],
            "log2fc": [float(x) for x in gdf["logfoldchanges"]],
            "pvals": [float(x) for x in gdf["pvals"]],
            "pvals_adj": [float(x) for x in gdf["pvals_adj"]],
        }

    group_sizes = adata.obs["leiden"].value_counts().to_dict()

    return {
        "method": method,
        "groupby": "leiden",
        "reference": "rest",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "group_sizes": {str(k): int(v) for k, v in group_sizes.items()},
        "groups": groups,
    }


def extract_small(full_results: dict) -> dict:
    """Extract top N genes per group for fast unit tests."""
    small = {k: v for k, v in full_results.items() if k != "groups"}
    small["top_n"] = TOP_N_SMALL
    small["groups"] = {}
    for gname, gdata in full_results["groups"].items():
        small["groups"][gname] = {
            k: v[:TOP_N_SMALL] for k, v in gdata.items()
        }
    return small


def extract_validation_data(adata) -> dict:
    """Extract log1p-normalized expression data for a subset of genes + leiden assignments.

    Uses adata.raw.X (log1p-normalized) since that's what rank_genes_groups uses
    by default. This allows the JS test to recompute DGE from scratch and compare.
    """
    # Use raw data (log1p-normalized) to match rank_genes_groups default behavior
    raw = adata.raw
    raw_gene_names = raw.var_names.tolist()

    # Select genes evenly spaced from the raw gene set, but only those that
    # appear in the rank_genes_groups results (i.e. all raw genes)
    selected_indices = np.linspace(0, len(raw_gene_names) - 1, VALIDATION_GENES, dtype=int)
    selected_genes = [raw_gene_names[i] for i in selected_indices]

    X = raw.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    expression = {}
    for idx, gene in zip(selected_indices, selected_genes):
        expression[gene] = [float(x) for x in X[:, idx]]

    leiden = adata.obs["leiden"].tolist()

    return {
        "n_cells": int(adata.n_obs),
        "n_genes_total": int(raw.n_vars),
        "n_genes_selected": len(selected_genes),
        "genes": selected_genes,
        "leiden": leiden,
        "expression": expression,
    }


def write_json(data: dict, path: str):
    with open(path, "w") as f:
        json.dump(data, f, allow_nan=False, separators=(",", ":"))
    size_mb = len(json.dumps(data, separators=(",", ":"))) / 1024 / 1024
    print(f"  Written {path} ({size_mb:.1f} MB)")


def main():
    print(f"Loading {H5AD_PATH}...")
    adata = sc.read_h5ad(H5AD_PATH)
    print(f"  {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"  Leiden clusters: {sorted(adata.obs['leiden'].unique())}")

    # --- T-test ---
    print("\nRunning rank_genes_groups (t-test)...")
    sc.tl.rank_genes_groups(adata, "leiden", method="t-test", reference="rest")
    ttest_full = extract_results(adata, "t-test")
    ttest_small = extract_small(ttest_full)

    write_json(ttest_full, f"{FIXTURE_DIR}/dge_reference_ttest.json")
    write_json(ttest_small, f"{FIXTURE_DIR}/dge_reference_ttest_small.json")

    # --- Wilcoxon ---
    print("\nRunning rank_genes_groups (wilcoxon)...")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", reference="rest")
    wilcoxon_full = extract_results(adata, "wilcoxon")
    wilcoxon_small = extract_small(wilcoxon_full)

    write_json(wilcoxon_full, f"{FIXTURE_DIR}/dge_reference_wilcoxon.json")
    write_json(wilcoxon_small, f"{FIXTURE_DIR}/dge_reference_wilcoxon_small.json")

    # --- Validation expression data ---
    print("\nExtracting validation expression data...")
    validation = extract_validation_data(adata)
    write_json(validation, f"{FIXTURE_DIR}/dge_validation_data.json")

    print("\nDone. Generated fixtures:")
    print(f"  {FIXTURE_DIR}/dge_reference_ttest.json")
    print(f"  {FIXTURE_DIR}/dge_reference_ttest_small.json")
    print(f"  {FIXTURE_DIR}/dge_reference_wilcoxon.json")
    print(f"  {FIXTURE_DIR}/dge_reference_wilcoxon_small.json")
    print(f"  {FIXTURE_DIR}/dge_validation_data.json")


if __name__ == "__main__":
    main()
