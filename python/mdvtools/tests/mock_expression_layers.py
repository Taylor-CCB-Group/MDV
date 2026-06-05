"""Optional extra AnnData.layers for MDV rows_as_columns subgroup testing."""

from __future__ import annotations

import numpy as np
import scanpy as sc


def add_synth_expression_layers(adata: sc.AnnData) -> None:
    """Add two layers alongside X so convert_scanpy_to_mdv creates multiple link subgroups."""
    x = np.asarray(adata.X, dtype=np.float32)
    adata.layers["synth_layer_a"] = np.asarray(x, copy=True)
    adata.layers["synth_layer_b"] = np.log1p(x + 1e-6).astype(np.float32)
