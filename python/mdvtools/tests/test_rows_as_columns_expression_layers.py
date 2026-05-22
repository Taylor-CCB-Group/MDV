"""
Extra AnnData.layers become multiple rows_as_columns subgroups (same path as
SpatialData merge → convert_scanpy_to_mdv). Synthetic SpatialData can add
layers via generate_synthetic_spatial_project --extra-expression-layers.
"""

import numpy as np
import scanpy as sc

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject


def test_annadata_layers_yield_multiple_subgroups(tmp_path):
    n_cells, n_genes = 30, 10
    rng = np.random.default_rng(0)
    x = rng.random((n_cells, n_genes)).astype(np.float32)
    adata = sc.AnnData(X=x)
    adata.layers["synth_layer_a"] = np.asarray(x, copy=True)
    adata.layers["synth_layer_b"] = np.log1p(x + 1e-6)

    out = tmp_path / "proj"
    convert_scanpy_to_mdv(str(out), adata, delete_existing=True)

    project = MDVProject(str(out))
    links = project.get_links("cells", "rows_as_columns")
    assert links
    rac = links[0]["link"]["rows_as_columns"]
    names = set((rac.get("subgroups") or {}).keys())
    assert "gs" in names
    assert "synth_layer_a" in names
    assert "synth_layer_b" in names
    assert len(names) >= 3
