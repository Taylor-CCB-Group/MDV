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


def test_reverse_link_stores_cells_as_gene_columns(tmp_path):
    import h5py
    import scipy.sparse

    n_cells, n_genes = 12, 5
    rng = np.random.default_rng(1)
    x = rng.random((n_cells, n_genes)).astype(np.float32)
    x[x < 0.5] = 0
    adata = sc.AnnData(X=scipy.sparse.csr_matrix(x))
    # a gene metadata column with the same name as the transformation label
    adata.var["raw"] = [f"r{i}" for i in range(n_genes)]

    out = tmp_path / "proj"
    convert_scanpy_to_mdv(
        str(out),
        adata,
        delete_existing=True,
        add_layer_data=False,
        add_reverse_link=True,
        x_data_transformation="raw",
    )

    project = MDVProject(str(out))
    links = project.get_links("genes", "rows_as_columns")
    assert len(links) == 1
    assert links[0]["datasource"] == "cells"
    rac = links[0]["link"]["rows_as_columns"]
    assert rac["name_column"] == "mdv_cell_id"
    assert set(rac["subgroups"]) == {"raw"}
    subgroup = rac["subgroups"]["raw"]
    assert subgroup["label"] == "raw"
    assert subgroup["name"] != "raw"
    assert list(project.get_column("genes", "raw")) == [f"r{i}" for i in range(n_genes)]

    # each cell's column on the genes datasource holds that cell's row of X
    with h5py.File(project.h5file, "r") as h5:
        gr = h5[f"genes/{subgroup['name']}"]
        assert isinstance(gr, h5py.Group)
        x_ds, i_ds, p_ds = gr["x"], gr["i"], gr["p"]
        assert isinstance(x_ds, h5py.Dataset)
        assert isinstance(i_ds, h5py.Dataset)
        assert isinstance(p_ds, h5py.Dataset)
        stored = scipy.sparse.csc_matrix(
            (x_ds[:], i_ds[:], p_ds[:]), shape=(n_genes, n_cells)
        )
    np.testing.assert_allclose(stored.toarray(), x.T)


def test_chunk_budget_counts_source_dtype_and_flattened_copy(tmp_path, monkeypatch):
    import psutil
    from types import SimpleNamespace

    available = 1024**3
    monkeypatch.setattr(psutil, "virtual_memory", lambda: SimpleNamespace(available=available))
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)

    rows = 65_536
    chunk = project._get_optimal_chunk_size(rows, 100_000, np.dtype(np.float64))
    # float64 slice + its flattened copy + the float32 copy h5py writes
    peak = chunk * rows * (8 + 8 + 4)
    assert peak <= available // 2
