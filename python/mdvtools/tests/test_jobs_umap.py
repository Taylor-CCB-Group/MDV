from typing import cast

import h5py
import numpy as np
import pandas as pd
import scipy.sparse

from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.workspace import Workspace, materialize_matrix_tray
from mdvtools.jobs.registry import ToolSpec, ParamSpec, OutputSpec

# minimal spec
UMAP_SPEC = ToolSpec(
    id = "umap",
    name = "UMAP",
    description = "test",
    params =
    [
        ParamSpec("datasource", "dropdown", "Datasource"),
        ParamSpec("layer", "dropdown", "Matrix", default="gs"),
        ParamSpec("output_name", "text", "New column base name", default="UMAP"),
    ],
    output=OutputSpec("column", "datasource", "output_name"),
    entrypoint="mdvtools.jobs.workers.umap_worker:run",
    input_shape="matrix",
)

def _make_matrix_project(tmp_path):
    project = MDVProject(tmp_path)
    project.add_datasource("cells", pd.DataFrame({"cell_id": ["c0", "c1", "c2"]}))
    project.add_datasource("genes", pd.DataFrame({"name": ["g0", "g1"]}))
    project.add_rows_as_columns_link("cells", "genes", "name", "Gene Expr")
    X = scipy.sparse.csc_matrix(np.array([[1.0, 0.0], [0.0, 2.0], [3.0, 4.0]], dtype=np.float32))  # 3 cells x 2 genes
    project.add_rows_as_columns_subgroup("cells", "genes", "gs", X, name="gene_scores")
    return project, X

def test_matrix_tray_roundtrips_stored_X(tmp_path):
    project, X = _make_matrix_project(tmp_path)
    ws = Workspace(tmp_path / "scratch", "job1")
    params = {"datasource": "cells", "layer": "gs", "output_name": "UMAP"}

    materialize_matrix_tray(project, UMAP_SPEC, params, ws)

    with h5py.File(ws.input / "tray.h5", "r") as f:
        n_cells, n_genes = int(cast(int, f.attrs["n_cells"])), int(cast(int, f.attrs["n_genes"]))
        x = cast(h5py.Dataset, f["x"])
        i = cast(h5py.Dataset, f["i"])
        p = cast(h5py.Dataset, f["p"])
        rebuilt = scipy.sparse.csc_matrix(
            (x[:], i[:], p[:]), shape=(n_cells, n_genes)
        )
        assert (n_cells, n_genes) == (3, 2)
        assert np.allclose(rebuilt.toarray(), X.toarray())
