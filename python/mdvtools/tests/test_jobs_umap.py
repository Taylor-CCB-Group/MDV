from typing import cast
import json

import h5py
import numpy as np
import pandas as pd
import scipy.sparse

from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.workspace import Workspace, materialize_matrix_tray
from mdvtools.jobs.registry import ToolSpec, ParamSpec, OutputSpec
from mdvtools.jobs.workers.umap_worker import run as umap_run

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

def _write_matrix_tray(ws, X, out_name="UMAP", neighbors_kwargs=None, umap_kwargs=None):
    """Hand-build the tray materialize_matrix_tray would produce, so this test crosses the courier
    boundary the same way the worker will in production - and stays MDV-free, like the worker."""
    X = scipy.sparse.csc_matrix(X)
    with h5py.File(ws.input / "tray.h5", "w") as f:
        f.create_dataset("x", data=X.data.astype(np.float32))
        f.create_dataset("i", data=X.indices.astype(np.uint32))
        f.create_dataset("p", data=X.indptr)
        n_cells, n_genes = cast(tuple[int, int], X.shape)
        f.attrs["n_cells"] = n_cells
        f.attrs["n_genes"] = n_genes
        f.attrs["sparse"] = True
        f.attrs["output_name"] = out_name
        if neighbors_kwargs is not None:
            f.attrs["kwargs.neighbors"] = json.dumps(neighbors_kwargs)
        if umap_kwargs is not None:
            f.attrs["kwargs.umap"] = json.dumps(umap_kwargs)


def test_umap_worker_embeds_and_survives_nonfinite(tmp_path):
    rng = np.random.default_rng(0)
    dense = rng.random((60, 10)).astype(np.float32)
    dense[dense < 0.6] = 0.0        # ~40% nonzero -> genuinely sparse
    dense[3, 4] = np.inf            # non-finite entries the worker must sanitize before neighbors
    dense[7, 1] = np.nan

    ws = Workspace(tmp_path / "scratch", "job1")
    _write_matrix_tray(ws, dense, out_name="UMAP")

    umap_run(str(ws.root))

    assert ws.read_marker() == "done"
    with h5py.File(ws.output / "result.h5", "r") as f:
        u1 = cast(h5py.Dataset, f["UMAP_1"])[:]
        u2 = cast(h5py.Dataset, f["UMAP_2"])[:]
    assert u1.shape == (60,) and u2.shape == (60,)
    assert np.isfinite(u1).all() and np.isfinite(u2).all()   # proves the sanitize worked

    manifest = json.loads((ws.output / "manifest.json").read_text())
    assert manifest["rows"] == 60

def test_matrix_tray_bundles_kwargs_by_applies_to(tmp_path):
    project, X = _make_matrix_project(tmp_path)
    spec = ToolSpec(
        id="umap", name="UMAP", description="test",
        params=[
            ParamSpec("datasource", "dropdown", "Datasource"),
            ParamSpec("layer", "dropdown", "Matrix", default="gs"),
            ParamSpec("output_name", "text", "Base name", default="UMAP"),
            ParamSpec("n_neighbors", "int", "Neighbors", applies_to="neighbors"),
            ParamSpec("min_dist", "float", "Min dist", applies_to="umap"),
            ParamSpec("n_components", "int", "Dims", applies_to="umap"),
        ],
        output=OutputSpec("column", "datasource", "output_name"),
        entrypoint="mdvtools.jobs.workers.umap_worker:run",
        input_shape="matrix",
    )
    ws = Workspace(tmp_path / "scratch", "job1")
    params = {"datasource": "cells", "layer": "gs", "output_name": "UMAP",
              "n_neighbors": 8, "min_dist": 0.2, "n_components": 3}

    materialize_matrix_tray(project, spec, params, ws)

    with h5py.File(ws.input / "tray.h5", "r") as f:
        assert f.attrs["output_name"] == "UMAP"
        assert json.loads(cast(str, f.attrs.get("kwargs.neighbors", "{}"))) == {"n_neighbors": 8}
        assert json.loads(cast(str, f.attrs.get("kwargs.umap", "{}"))) == {"min_dist": 0.2, "n_components": 3}

def test_umap_worker_honors_n_components(tmp_path):
    rng = np.random.default_rng(0)
    dense = rng.random((60, 10)).astype(np.float32)
    dense[dense < 0.6] = 0.0

    ws = Workspace(tmp_path / "scratch", "job1")
    _write_matrix_tray(ws, dense, out_name="UMAP", umap_kwargs={"n_components": 3})

    umap_run(str(ws.root))

    assert ws.read_marker() == "done"
    with h5py.File(ws.output / "result.h5", "r") as f:
        assert "UMAP_1" in f and "UMAP_2" in f and "UMAP_3" in f   # 3 dims -> 3 columns
    manifest = json.loads((ws.output / "manifest.json").read_text())
    assert manifest["columns"] == 3
