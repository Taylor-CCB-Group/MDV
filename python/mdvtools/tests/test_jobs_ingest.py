import h5py
import numpy as np
import pandas as pd

from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.workspace import Workspace
from mdvtools.jobs.ingest import ingest_column_output


def _make_project(tmp_path):
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    df = pd.DataFrame({"sample": ["s1", "s2", "s3"], "cluster": ["a", "b", "a"]})
    project.add_datasource("cells", df)
    return project


def _write_result(ws, output_name, values):
    s = h5py.string_dtype()
    with h5py.File(ws.output / "result.h5", "w") as f:
        f.create_dataset(output_name, data=np.array(values, dtype=s))
        f.attrs["output_name"] = output_name


def test_ingest_add_column(tmp_path):
    project = _make_project(tmp_path)
    ws = Workspace(tmp_path / "jobs", "job1")
    params = {"datasource": "cells", "output_name": "donor_tissue"}

    _write_result(ws, "donor_tissue", ["s1_a", "s2_b", "s3_a"])
    ingest_column_output(project, params, ws)

    fields = [c["field"] for c in project.get_datasource_metadata("cells")["columns"]]
    assert fields.count("donor_tissue") == 1
    assert project.get_column("cells", "donor_tissue") == ["s1_a", "s2_b", "s3_a"]


def test_ingest_is_idempotent_replaces_not_duplicates(tmp_path):
    project = _make_project(tmp_path)
    ws = Workspace(tmp_path / "jobs", "job1")
    params = {"datasource": "cells", "output_name": "donor_tissue"}

    _write_result(ws, "donor_tissue", ["s1_a", "s2_b", "s3_a"])
    ingest_column_output(project, params, ws)

    # re-ingest: with same output, different values
    _write_result(ws, "donor_tissue", ["s1_b", "s2_a", "s3_b"])
    ingest_column_output(project, params, ws)

    fields = [c["field"] for c in project.get_datasource_metadata("cells")["columns"]]
    assert fields.count("donor_tissue") == 1  # replaced column not a second column
    assert project.get_column("cells", "donor_tissue") == ["s1_b", "s2_a", "s3_b"]


def test_ingest_lands_multiple_numeric_columns(tmp_path):
    project = _make_project(tmp_path)
    ws = Workspace(tmp_path / "jobs", "job1")
    params = {"datasource": "cells", "output_name": "UMAP"}

    u1 = np.array([0.1, 0.2, 0.3], dtype=np.float64)
    u2 = np.array([1.5, 2.5, 3.5], dtype=np.float64)
    with h5py.File(ws.output / "result.h5", "w") as f:
        f.create_dataset("UMAP_1", data=u1)
        f.create_dataset("UMAP_2", data=u2)
        f.attrs["output_name"] = "UMAP"          # rides as attr, not a column key

    result = ingest_column_output(project, params, ws)

    assert result["outputs"] == [("cells", "UMAP_1"), ("cells", "UMAP_2")]   # N outputs, not 1
    cols = {c["field"]: c["datatype"] for c in project.get_datasource_metadata("cells")["columns"]}
    assert cols["UMAP_1"] == "double" and cols["UMAP_2"] == "double"          # numeric, NOT text
    assert np.allclose(project.get_column("cells", "UMAP_1"), u1, atol=1e-4)
    assert np.allclose(project.get_column("cells", "UMAP_2"), u2, atol=1e-4)
