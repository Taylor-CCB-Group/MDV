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
