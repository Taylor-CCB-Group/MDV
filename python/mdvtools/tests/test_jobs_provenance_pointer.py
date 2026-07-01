import h5py
import numpy as np
import pandas as pd

from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.workspace import Workspace
from mdvtools.jobs.ingest import ingest_column_output


# define a project
def _make_project(tmp_path):
    project = MDVProject(tmp_path)
    df = pd.DataFrame({"sample": ["s1", "s2", "s3"], "cluster": ["a", "b", "a"]})
    project.add_datasource("cells", df)
    return project


def _write_results(ws, output_name, values):
    s = h5py.string_dtype()
    with h5py.File(ws.output / "result.h5", "w") as f:
        f.create_dataset(output_name, data=np.array(values, dtype=s))
        f.attrs["output_name"] = output_name


def test_ingester_reports_outputs(tmp_path):
    project = _make_project(tmp_path)
    ws = Workspace(tmp_path / "jobs", "job1")
    params = {"datasource": "cells", "output_name": "donor_tissue"}
    _write_results(ws, "donor_tissue", ["s1_a", "s2_b", "s3_a"])

    result = ingest_column_output(project, params, ws)
    assert result is not None  # to avoid basepyright error
    assert result["outputs"] == [("cells", "donor_tissue")]
    assert "manifest" in result  # present even when no manifest.json was written (None)
