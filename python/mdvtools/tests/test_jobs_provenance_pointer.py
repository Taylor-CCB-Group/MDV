import h5py
import time
import numpy as np
import pandas as pd

from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.workspace import Workspace
from mdvtools.jobs.ingest import ingest_column_output
from mdvtools.jobs.manager import JobManager
from mdvtools.jobs.jobstore import Status


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


def _drive(mgr, timeout=60):
    deadline = time.time() + timeout
    while time.time() < deadline:
        mgr.tick()
        statuses = [r.status for r in mgr.store.load_all()]
        if all(s in (Status.DONE.value, Status.FAILED.value) for s in statuses):
            return
        time.sleep(0.1)
    raise AttributeError("jobs did not finish in time.")


def _submit_concat(mgr, output_name="sample_cluster"):
    return mgr.submit(
        "concat_columns",
        {
            "datasource": "cells",
            "column_a": "sample",
            "column_b": "cluster",
            "separator": "_",
            "output_name": output_name,
        },
    )


def test_pointer_stamped_on_column(tmp_path):
    project = _make_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch", max_concurrent=2)
    job_id = _submit_concat(mgr)
    _drive(mgr)

    pointer = project.get_column_metadata("cells", "sample_cluster")["provenance"]
    assert pointer["kind"] == "job"
    assert pointer["job_id"] == job_id
    assert pointer["tool_id"] == "concat_columns"
    assert len(pointer["content_hash"]) == 16


def test_pointer_hash_matches_record(tmp_path):
    # the same content_hash on the column and the record - one hash function (ADR0009)
    project = _make_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch", max_concurrent=2)
    job_id = _submit_concat(mgr)
    _drive(mgr)

    pointer = project.get_column_metadata("cells", "sample_cluster")["provenance"]
    rec = {r.job_id: r for r in mgr.store.load_all()}[job_id]
    assert pointer["content_hash"] == rec.provenance["content_hash"]


def test_ingester_reports_outputs(tmp_path):
    project = _make_project(tmp_path)
    ws = Workspace(tmp_path / "jobs", "job1")
    params = {"datasource": "cells", "output_name": "donor_tissue"}
    _write_results(ws, "donor_tissue", ["s1_a", "s2_b", "s3_a"])

    result = ingest_column_output(project, params, ws)
    assert result is not None  # to avoid basepyright error
    assert result["outputs"] == [("cells", "donor_tissue")]
    assert "manifest" in result  # present even when no manifest.json was written (None)
