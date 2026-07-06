import h5py
import time
import numpy as np
import pandas as pd
from pathlib import Path

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


def test_get_column_provenance_resolves_record(tmp_path):
    project = _make_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch", max_concurrent=2)
    job_id = _submit_concat(mgr)
    _drive(mgr)

    prov = project.get_column_provenance("cells", "sample_cluster")
    assert prov is not None
    assert prov["_resolved"] is True
    assert prov["job_id"] == job_id
    assert prov["params"]["output_name"] == "sample_cluster"


def test_get_column_provenance_absent_returns_none(tmp_path):
    project = _make_project(tmp_path)
    # 'sample' came from add_datasource, never produced by a job -> no pointer
    assert project.get_column_provenance("cells", "sample") is None


def test_get_column_provenance_dangling_when_record_purged(tmp_path):
    project = _make_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch", max_concurrent=2)
    job_id = _submit_concat(mgr)
    _drive(mgr)

    # purge the durable record, leaving the column pointer dangling
    (Path(project.dir) / "jobs" / "records" / f"{job_id}.json").unlink()

    prov = project.get_column_provenance("cells", "sample_cluster")
    assert prov is not None
    assert prov["_resolved"] is False
    # denormalized still describes the column when the record is gone
    assert prov["job_id"] == job_id
    assert prov["tool_id"] == "concat_columns"
    assert len(prov["content_hash"]) == 16
    # ...but the full-record-only fields are not present
    assert "params" not in prov


def test_ingest_points_to_newest_job(tmp_path):
    project = _make_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch", max_concurrent=2)
    first = _submit_concat(mgr)
    _drive(mgr)
    second = _submit_concat(mgr)
    _drive(mgr)

    assert first != second
    # the column pointer is the newest run (last writer wins)
    pointer = project.get_column_metadata("cells", "sample_cluster")["provenance"]
    assert pointer["job_id"] == second
    second_prov = project.get_column_provenance("cells", "sample_cluster")
    assert second_prov is not None
    assert second_prov["job_id"] == second

    # history is preserved: both records still exist as separate files (ADR0009 -
    # records are append-only history; the column pointer is the current truth)
    records = Path(project.dir) / "jobs" / "records"
    assert (records / f"{first}.json").exists()
    assert (records / f"{second}.json").exists()
