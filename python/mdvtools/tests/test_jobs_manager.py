import time
import pandas as pd

from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.manager import JobManager
from mdvtools.jobs.jobstore import Status


def _make_project(tmp_path):
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    df = pd.DataFrame({"sample": ["s1", "s2", "s3"], "cluster": ["a", "b", "a"]})
    project.add_datasource("cells", df)
    return project


def _drive(mgr, timeout=60):
    deadline = time.time() + timeout
    while time.time() < deadline:
        mgr.tick()
        statuses = [r.status for r in mgr.store.load_all()]
        if all(s in (Status.DONE.value, Status.FAILED.value) for s in statuses):
            return
        time.sleep(0.1)
    raise AttributeError("jobs did not finish in time")


def test_end_to_end_single_job(tmp_path):
    project = _make_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch", max_concurrent=2)
    mgr.submit(
        "concat_columns",
        {
            "datasource": "cells",
            "column_a": "sample",
            "column_b": "cluster",
            "separator": "_",
            "output_name": "sample_cluster",
        },
    )
    _drive(mgr)

    assert [r.status for r in mgr.store.load_all()] == [Status.DONE.value]
    assert project.get_column("cells", "sample_cluster") == ["s1_a", "s2_b", "s3_a"]


def test_max_concurrent_holds_extra_jobs_queued(tmp_path):
    project = _make_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch", max_concurrent=2)

    for i in range(3):
        mgr.submit(
            "concat_columns",
            {
                "datasource": "cells",
                "column_a": "sample",
                "column_b": "cluster",
                "separator": "_",
                "output_name": f"out_{i}",
            },
        )

    # right after submitting 3 with a bound of 2, before any tick: 2 active, 1 queued
    statuses = sorted(r.status for r in mgr.store.load_all())
    assert statuses.count(Status.RUNNING.value) == 2
    assert statuses.count(Status.QUEUED.value) == 1

    _drive(mgr)

    # all three columns landed
    for i in range(3):
        assert project.get_column("cells", f"out_{i}") == ["s1_a", "s2_b", "s3_a"]


def test_provenance_promoted_and_workspace_cleaned(tmp_path):
    project = _make_project(tmp_path)
    workspace_root = tmp_path / "scratch"
    mgr = JobManager(project, workspace_root=workspace_root, max_concurrent=2)
    job_id = mgr.submit(
        "concat_columns",
        {
            "datasource": "cells",
            "column_a": "sample",
            "column_b": "cluster",
            "separator": "_",
            "output_name": "sample_cluster",
        },
    )
    _drive(mgr)

    rec = {r.job_id: r for r in mgr.store.load_all()}[job_id]
    assert rec.status == Status.DONE.value
    prov = rec.provenance
    assert prov is not None
    assert prov["tool_id"] == "concat_columns"
    assert prov["params"]["output_name"] == "sample_cluster"
    assert prov["output"] == {"rows": 3}
    assert len(prov["content_hash"]) == 16
    # workspace scratch (keyed by job_id, outside the project) is GC'd on success
    assert not (workspace_root / job_id).exists()
