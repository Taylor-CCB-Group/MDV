import time
from pathlib import Path

import pandas as pd
import numpy as np
import scipy.sparse
from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.manager import JobManager
from mdvtools.jobs.jobstore import Status
from mdvtools.jobs.executor import Handle

class _FakeExecutor:
    """Backend stand-in: submit() records the call but runs NO worker (no marker is ever
    written); poll() returns a state the test controls. Lets us drive the manager's
    marker-absent → poll fallback deterministically, without caring about job output."""

    def __init__(self, poll_result="running"):
        self._poll = poll_result
        self.submits = 0

    def submit(self, entrypoint, workspace):
        self.submits += 1
        return Handle("fake", str(self.submits))

    def poll(self, handle):
        return self._poll

    def locate_result(self, handle, workspace):
        return Path(workspace) / "output"


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

def _make_matrix_project(tmp_path, n_cells=60, n_genes=10, seed=0):
    """A stored expression matrix (gs subgroup); ~60 cells so scanpy's default n_neighbors=15 fits."""
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    project.add_datasource("cells", pd.DataFrame({"cell_id": [f"c{i}" for i in range(n_cells)]}))
    project.add_datasource("genes", pd.DataFrame({"name": [f"g{j}" for j in range(n_genes)]}))
    project.add_rows_as_columns_link("cells", "genes", "name", "Gene Expr")
    rng = np.random.default_rng(seed)
    dense = rng.random((n_cells, n_genes)).astype(np.float32)
    dense[dense < 0.6] = 0.0                                  # genuinely sparse
    X = scipy.sparse.csc_matrix(dense)
    project.add_rows_as_columns_subgroup("cells", "genes", "gs", X, name="gene_scores")
    return project

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

def test_umap_job_end_to_end_lands_two_columns_with_provenance(tmp_path):
    project = _make_matrix_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch", max_concurrent=2)
    job_id = mgr.submit("umap", {"datasource": "cells", "layer": "gs", "output_name": "UMAP"})

    _drive(mgr, timeout=180)   # real scanpy UMAP in a subprocess - allow first-run import/JIT

    rec = {r.job_id: r for r in mgr.store.load_all()}[job_id]
    assert rec.status == Status.DONE.value

    cols = {c["field"]: c["datatype"] for c in project.get_datasource_metadata("cells")["columns"]}
    assert cols.get("UMAP_1") == "double" and cols.get("UMAP_2") == "double"   # numeric embedding
    assert len(project.get_column("cells", "UMAP_1")) == 60                    # one coord per cell

    for col in ("UMAP_1", "UMAP_2"):                                          # provenance on each
        prov = project.get_column_provenance("cells", col)
        assert prov is not None and prov["job_id"] == job_id and prov["tool_id"] == "umap"


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
    # temp fix
    assert prov["output"] == {"rows": 3}
    assert len(prov["content_hash"]) == 16

    # workspace scratch (keyed by job_id, outside the project) is GC'd on success
    assert not (workspace_root / job_id).exists()

def test_umap_job_end_to_end_honors_params(tmp_path):
    project = _make_matrix_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch", max_concurrent=2)
    job_id = mgr.submit(
        "umap",
        {"datasource": "cells", "layer": "gs", "output_name": "UMAP",
         "n_neighbors": 10, "n_components": 3},
    )

    _drive(mgr, timeout=180)   # real scanpy UMAP in a subprocess

    rec = {r.job_id: r for r in mgr.store.load_all()}[job_id]
    assert rec.status == Status.DONE.value

    cols = {c["field"]: c["datatype"] for c in project.get_datasource_metadata("cells")["columns"]}
    assert cols.get("UMAP_1") == "double"
    assert cols.get("UMAP_2") == "double"
    assert cols.get("UMAP_3") == "double"                 # n_components=3 reached scanpy end-to-end

    for col in ("UMAP_1", "UMAP_2", "UMAP_3"):            # provenance on each output
        prov = project.get_column_provenance("cells", col)
        assert prov is not None
        assert prov["job_id"] == job_id
        assert prov["params"]["n_components"] == 3        # params feed provenance identity

def test_job_that_vanishes_without_marker_is_failed(tmp_path):
    project = _make_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch",
                     executor=_FakeExecutor(poll_result="lost"))
    mgr.submit("concat_columns",
               {"datasource": "cells", "column_a": "sample", "column_b": "cluster",
                "output_name": "out"})

    mgr.tick()   # marker never written; poll says lost → the manager gives up, not waits forever

    assert [r.status for r in mgr.store.load_all()] == [Status.FAILED.value]


def test_running_job_without_marker_stays_running(tmp_path):
    project = _make_project(tmp_path)
    mgr = JobManager(project, workspace_root=tmp_path / "scratch",
                     executor=_FakeExecutor(poll_result="running"))
    mgr.submit("concat_columns",
               {"datasource": "cells", "column_a": "sample", "column_b": "cluster",
                "output_name": "out"})

    mgr.tick()   # no marker yet, but the executor says it's alive → wait, don't fail

    assert [r.status for r in mgr.store.load_all()] == [Status.RUNNING.value]


def test_unbounded_concurrency_submits_all_queued_at_once(tmp_path):
    project = _make_project(tmp_path)
    executor = _FakeExecutor(poll_result="running")     # jobs never complete
    mgr = JobManager(project, workspace_root=tmp_path / "scratch",
                     executor=executor, max_concurrent=None)

    for i in range(3):
        mgr.submit("concat_columns",
                   {"datasource": "cells", "column_a": "sample", "column_b": "cluster",
                    "output_name": f"out_{i}"})

    # no manager bound → all three submitted immediately, none held QUEUED
    statuses = [r.status for r in mgr.store.load_all()]
    assert statuses.count(Status.RUNNING.value) == 3
    assert statuses.count(Status.QUEUED.value) == 0
    assert executor.submits == 3
