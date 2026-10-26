import pandas as pd
import pytest

from mdvtools.mdvproject import MDVProject
from mdvtools.server import build_app
from mdvtools.server_extension import MDVServerOptions
from mdvtools.jobs.registry import serialize_registry

import mdvtools.server as server_module
from mdvtools.jobs.service import JobService

@pytest.fixture(autouse=True)
def _fresh_job_service():
    # project.id is the dir basename, so tests reusing "proj" collide through the
    # module-global service; a fresh service per test clears its manager cache.
    server_module.job_service = JobService()
    yield

def _project_with_cells(tmp_path):
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    df = pd.DataFrame({"sample": ["s1", "s2", "s3"], "cluster": ["a", "b", "a"]})
    project.add_datasource("cells", df)
    return project


def test_post_jobs_queues_and_returns_202(tmp_path):
    """POST /jobs writes a QUEUED record and returns 202 {job_id}. The driver isn't
    started in this test, so the record stays QUEUED (submit is write-ahead only)."""
    from pathlib import Path
    from mdvtools.jobs.jobstore import JobStore, Status
    from mdvtools.jobs import JOBS_DIRNAME

    project = _project_with_cells(tmp_path)
    app = build_app(project, MDVServerOptions(open_browser=False, websocket=False))

    resp = app.test_client().post(
        "/jobs",
        json={
            "tool_id": "concat_columns",
            "params": {
                "datasource": "cells",
                "column_a": "sample",
                "column_b": "cluster",
                "output_name": "out",
            },
        },
    )

    assert resp.status_code == 202
    job_id = resp.get_json()["job_id"]
    assert job_id

    recs = JobStore(Path(project.dir) / JOBS_DIRNAME).load_all()
    assert [r.status for r in recs] == [Status.QUEUED.value]
    assert recs[0].job_id == job_id


def test_post_jobs_rejects_unknown_tool_with_400(tmp_path):
    """Backend re-validates (ADR-0006). An unknown tool raises KeyError in submit,
    which the route turns into 400, not a 500."""
    project = _project_with_cells(tmp_path)
    app = build_app(project, MDVServerOptions(open_browser=False, websocket=False))

    resp = app.test_client().post(
        "/jobs", json={"tool_id": "does_not_exist", "params": {}}
    )

    assert resp.status_code == 400

def test_jobs_tools_route_answers_over_http(tmp_path):
    """build_app registers routes and returns the app without serving, so a Flask
    test client can check the route."""
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    app = build_app(project, MDVServerOptions(open_browser=False, websocket=False))

    resp = app.test_client().get("/jobs/tools")

    assert resp.status_code == 200
    assert resp.get_json() == serialize_registry()

def test_start_driver_reconciles_before_starting():
    """Boot wiring seam: recovery_scan must run before start(), so no tick fires
    against a project whose manager the scan is still building (ADR-0012)."""
    from mdvtools.server import start_driver

    calls = []

    class SpyService:
        def recovery_scan(self, projects):
            calls.append(("scan", list(projects)))

        def start(self):
            calls.append(("start",))

    start_driver(SpyService(), ["p1", "p2"])

    assert calls == [("scan", ["p1", "p2"]), ("start",)]


def test_get_jobs_lists_records_without_handle(tmp_path):
    from mdvtools.jobs.jobstore import Status

    project = _project_with_cells(tmp_path)
    app = build_app(project, MDVServerOptions(open_browser=False, websocket=False))
    client = app.test_client()
    job_id = client.post(
        "/jobs",
        json={
            "tool_id": "concat_columns",
            "params": {
                "datasource": "cells",
                "column_a": "sample",
                "column_b": "cluster",
                "output_name": "out",
            },
        },
    ).get_json()["job_id"]

    # stand in for the driver having dispatched it, so the record carries an internal handle
    store = server_module.job_service.get_or_create(project).store
    rec = {r.job_id: r for r in store.load_all()}[job_id]
    store.set(rec, Status.RUNNING, handle={"kind": "local", "ref": "123"})

    resp = client.get("/jobs")

    assert resp.status_code == 200
    jobs = resp.get_json()
    assert [j["job_id"] for j in jobs] == [job_id]
    assert jobs[0]["status"] == Status.RUNNING.value
    assert "handle" not in jobs[0]


def test_get_job_by_id_returns_client_safe_record(tmp_path):
    from mdvtools.jobs.jobstore import Status

    project = _project_with_cells(tmp_path)
    app = build_app(project, MDVServerOptions(open_browser=False, websocket=False))
    client = app.test_client()
    job_id = client.post(
        "/jobs",
        json={
            "tool_id": "concat_columns",
            "params": {
                "datasource": "cells",
                "column_a": "sample",
                "column_b": "cluster",
                "output_name": "out",
            },
        },
    ).get_json()["job_id"]

    resp = client.get(f"/jobs/{job_id}")

    assert resp.status_code == 200
    job = resp.get_json()
    assert job["job_id"] == job_id
    assert job["status"] == Status.QUEUED.value
    assert "handle" not in job


def test_get_job_by_unknown_id_returns_404(tmp_path):
    project = _project_with_cells(tmp_path)
    app = build_app(project, MDVServerOptions(open_browser=False, websocket=False))

    resp = app.test_client().get("/jobs/does_not_exist")

    assert resp.status_code == 404
    assert resp.get_json()["error"]


def test_jobs_health_reports_quarantined_count(tmp_path):
    from pathlib import Path
    from mdvtools.jobs import JOBS_DIRNAME

    project = _project_with_cells(tmp_path)
    records = Path(project.dir) / JOBS_DIRNAME / "records"
    records.mkdir(parents=True, exist_ok=True)
    (records / "broken.json").write_text("{not json")
    app = build_app(project, MDVServerOptions(open_browser=False, websocket=False))

    resp = app.test_client().get("/jobs/health")

    assert resp.status_code == 200
    assert resp.get_json() == {"quarantined": 1}


def _app_offering_slurm(tmp_path):
    """An app whose job service offers local and slurm, as MDV_JOBS_SLURM=cli would (ADR-0016)."""
    from mdvtools.jobs.service import manager_factory_from_env

    server_module.job_service = JobService(manager_factory=manager_factory_from_env(
        {"MDV_JOBS_SLURM": "cli", "MDV_JOBS_SLURM_WORKSPACE_ROOT": str(tmp_path / "shared")}
    ))
    project = _project_with_cells(tmp_path)
    return build_app(project, MDVServerOptions(open_browser=False, websocket=False))


_CONCAT_PARAMS = {"datasource": "cells", "column_a": "sample", "column_b": "cluster", "output_name": "out"}


def test_get_jobs_executors_lists_the_offered_names(tmp_path):
    app = _app_offering_slurm(tmp_path)

    resp = app.test_client().get("/jobs/executors")

    assert resp.status_code == 200
    assert resp.get_json() == ["local", "slurm"]    # first is the default


def test_post_jobs_with_executor_stores_it_and_the_client_view_shows_it(tmp_path):
    client = _app_offering_slurm(tmp_path).test_client()

    resp = client.post("/jobs", json={
        "tool_id": "concat_columns", "params": _CONCAT_PARAMS, "executor": "slurm",
    })

    assert resp.status_code == 202
    job = client.get(f"/jobs/{resp.get_json()['job_id']}").get_json()
    assert job["executor"] == "slurm"
    assert "handle" not in job


def test_post_jobs_with_unoffered_executor_returns_400(tmp_path):
    client = _app_offering_slurm(tmp_path).test_client()

    resp = client.post("/jobs", json={
        "tool_id": "concat_columns", "params": _CONCAT_PARAMS, "executor": "k8s",
    })

    assert resp.status_code == 400
    assert "k8s" in resp.get_json()["error"]
    assert client.get("/jobs").get_json() == []     # nothing queued


def test_post_jobs_with_non_string_executor_returns_400(tmp_path):
    client = _app_offering_slurm(tmp_path).test_client()

    resp = client.post("/jobs", json={
        "tool_id": "concat_columns", "params": _CONCAT_PARAMS, "executor": 5,
    })

    assert resp.status_code == 400


def _app_with_matrix(tmp_path):
    import numpy as np
    import scipy.sparse

    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    project.add_datasource("cells", pd.DataFrame({"cell_id": ["c0", "c1", "c2"]}))
    project.add_datasource("genes", pd.DataFrame({"name": ["g0", "g1"]}))
    project.add_rows_as_columns_link("cells", "genes", "name", "Gene Expr")
    X = scipy.sparse.csc_matrix(np.array([[1.0, 0.0], [0.0, 2.0], [3.0, 4.0]], dtype=np.float32))
    project.add_rows_as_columns_subgroup("cells", "genes", "gs", X, name="gene_scores")
    return build_app(project, MDVServerOptions(open_browser=False, websocket=False))


_UMAP_ON_SELECTION = {"datasource": "cells", "layer": "gs", "selection": "T cells", "output_name": "UMAP_T_cells"}


def test_post_jobs_with_rows_records_the_selection_and_the_client_view_shows_its_size(tmp_path):
    client = _app_with_matrix(tmp_path).test_client()

    resp = client.post("/jobs", json={"tool_id": "umap", "params": _UMAP_ON_SELECTION, "rows": [2, 0]})

    assert resp.status_code == 202
    job = client.get(f"/jobs/{resp.get_json()['job_id']}").get_json()
    assert job["selection_rows"] == 2
    assert job["input_filter_hash"]
    assert "rows" not in job      # the indices stay in the sidecar


@pytest.mark.parametrize("rows", ["0,2", {"0": 1}, 3])
def test_post_jobs_with_rows_that_are_not_a_list_returns_400(tmp_path, rows):
    client = _app_with_matrix(tmp_path).test_client()

    resp = client.post("/jobs", json={"tool_id": "umap", "params": _UMAP_ON_SELECTION, "rows": rows})

    assert resp.status_code == 400
    assert "rows must be a list" in resp.get_json()["error"]


def test_post_jobs_with_rows_outside_the_datasource_returns_400_and_queues_nothing(tmp_path):
    client = _app_with_matrix(tmp_path).test_client()

    resp = client.post("/jobs", json={"tool_id": "umap", "params": _UMAP_ON_SELECTION, "rows": [0, 3]})

    assert resp.status_code == 400
    assert "outside" in resp.get_json()["error"]
    assert client.get("/jobs").get_json() == []
