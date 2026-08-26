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
