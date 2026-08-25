from mdvtools.mdvproject import MDVProject
from mdvtools.server import build_app
from mdvtools.server_extension import MDVServerOptions
from mdvtools.jobs.registry import serialize_registry


def test_jobs_tools_route_answers_over_http(tmp_path):
    """build_app registers routes and returns the app without serving, so a Flask
    test client can check the route."""
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    app = build_app(project, MDVServerOptions(open_browser=False, websocket=False))

    resp = app.test_client().get("/jobs/tools")

    assert resp.status_code == 200
    assert resp.get_json() == serialize_registry()
