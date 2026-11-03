"""GET /projects/<id>: one project, in the same shape as an entry of GET /projects."""

import json

import pytest
from flask import Flask

from mdvtools.auth import authutils
from mdvtools.dbutils.dbmodels import Project, db
from mdvtools.dbutils.dbservice import ProjectService
from mdvtools.dbutils.routes import register_routes


def write_project_files(path, view_image=None):
    path.mkdir()
    (path / "state.json").write_text(json.dumps({"all_views": ["main"]}))
    views = {"main": {"viewImage": view_image}} if view_image else {"main": {}}
    (path / "views.json").write_text(json.dumps(views))
    (path / "datasources.json").write_text("[]")


@pytest.fixture()
def db_app(tmp_path, monkeypatch):
    # class-level list that other tests leave entries in
    monkeypatch.setattr(ProjectService, "failed_projects", [])
    app = Flask(__name__)
    app.config.update(
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    db.init_app(app)
    with app.app_context():
        db.create_all()
        for project_id, deleted in ((1, False), (2, True)):
            path = tmp_path / str(project_id)
            write_project_files(path, view_image=f"img-{project_id}")
            db.session.add(
                Project(id=project_id, name=f"project-{project_id}", path=str(path), is_deleted=deleted)
            )
        db.session.commit()
        yield app
        db.session.remove()
        db.drop_all()


def test_get_active_project_returns_one_entry(db_app):
    entry = ProjectService.get_active_project(1)
    assert entry is not None
    assert entry["name"] == "project-1"
    assert entry["thumbnail"] == "img-1"
    assert entry == next(p for p in ProjectService.get_active_projects() if p["id"] == 1)


def test_get_active_project_skips_deleted_and_missing(db_app):
    assert ProjectService.get_active_project(2) is None
    assert ProjectService.get_active_project(99) is None


def test_project_route_without_auth(db_app):
    register_routes(db_app, ENABLE_AUTH=False)
    client = db_app.test_client()

    response = client.get("/projects/1")
    assert response.status_code == 200
    body = response.get_json()
    assert body["id"] == 1
    assert body["name"] == "project-1"
    assert body["thumbnail"] == "img-1"
    assert body == client.get("/projects").get_json()[0]

    assert client.get("/projects/2").status_code == 404
    assert client.get("/projects/99").status_code == 404


@pytest.fixture()
def auth_app(tmp_path, monkeypatch):
    def cached(project_id):
        path = tmp_path / str(project_id)
        write_project_files(path, view_image=f"img-{project_id}")
        return {
            "id": project_id,
            "name": f"project-{project_id}",
            "path": str(path),
            "lastModified": "2026-11-03 09:00:00",
            "thumbnail": None,
            "readme": None,
        }

    read_only = {"can_read": True, "can_write": False, "is_owner": False}
    no_access = {"can_read": False, "can_write": False, "is_owner": False}
    monkeypatch.setattr(authutils, "active_projects_cache", [cached(1), cached(2)])
    monkeypatch.setattr(authutils, "user_project_cache", {7: {1: read_only, 2: no_access}})
    monkeypatch.setattr(authutils, "all_users_cache", [])

    app = Flask(__name__)
    app.secret_key = "test"
    register_routes(app, ENABLE_AUTH=True)
    return app


def test_project_route_with_auth_checks_permissions(auth_app):
    client = auth_app.test_client()
    with client.session_transaction() as session:
        session["user"] = {"id": 7}

    response = client.get("/projects/1")
    assert response.status_code == 200
    body = response.get_json()
    assert body["name"] == "project-1"
    assert body["thumbnail"] == "img-1"
    assert body["permissions"]["can_read"] is True
    assert body == client.get("/projects").get_json()[0]

    # listed in the cache, but this user can't see it
    assert client.get("/projects/2").status_code == 404
    assert client.get("/projects/99").status_code == 404


def test_project_route_with_auth_needs_a_session_user(auth_app):
    assert auth_app.test_client().get("/projects/1").status_code == 500
