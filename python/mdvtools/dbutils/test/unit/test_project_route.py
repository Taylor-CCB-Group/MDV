"""GET /projects/<id>: one project, in the same shape as an entry of GET /projects."""

import json

import pytest
from flask import Flask

from mdvtools.dbutils.dbmodels import Project, User, UserProject, db
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
def auth_app(db_app, tmp_path):
    """User 7 can read project 1 and holds a role without read access on project 3.
    Project 1's owner is someone else."""
    path = tmp_path / "3"
    write_project_files(path, view_image="img-3")
    db.session.add(Project(id=3, name="project-3", path=str(path)))
    db.session.add_all([
        User(id=7, email="reader@example.org", auth_id="reader"),
        User(id=8, email="owner@example.org", auth_id="owner"),
    ])
    db.session.add_all([
        UserProject(7, 1, can_read=True),
        UserProject(7, 3, can_read=False),
        UserProject(8, 1, is_owner=True),
    ])
    db.session.commit()

    db_app.secret_key = "test"
    register_routes(db_app, ENABLE_AUTH=True)
    return db_app


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
    assert body["owner"] == ["owner@example.org"]
    assert client.get("/projects").get_json() == [body]

    # a role without read access, a deleted project and a missing one
    assert client.get("/projects/3").status_code == 404
    assert client.get("/projects/2").status_code == 404
    assert client.get("/projects/99").status_code == 404


def test_project_route_with_auth_needs_a_session_user(auth_app):
    assert auth_app.test_client().get("/projects/1").status_code == 401
