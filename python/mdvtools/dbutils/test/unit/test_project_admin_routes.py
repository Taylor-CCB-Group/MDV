"""Tests for administrator access to the project management routes.

A project that arrives by file copy and is picked up at boot has no owner, so
without an administrator bypass its name and access level can only be changed
through direct database access.
"""

import json

import pytest
from flask import Flask

from mdvtools.dbutils.dbmodels import Project, User, db
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension


@pytest.fixture()
def app(tmp_path):
    app = Flask(__name__)
    app.config.update(
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
        ENABLE_AUTH=True,
        SECRET_KEY="test",
        projects_base_dir=str(tmp_path),
    )
    db.init_app(app)
    with app.app_context():
        db.create_all()
        yield app
        db.session.remove()
        db.drop_all()


@pytest.fixture()
def client(app, monkeypatch):
    """The routes read the permission caches they captured at registration, so
    the caches are replaced before the routes are registered."""
    from mdvtools.auth import authutils

    monkeypatch.setattr(authutils, "user_project_cache", {})
    monkeypatch.setattr(authutils, "active_projects_cache", [])
    monkeypatch.setattr(authutils, "user_cache", {})
    monkeypatch.setattr(authutils, "all_users_cache", [])
    with app.app_context():
        ProjectManagerExtension().register_global_routes(app, app.config)
    return app.test_client()


def add_unowned_project(tmp_path, name="orphan"):
    """A catalog row with no user_projects row, as a boot-time rescan creates."""
    path = tmp_path / name
    path.mkdir()
    (path / "datasources.json").write_text(json.dumps([]))
    (path / "views.json").write_text(json.dumps({}))
    (path / "state.json").write_text(json.dumps({"all_views": [], "permission": "edit"}))
    project = Project()
    project.name = name
    project.path = str(path)
    db.session.add(project)
    db.session.commit()
    return project


def add_user(email, auth_id, *, is_admin=False):
    user = User(email=email, auth_id=auth_id, is_admin=is_admin)
    db.session.add(user)
    db.session.commit()
    return user


def sign_in(client, user_id, *, is_admin):
    with client.session_transaction() as session:
        session["user"] = {"id": user_id, "auth_id": f"user-{user_id}", "is_admin": is_admin}


def test_administrator_renames_a_project_they_do_not_own(app, client, tmp_path):
    with app.app_context():
        project_id = add_unowned_project(tmp_path).id
        administrator_id = add_user("admin@example.com", "admin", is_admin=True).id

    sign_in(client, administrator_id, is_admin=True)

    response = client.put(f"/projects/{project_id}/rename", data={"name": "Tumour atlas"})

    assert response.status_code == 200, response.get_data(as_text=True)
    with app.app_context():
        assert db.session.get(Project, project_id).name == "Tumour atlas"


def test_user_who_is_neither_owner_nor_administrator_cannot_rename(app, client, tmp_path):
    with app.app_context():
        project_id = add_unowned_project(tmp_path).id
        user_id = add_user("reader@example.com", "reader").id

    sign_in(client, user_id, is_admin=False)

    response = client.put(f"/projects/{project_id}/rename", data={"name": "Tumour atlas"})

    assert response.status_code == 403
    with app.app_context():
        assert db.session.get(Project, project_id).name == "orphan"


def test_administrator_changes_the_access_level_of_a_project_they_do_not_own(app, client, tmp_path):
    with app.app_context():
        project_id = add_unowned_project(tmp_path).id
        administrator_id = add_user("admin@example.com", "admin", is_admin=True).id

    sign_in(client, administrator_id, is_admin=True)

    response = client.put(f"/projects/{project_id}/access", data={"type": "read-only"})

    assert response.status_code == 200, response.get_data(as_text=True)
    with app.app_context():
        assert db.session.get(Project, project_id).access_level == "read-only"


def test_user_who_is_neither_owner_nor_administrator_cannot_change_access(app, client, tmp_path):
    with app.app_context():
        project_id = add_unowned_project(tmp_path).id
        user_id = add_user("reader@example.com", "reader").id

    sign_in(client, user_id, is_admin=False)

    response = client.put(f"/projects/{project_id}/access", data={"type": "read-only"})

    assert response.status_code == 403
    with app.app_context():
        assert db.session.get(Project, project_id).access_level == "editable"
