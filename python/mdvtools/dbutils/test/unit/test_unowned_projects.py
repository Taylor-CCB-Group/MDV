"""Tests for a project that arrives by file copy and has no owner.

A scan has no user to attribute a project to, and the project list only shows a
user the projects they hold a permission row for, so an unowned project would
otherwise be invisible to everyone.
"""

import json

import pytest
from flask import Flask

from mdvtools.dbutils.dbmodels import Project, User, UserProject, db
from mdvtools.dbutils.mdv_server_app import grant_admins_ownership_of_unowned_projects
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


def test_every_administrator_becomes_an_owner_of_an_unowned_project(app, tmp_path):
    with app.app_context():
        project_id = add_unowned_project(tmp_path).id
        first_admin_id = add_user("admin1@example.com", "admin1", is_admin=True).id
        second_admin_id = add_user("admin2@example.com", "admin2", is_admin=True).id
        add_user("reader@example.com", "reader")

        grant_admins_ownership_of_unowned_projects()

        owners = {row.user_id for row in UserProject.query.filter_by(project_id=project_id)}
        assert owners == {first_admin_id, second_admin_id}
        assert all(row.is_owner for row in UserProject.query.all())


def test_a_project_that_already_has_an_owner_is_left_alone(app, tmp_path):
    with app.app_context():
        project_id = add_unowned_project(tmp_path).id
        owner_id = add_user("owner@example.com", "owner").id
        add_user("admin@example.com", "admin", is_admin=True)
        db.session.add(UserProject(owner_id, project_id, is_owner=True))
        db.session.commit()

        grant_admins_ownership_of_unowned_projects()

        assert {row.user_id for row in UserProject.query.all()} == {owner_id}


def test_a_project_in_the_recycle_bin_is_not_given_to_administrators(app, tmp_path):
    with app.app_context():
        project = add_unowned_project(tmp_path)
        project.is_deleted = True
        db.session.commit()
        add_user("admin@example.com", "admin", is_admin=True)

        grant_admins_ownership_of_unowned_projects()

        assert UserProject.query.count() == 0


def test_rescan_gives_administrators_a_project_that_was_left_without_an_owner(app, tmp_path, monkeypatch):
    """A project that got its catalog row at a boot with no administrator present
    keeps no owner, and the scan creates no row for it the second time, so a
    rescan has to pick it up on its own."""
    from mdvtools.auth import authutils
    from mdvtools.dbutils.routes import register_routes

    with app.app_context():
        project_id = add_unowned_project(tmp_path).id
        administrator_id = add_user("admin@example.com", "admin", is_admin=True).id

    monkeypatch.setattr(authutils, "user_project_cache", {})
    monkeypatch.setattr(authutils, "active_projects_cache", [])
    monkeypatch.setattr(authutils, "user_cache", {})
    monkeypatch.setattr(authutils, "all_users_cache", [])
    with app.app_context():
        register_routes(app, True)

    client = app.test_client()
    sign_in(client, administrator_id, is_admin=True)

    response = client.get("/rescan_projects")

    assert response.status_code == 302, response.get_data(as_text=True)
    with app.app_context():
        assignment = UserProject.query.filter_by(project_id=project_id, user_id=administrator_id).one()
        assert assignment.is_owner


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


def test_administrator_manages_sharing_on_a_project_they_do_not_own(app, client, tmp_path):
    with app.app_context():
        project_id = add_unowned_project(tmp_path).id
        administrator_id = add_user("admin@example.com", "admin", is_admin=True).id
        member_id = add_user("member@example.com", "member").id

    sign_in(client, administrator_id, is_admin=True)

    listing = client.get(f"/projects/{project_id}/share")
    assert listing.status_code == 200, listing.get_data(as_text=True)
    assert listing.json["shared_users"] == []

    added = client.post(f"/projects/{project_id}/share", json={"user_id": member_id, "permission": "owner"})
    assert added.status_code == 200, added.get_data(as_text=True)
    with app.app_context():
        assert UserProject.query.filter_by(project_id=project_id, user_id=member_id).one().is_owner

    edited = client.post(f"/projects/{project_id}/share/{member_id}/edit", json={"permission": "view"})
    assert edited.status_code == 200, edited.get_data(as_text=True)
    with app.app_context():
        assert not UserProject.query.filter_by(project_id=project_id, user_id=member_id).one().is_owner

    removed = client.post(f"/projects/{project_id}/share/{member_id}/delete")
    assert removed.status_code == 200, removed.get_data(as_text=True)
    with app.app_context():
        assert UserProject.query.filter_by(project_id=project_id, user_id=member_id).count() == 0


def test_user_who_is_neither_owner_nor_administrator_cannot_manage_sharing(app, client, tmp_path):
    with app.app_context():
        project_id = add_unowned_project(tmp_path).id
        user_id = add_user("reader@example.com", "reader").id

    sign_in(client, user_id, is_admin=False)

    assert client.get(f"/projects/{project_id}/share").status_code == 403
    assert client.post(
        f"/projects/{project_id}/share", json={"user_id": user_id, "permission": "owner"}
    ).status_code == 403
    assert client.post(
        f"/projects/{project_id}/share/{user_id}/edit", json={"permission": "owner"}
    ).status_code == 403
    assert client.post(f"/projects/{project_id}/share/{user_id}/delete").status_code == 403
    with app.app_context():
        assert UserProject.query.count() == 0
