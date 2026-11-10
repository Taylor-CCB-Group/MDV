"""The share dialog leaves out Deactivated users (ADR 0020).

A Deactivated user keeps their Project roles, so reactivating them restores
access. Until then they appear neither among a Project's members nor in the list
of people it can be shared with, and the share routes refuse to add them or
change their role.

Requests go to ProjectManagerExtension's real routes against an in-memory
database, with the signed-in Owner put straight into the session. The auth hook
is left out because these tests do not change sign-in.
"""

import pytest
from flask import Flask

from mdvtools.auth.authutils import cache_user_projects
from mdvtools.dbutils.dbmodels import Project, User, UserProject, db
from mdvtools.dbutils.dbservice import UserProjectService
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension


@pytest.fixture()
def app():
    app = Flask(__name__)
    app.secret_key = "test"
    app.config.update(
        ENABLE_AUTH=True,
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    db.init_app(app)
    ProjectManagerExtension().register_global_routes(app, {})

    with app.app_context():
        db.create_all()
        yield app
        db.session.remove()
        db.drop_all()


def add_user(email, *, is_active=True):
    user = User(email=email, auth_id=f"auth0|{email}", is_active=is_active)
    db.session.add(user)
    db.session.commit()
    return user


@pytest.fixture()
def owner(app):
    return add_user("owner@example.org")


@pytest.fixture()
def project(app, owner):
    project = Project(name="Atlas", path="/projects/atlas")
    db.session.add(project)
    db.session.commit()
    UserProjectService.add_or_update_user_project(owner.id, project.id, is_owner=True, can_write=True)
    return project


def client_for(app, user):
    # The share routes read the Owner's role from the permission cache.
    cache_user_projects()
    client = app.test_client()
    with client.session_transaction() as session:
        session["user"] = {"id": user.id, "email": user.email, "is_admin": user.is_admin}
    return client


def test_members_leave_out_a_deactivated_user(app, owner, project):
    active = add_user("ann@example.org")
    deactivated = add_user("bo@example.org", is_active=False)
    for user in (active, deactivated):
        UserProjectService.add_or_update_user_project(user.id, project.id, is_owner=False, can_write=False)

    response = client_for(app, owner).get(f"/projects/{project.id}/share")

    assert response.status_code == 200, response.get_data(as_text=True)
    members = {u["email"] for u in response.get_json()["shared_users"]}
    assert members == {"owner@example.org", "ann@example.org"}


def test_people_to_share_with_leave_out_a_deactivated_user(app, owner, project):
    add_user("ann@example.org")
    add_user("bo@example.org", is_active=False)

    response = client_for(app, owner).get(f"/projects/{project.id}/share")

    assert response.status_code == 200, response.get_data(as_text=True)
    assert [u["email"] for u in response.get_json()["all_users"]] == ["ann@example.org"]


def test_sharing_with_a_deactivated_user_is_refused(app, owner, project):
    deactivated = add_user("bo@example.org", is_active=False)

    response = client_for(app, owner).post(
        f"/projects/{project.id}/share", json={"user_id": deactivated.id, "permission": "view"}
    )

    assert response.status_code == 400
    assert "deactivated" in response.get_json()["error"].lower()
    assert UserProject.query.filter_by(user_id=deactivated.id).first() is None


def test_sharing_with_an_active_user_grants_the_role(app, owner, project):
    active = add_user("ann@example.org")

    response = client_for(app, owner).post(
        f"/projects/{project.id}/share", json={"user_id": active.id, "permission": "view"}
    )

    assert response.status_code == 200, response.get_data(as_text=True)
    granted = UserProject.query.filter_by(user_id=active.id, project_id=project.id).one()
    assert (granted.can_read, granted.can_write, granted.is_owner) == (True, False, False)


def test_changing_a_deactivated_users_role_is_refused(app, owner, project):
    deactivated = add_user("bo@example.org", is_active=False)
    UserProjectService.add_or_update_user_project(deactivated.id, project.id, is_owner=False, can_write=False)

    response = client_for(app, owner).post(
        f"/projects/{project.id}/share/{deactivated.id}/edit", json={"permission": "owner"}
    )

    assert response.status_code == 400
    assert "deactivated" in response.get_json()["error"].lower()
    db.session.expire_all()
    kept = UserProject.query.filter_by(user_id=deactivated.id, project_id=project.id).one()
    assert (kept.can_write, kept.is_owner) == (False, False)


def test_changing_an_active_users_role_updates_it(app, owner, project):
    active = add_user("ann@example.org")
    UserProjectService.add_or_update_user_project(active.id, project.id, is_owner=False, can_write=False)

    response = client_for(app, owner).post(
        f"/projects/{project.id}/share/{active.id}/edit", json={"permission": "edit"}
    )

    assert response.status_code == 200, response.get_data(as_text=True)
    db.session.expire_all()
    changed = UserProject.query.filter_by(user_id=active.id, project_id=project.id).one()
    assert (changed.can_write, changed.is_owner) == (True, False)
