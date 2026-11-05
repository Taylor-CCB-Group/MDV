"""A sign-in with no User accepts an open Invitation (ADRs 0019 and 0021).

The Invitation is accepted only when the email matches, the sign-in came through
the Invitation's Sign-in method, and Auth0 vouches for the email: email_verified
for a password connection, Microsoft's xms_edov claim for Azure AD. A post-login
Action adds the connection and its strategy to what Auth0 returns.

Requests go through the real auth hook and Auth0Provider against an in-memory
database. Only the calls that would reach Auth0 are stubbed.
"""

from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from flask import Flask

import mdvtools
from mdvtools.auth.auth0_provider import Auth0Provider
from mdvtools.auth.authutils import is_authenticated, register_before_request_auth
from mdvtools.dbutils.dbmodels import Invitation, InvitationProject, Project, User, UserProject, db
from mdvtools.dbutils.dbservice import UserProjectService

PASSWORD = "Username-Password-Authentication"
MICROSOFT = "acme-entra"


@pytest.fixture()
def app():
    app = Flask(__name__, template_folder=str(Path(mdvtools.__file__).parent / "templates"))
    app.secret_key = "test"
    app.config.update(
        ENABLE_AUTH=True,
        DEFAULT_AUTH_METHOD="auth0",
        LOGIN_REDIRECT_URL="/login_dev",
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    db.init_app(app)
    register_before_request_auth(app)

    @app.route("/projects/<int:project_id>/open")
    def open_project(project_id):
        # Stands in for the project routes, which check access the same way.
        ok, reason = is_authenticated(project_id=project_id)
        return ("opened", 200) if ok else (str(reason), 403)

    with app.app_context():
        db.create_all()
        yield app
        db.session.remove()
        db.drop_all()


@pytest.fixture()
def provider(app):
    with patch("mdvtools.auth.auth0_provider.requests.get") as metadata:
        metadata.return_value.status_code = 200
        metadata.return_value.json.return_value = {
            "jwks_uri": "https://test.auth0.com/.well-known/jwks.json"
        }
        return Auth0Provider(
            app,
            oauth=MagicMock(),
            client_id="client",
            client_secret="secret",
            domain="test.auth0.com",
        )


@pytest.fixture()
def admin(app):
    user = User(email="admin@example.org", auth_id="auth0|admin", is_active=True, is_admin=True)
    db.session.add(user)
    db.session.commit()
    return user


@pytest.fixture()
def project(app):
    project = Project(name="Atlas", path="/projects/atlas")
    db.session.add(project)
    db.session.commit()
    return project


def invite(admin, email, connection, roles, **fields):
    """Write an Invitation the way MDV-Admin does: straight to the database."""
    invitation = Invitation(
        email=email, connection=connection, invited_by=admin.id, sent_at=datetime.now(), **fields
    )
    for project, role in roles.items():
        invitation.projects.append(InvitationProject(project_id=project.id, role=role))
    db.session.add(invitation)
    db.session.commit()
    return invitation


def password_sign_in(sub, email, *, verified=True, connection=PASSWORD):
    return {
        "sub": sub,
        "email": email,
        "email_verified": verified,
        "connection": connection,
        "connection_strategy": "auth0",
    }


def get_as(app, provider, path, sign_in):
    """Request path as someone Auth0 has just signed in with these claims."""
    with patch.multiple(
        provider,
        get_token=MagicMock(return_value="token"),
        is_token_valid=MagicMock(return_value=True),
        get_user=MagicMock(return_value=sign_in),
    ), patch("mdvtools.auth.authutils.get_auth_provider", return_value=provider):
        return app.test_client().get(path)


def test_new_person_with_a_password_invitation_opens_the_invited_project(app, provider, admin, project):
    invite(admin, "ann@example.org", PASSWORD, {project: "viewer"})

    response = get_as(
        app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|ann", "Ann@Example.org")
    )

    assert response.status_code == 200, response.get_data(as_text=True)


def assert_refused(response, invitation):
    assert response.status_code == 403
    assert "No access to this MDV site" in response.get_data(as_text=True)
    db.session.expire_all()
    assert invitation.is_open()


def test_password_sign_in_with_an_unverified_email_is_refused(app, provider, admin, project):
    # ADR 0021: on a password connection that allows self-signup, anyone can
    # register an account with someone else's email.
    invitation = invite(admin, "ann@example.org", PASSWORD, {project: "viewer"})

    response = get_as(
        app,
        provider,
        f"/projects/{project.id}/open",
        password_sign_in("auth0|imposter", "ann@example.org", verified=False),
    )

    assert_refused(response, invitation)


def test_sign_in_through_another_connection_is_refused(app, provider, admin, project):
    invitation = invite(admin, "ann@example.org", MICROSOFT, {project: "viewer"})

    response = get_as(
        app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|ann", "ann@example.org")
    )

    assert_refused(response, invitation)


def microsoft_sign_in(sub, email, *, xms_edov=True):
    return {
        "sub": sub,
        "email": email,
        # Entra's email claim is unverified (ADR 0021), so this says nothing.
        "email_verified": False,
        "connection": MICROSOFT,
        "connection_strategy": "waad",
        "xms_edov": xms_edov,
    }


def test_microsoft_sign_in_from_a_verified_domain_opens_the_invited_project(app, provider, admin, project):
    invite(admin, "bo@acme.example", MICROSOFT, {project: "viewer"})

    response = get_as(
        app, provider, f"/projects/{project.id}/open", microsoft_sign_in("waad|bo", "bo@acme.example")
    )

    assert response.status_code == 200, response.get_data(as_text=True)


def test_microsoft_sign_in_without_a_verified_domain_is_refused(app, provider, admin, project):
    invitation = invite(admin, "bo@acme.example", MICROSOFT, {project: "viewer"})

    response = get_as(
        app,
        provider,
        f"/projects/{project.id}/open",
        microsoft_sign_in("waad|bo", "bo@acme.example", xms_edov=False),
    )

    assert_refused(response, invitation)


def test_sign_in_through_a_connection_of_another_kind_is_refused(app, provider, admin, project):
    # Only password and Azure AD connections have a rule for vouching for the email.
    invitation = invite(admin, "ann@example.org", "google-oauth2", {project: "viewer"})
    sign_in = password_sign_in("google-oauth2|ann", "ann@example.org", connection="google-oauth2")
    sign_in["connection_strategy"] = "google-oauth2"

    response = get_as(app, provider, f"/projects/{project.id}/open", sign_in)

    assert_refused(response, invitation)


@pytest.fixture()
def second_project(app):
    project = Project(name="Survey", path="/projects/survey")
    db.session.add(project)
    db.session.commit()
    return project


def test_invitation_naming_a_user_moves_them_to_the_new_sign_in_method(
    app, provider, admin, project, second_project
):
    # ADR 0019: an Invitation sent to an existing User moves them to a new
    # Sign-in method and can carry a new email.
    user = User(email="bo@old.example", auth_id="auth0|bo", is_active=True)
    db.session.add(user)
    db.session.commit()
    db.session.add(UserProject(user_id=user.id, project_id=project.id, is_owner=True))
    db.session.commit()
    invite(admin, "bo@acme.example", MICROSOFT, {second_project: "viewer"}, user_id=user.id)

    new_account = microsoft_sign_in("waad|bo", "bo@acme.example")
    assert get_as(app, provider, f"/projects/{project.id}/open", new_account).status_code == 200
    assert get_as(app, provider, f"/projects/{second_project.id}/open", new_account).status_code == 200

    old_account = password_sign_in("auth0|bo", "bo@old.example")
    assert get_as(app, provider, f"/projects/{project.id}/open", old_account).status_code == 403


def test_deactivated_user_with_an_invitation_is_reactivated_with_its_roles(
    app, provider, admin, project, second_project
):
    # ADR 0019: inviting a Deactivated user reactivates them with their earlier
    # Project roles, and a role chosen in the Invitation replaces the earlier one.
    user = User(email="cy@example.org", auth_id="auth0|cy", is_active=False)
    db.session.add(user)
    db.session.commit()
    db.session.add(UserProject(user_id=user.id, project_id=project.id, is_owner=True))
    db.session.add(UserProject(user_id=user.id, project_id=second_project.id, can_read=True))
    db.session.commit()
    invite(admin, "cy@example.org", PASSWORD, {project: "viewer"}, user_id=user.id)

    sign_in = password_sign_in("auth0|cy", "cy@example.org")
    assert get_as(app, provider, f"/projects/{project.id}/open", sign_in).status_code == 200
    assert get_as(app, provider, f"/projects/{second_project.id}/open", sign_in).status_code == 200

    viewer = {"can_read": True, "can_write": False, "is_owner": False}
    assert UserProjectService.get_user_project_permissions(user.id, project.id) == viewer


def test_deactivated_user_cannot_accept_an_invitation_that_does_not_name_them(
    app, provider, admin, project
):
    db.session.add(User(email="cy@example.org", auth_id="auth0|cy", is_active=False))
    db.session.commit()
    invitation = invite(admin, "cy@example.org", PASSWORD, {project: "viewer"})

    response = get_as(
        app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|cy", "cy@example.org")
    )

    assert_refused(response, invitation)


def test_new_account_cannot_take_the_email_of_a_user_the_invitation_does_not_name(
    app, provider, admin, project
):
    # MDV-Admin names the User when it invites an email that already has one,
    # so an Invitation that names nobody must not reach that User.
    db.session.add(User(email="Dee@Example.org", auth_id="auth0|dee", is_active=True))
    db.session.commit()
    invitation = invite(admin, "dee@example.org", MICROSOFT, {project: "viewer"})

    response = get_as(
        app, provider, f"/projects/{project.id}/open", microsoft_sign_in("waad|dee", "dee@example.org")
    )

    assert_refused(response, invitation)
