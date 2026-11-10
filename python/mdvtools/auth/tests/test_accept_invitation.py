"""A sign-in accepts an open Invitation for its account (ADRs 0019 and 0021).

The Invitation is accepted only when the email matches, the sign-in came through
the Invitation's Sign-in method, and Auth0 vouches for the email: email_verified
for a password connection, Microsoft's xms_edov claim for a Microsoft
connection. A post-login
Action adds the connection and its strategy to what Auth0 returns.

Requests go through the real auth hook and Auth0Provider against an in-memory
database. Only the calls that would reach Auth0 are stubbed.
"""

from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from flask import Flask
from sqlalchemy import update

import mdvtools
from mdvtools.auth.auth0_provider import Auth0Provider
from mdvtools.auth.authutils import is_authenticated, register_before_request_auth
from mdvtools.dbutils.dbmodels import Invitation, InvitationProject, Project, User, UserProject, db
from mdvtools.dbutils.dbservice import InvitationService, UserProjectService

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


@pytest.mark.parametrize(
    "claims, expected",
    [
        pytest.param(
            {
                "https://mdv/connection": MICROSOFT,
                "https://mdv/connection_strategy": "oidc",
                "https://mdv/xms_edov": True,
            },
            (MICROSOFT, "oidc", True),
            id="claims from the Action",
        ),
        pytest.param({"https://mdv/xms_edov": "true"}, (None, None, False), id="xms_edov not a boolean"),
        pytest.param({}, (None, None, False), id="no Action"),
    ],
)
def test_sign_in_claims_are_read_from_userinfo(provider, claims, expected):
    # The other tests stub get_user, so this is what checks the claim names
    # MDV reads against the ones the post-login Action writes.
    userinfo = {"sub": "waad|ann", "email": "ann@acme.example", **claims}
    with patch("mdvtools.auth.auth0_provider.requests.get") as get:
        get.return_value.status_code = 200
        get.return_value.json.return_value = userinfo
        user = provider.get_user({"access_token": "token"})

    assert user is not None
    assert (user["connection"], user["connection_strategy"], user["xms_edov"]) == expected


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
    invitation = invite(admin, "ann@example.org", PASSWORD, {project: "viewer"})

    response = get_as(
        app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|ann", "Ann@Example.org")
    )

    assert response.status_code == 200, response.get_data(as_text=True)
    db.session.expire_all()
    assert invitation.accepted_at is not None
    assert invitation.accepted_user_id == User.query.filter_by(auth_id="auth0|ann").one().id
    assert Invitation.open_for_email("ann@example.org") == []


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


def microsoft_sign_in(sub, email, *, xms_edov=True, strategy="waad"):
    return {
        "sub": sub,
        "email": email,
        # Entra's email claim is unverified (ADR 0021), so this says nothing.
        "email_verified": False,
        "connection": MICROSOFT,
        "connection_strategy": strategy,
        "xms_edov": xms_edov,
    }


# A Microsoft connection is either Auth0's native Azure AD connection or an
# Enterprise OpenID Connect connection to Entra (ADR 0021).
microsoft_strategies = pytest.mark.parametrize("strategy", ["waad", "oidc"])


@microsoft_strategies
def test_microsoft_sign_in_from_a_verified_domain_opens_the_invited_project(
    app, provider, admin, project, strategy
):
    invite(admin, "bo@acme.example", MICROSOFT, {project: "viewer"})

    response = get_as(
        app,
        provider,
        f"/projects/{project.id}/open",
        microsoft_sign_in(f"{strategy}|bo", "bo@acme.example", strategy=strategy),
    )

    assert response.status_code == 200, response.get_data(as_text=True)


@microsoft_strategies
def test_microsoft_sign_in_without_a_verified_domain_is_refused(app, provider, admin, project, strategy):
    invitation = invite(admin, "bo@acme.example", MICROSOFT, {project: "viewer"})

    response = get_as(
        app,
        provider,
        f"/projects/{project.id}/open",
        microsoft_sign_in(f"{strategy}|bo", "bo@acme.example", xms_edov=False, strategy=strategy),
    )

    assert_refused(response, invitation)


def test_sign_in_through_a_connection_of_another_kind_is_refused(app, provider, admin, project):
    # Only password and Microsoft connections have a rule for vouching for the email.
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


def test_user_reactivated_before_signing_in_still_accepts_the_invitation(
    app, provider, admin, project
):
    # An Admin invites a Deactivated user, then reactivates them before they
    # sign in. Their usual account still accepts the Invitation that names them.
    user = User(email="eve@example.org", auth_id="auth0|eve", is_active=False)
    db.session.add(user)
    db.session.commit()
    db.session.add(UserProject(user_id=user.id, project_id=project.id, can_read=True))
    db.session.commit()
    invitation = invite(admin, "eve@example.org", PASSWORD, {project: "editor"}, user_id=user.id)
    user.is_active = True
    db.session.commit()

    response = get_as(
        app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|eve", "eve@example.org")
    )

    assert response.status_code == 200, response.get_data(as_text=True)
    editor = {"can_read": True, "can_write": True, "is_owner": False}
    assert UserProjectService.get_user_project_permissions(user.id, project.id) == editor
    db.session.expire_all()
    assert invitation.accepted_user_id == user.id
    assert Invitation.open_for_email("eve@example.org") == []


def test_active_user_signing_in_leaves_an_invitation_that_does_not_name_them(
    app, provider, admin, project
):
    user = User(email="fay@example.org", auth_id="auth0|fay", is_active=True)
    db.session.add(user)
    db.session.commit()
    db.session.add(UserProject(user_id=user.id, project_id=project.id, can_read=True))
    db.session.commit()
    invitation = invite(admin, "fay@example.org", PASSWORD, {project: "owner"})

    response = get_as(
        app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|fay", "fay@example.org")
    )

    assert response.status_code == 200, response.get_data(as_text=True)
    viewer = {"can_read": True, "can_write": False, "is_owner": False}
    assert UserProjectService.get_user_project_permissions(user.id, project.id) == viewer
    db.session.expire_all()
    assert invitation.accepted_at is None


def test_deactivated_admin_accepting_an_invitation_is_reactivated_without_admin_rights(
    app, provider, admin, project
):
    # Deactivating an Admin leaves is_admin set. An Invitation grants its
    # projects, so Admin rights come back only through an Admin's choice.
    user = User(email="gil@example.org", auth_id="auth0|gil", is_active=False, is_admin=True)
    db.session.add(user)
    db.session.commit()
    invite(admin, "gil@example.org", PASSWORD, {project: "viewer"}, user_id=user.id)

    response = get_as(
        app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|gil", "gil@example.org")
    )

    assert response.status_code == 200, response.get_data(as_text=True)
    db.session.expire_all()
    assert (user.is_active, user.is_admin) == (True, False)


def test_admin_reactivated_before_signing_in_keeps_admin_rights(app, provider, admin, project):
    # Reactivating them was an Admin's choice, which includes their Admin rights.
    user = User(email="hal@example.org", auth_id="auth0|hal", is_active=True, is_admin=True)
    db.session.add(user)
    db.session.commit()
    invitation = invite(admin, "hal@example.org", PASSWORD, {project: "viewer"}, user_id=user.id)

    response = get_as(
        app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|hal", "hal@example.org")
    )

    assert response.status_code == 200, response.get_data(as_text=True)
    db.session.expire_all()
    assert invitation.accepted_user_id == user.id
    assert user.is_admin is True


def test_invitation_is_refused_while_another_user_holds_its_email_in_another_case(
    app, provider, admin, project
):
    # Older Users' emails were stored as typed, so two Users can hold the same
    # email in different cases. Whichever row the database returns first, the
    # Invitation must not reach the User it does not name.
    named = User(email="ivy@example.org", auth_id="auth0|ivy", is_active=False)
    db.session.add_all([named, User(email="Ivy@Example.org", auth_id="auth0|ivy-old", is_active=True)])
    db.session.commit()
    invitation = invite(admin, "ivy@example.org", PASSWORD, {project: "viewer"}, user_id=named.id)

    response = get_as(
        app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|ivy", "ivy@example.org")
    )

    assert_refused(response, invitation)


@pytest.mark.parametrize("winner_active, same_account", [(True, True), (False, True), (True, False)])
def test_second_acceptance_after_first_commits(app, provider, admin, project, winner_active, same_account):
    invitation = invite(admin, "jo@example.org", PASSWORD, {project: "viewer"})
    sign_in = password_sign_in("auth0|jo", "jo@example.org")
    winner_sign_in = sign_in if same_account else password_sign_in("auth0|other", "jo@example.org")
    matches = InvitationService._sign_in_matches
    winners = []

    def accepted_once_matched(candidate, claims):
        matched = matches(candidate, claims)
        if matched:
            # Both requests picked the open invitation. Finish the first
            # acceptance before the second request tries its conditional update.
            with patch.object(InvitationService, "_sign_in_matches", staticmethod(matches)):
                winner = InvitationService.accept_for_sign_in(winner_sign_in)
            assert winner is not None
            winners.append(winner.id)
            if not winner_active:
                winner.is_active = False
                db.session.commit()
        return matched

    with patch.object(InvitationService, "_sign_in_matches", staticmethod(accepted_once_matched)):
        response = get_as(app, provider, f"/projects/{project.id}/open", sign_in)

    assert response.status_code == (200 if winner_active and same_account else 403)
    db.session.expire_all()
    assert invitation.accepted_at is not None
    assert invitation.accepted_user_id == winners[0]
    assert invitation.revoked_at is None
    assert User.query.filter_by(email="jo@example.org").count() == 1
    assert UserProject.query.filter_by(project_id=project.id, user_id=winners[0]).count() == 1
    assert User.query.filter_by(id=winners[0]).one().is_active == winner_active


def test_invitation_revoked_while_its_sign_in_is_accepted_grants_nothing(app, provider, admin, project):
    # An Admin revokes the Invitation after MDV has picked it but before MDV
    # writes the acceptance. The plain SQL update commits the revoke without
    # changing the copy MDV already loaded, as a revoke from MDV-Admin would.
    invitation = invite(admin, "jo@example.org", PASSWORD, {project: "viewer"})
    matches = InvitationService._sign_in_matches

    def revoked_once_matched(candidate, sign_in):
        matched = matches(candidate, sign_in)
        if matched:
            db.session.execute(
                update(Invitation).where(Invitation.id == candidate.id).values(revoked_at=datetime.now())
            )
            db.session.commit()
        return matched

    with patch.object(InvitationService, "_sign_in_matches", staticmethod(revoked_once_matched)):
        response = get_as(
            app, provider, f"/projects/{project.id}/open", password_sign_in("auth0|jo", "jo@example.org")
        )

    assert response.status_code == 403
    db.session.expire_all()
    assert User.query.filter_by(email="jo@example.org").first() is None
    assert UserProject.query.filter_by(project_id=project.id).count() == 0
    assert (invitation.accepted_at, invitation.revoked_at is not None) == (None, True)
