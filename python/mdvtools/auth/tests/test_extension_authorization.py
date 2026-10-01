"""Tests for the admin gate that extensions sit behind.

Run against a real in-memory database rather than a mocked one, because the
behaviour under test is precisely that the gate asks the database instead of
believing the session cookie.
"""

import pytest
from flask import Flask

from mdvtools.auth.authutils import admin_required
from mdvtools.dbutils.dbmodels import User, db


def create_app(enable_auth: bool) -> Flask:
    app = Flask(__name__)
    app.secret_key = "test"
    app.config.update(
        ENABLE_AUTH=enable_auth,
        LOGIN_REDIRECT_URL="/login_dev",
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    db.init_app(app)

    @app.get("/protected-extension-route")
    @admin_required
    def protected_extension_route():
        return {"status": "allowed"}

    return app


@pytest.fixture()
def app():
    app = create_app(enable_auth=True)
    with app.app_context():
        db.create_all()
        yield app
        db.session.remove()
        db.drop_all()


def add_user(user_id=1, *, is_admin=False, is_active=True):
    user = User(
        id=user_id,
        email=f"user{user_id}@example.org",
        auth_id=f"auth0|{user_id}",
        is_admin=is_admin,
        is_active=is_active,
        password="",
    )
    db.session.add(user)
    db.session.commit()
    return user


def signed_in_as(app, user_id, **claims):
    """A client holding a session cookie, whose contents the caller chooses.

    Claims are deliberately settable so a test can assert what happens when the
    cookie says one thing and the database says another - which is the whole
    point of this gate, and cannot be expressed by logging in for real.
    """
    client = app.test_client()
    with client.session_transaction() as flask_session:
        flask_session["user"] = {"id": user_id, **claims}
    return client


def test_admin_required_allows_auth_disabled_development() -> None:
    response = (
        create_app(enable_auth=False).test_client().get("/protected-extension-route")
    )

    assert response.status_code == 200
    assert response.get_json() == {"status": "allowed"}


def test_admin_required_redirects_anonymous_user_to_mdv_login(app) -> None:
    response = app.test_client().get("/protected-extension-route")

    assert response.status_code == 302
    assert response.headers["Location"].endswith("/login_dev")


def test_admin_required_forbids_authenticated_non_admin(app) -> None:
    add_user(1, is_admin=False)
    client = signed_in_as(app, 1, is_admin=False)

    assert client.get("/protected-extension-route").status_code == 403


def test_admin_required_allows_authenticated_admin(app) -> None:
    add_user(1, is_admin=True)
    client = signed_in_as(app, 1, is_admin=True)

    assert client.get("/protected-extension-route").status_code == 200


class TestTheSessionIsNotTheAuthority:
    """A signed cookie cannot be forged, but it can be out of date.

    Everything here describes one person: somebody who was an administrator when
    they signed in, and is not one now. Without these checks they keep working
    until the cookie goes away, which means revoking access does not revoke it.
    """

    def test_a_demoted_admin_is_refused_despite_the_cookie(self, app) -> None:
        add_user(1, is_admin=False)
        client = signed_in_as(app, 1, is_admin=True)

        assert client.get("/protected-extension-route").status_code == 403

    def test_a_deactivated_user_is_treated_as_signed_out(self, app) -> None:
        add_user(1, is_admin=True, is_active=False)
        client = signed_in_as(app, 1, is_admin=True)

        response = client.get("/protected-extension-route")
        assert response.status_code == 302
        assert response.headers["Location"].endswith("/login_dev")

    def test_a_deleted_user_is_treated_as_signed_out(self, app) -> None:
        """The bootstrap rollback deletes a user who may already hold a session."""
        client = signed_in_as(app, 404, is_admin=True)

        response = client.get("/protected-extension-route")
        assert response.status_code == 302
        assert response.headers["Location"].endswith("/login_dev")

    def test_a_session_naming_nobody_is_refused(self, app) -> None:
        client = app.test_client()
        with client.session_transaction() as flask_session:
            flask_session["user"] = {"is_admin": True}

        assert client.get("/protected-extension-route").status_code == 302

    def test_a_dummy_provider_admin_is_left_alone(self, app) -> None:
        """The dummy provider invents a user with no row behind it. Looking its id
        up would lock local development out, or silently borrow the privileges of
        whoever really holds that id."""
        add_user(1, is_admin=False)
        app.config["DEFAULT_AUTH_METHOD"] = "dummy"
        client = signed_in_as(app, 1, is_admin=True)

        assert client.get("/protected-extension-route").status_code == 200

    def test_a_promoted_user_does_not_have_to_sign_in_again(self, app) -> None:
        """The refresh runs in both directions - it is not a deny-list."""
        add_user(1, is_admin=True)
        client = signed_in_as(app, 1, is_admin=False)

        assert client.get("/protected-extension-route").status_code == 200
