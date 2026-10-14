"""Regression tests for mdv-roadmap issue 29.

A session cookie issued by one MDV deployment must not carry that deployment's
local user id or admin flag into another deployment on the same host. The auth
hook re-resolves session["user"] from this deployment's own users on every
request, keyed by the provider's auth_id.
"""

from types import SimpleNamespace
from unittest.mock import patch

import pytest
from flask import Flask, jsonify, session

from mdvtools.auth.authutils import register_before_request_auth
from mdvtools.auth.dummy_provider import DummyAuthProvider

# The same Auth0 account as recorded by the deployment that issued the cookie,
# where it happens to be user 2 and an administrator.
FOREIGN_SESSION_USER = {"id": 2, "auth_id": "auth0|someone", "email": "someone@example.com", "is_admin": True}
# That account's row in this deployment's database.
LOCAL_USER = {"id": 7, "auth_id": "auth0|someone", "email": "someone@example.com", "is_admin": False}


@pytest.fixture
def app():
    app = Flask("test_session_identity")
    app.config.update(
        TESTING=True,
        ENABLE_AUTH=True,
        DEFAULT_AUTH_METHOD="auth0",
        LOGIN_REDIRECT_URL="/login_dev",
    )
    app.secret_key = "test_secret"
    register_before_request_auth(app)

    @app.route("/whoami")
    def whoami():
        # Stands in for every route that reads session["user"] after the auth hook.
        return jsonify(session["user"])

    return app


@pytest.fixture
def client(app):
    client = app.test_client()
    with client.session_transaction() as sess:
        sess["user"] = dict(FOREIGN_SESSION_USER)
    return client


def test_foreign_session_user_is_replaced_by_local_user(client):
    with patch("mdvtools.auth.authutils.user_cache", {LOCAL_USER["auth_id"]: dict(LOCAL_USER)}):
        rv = client.get("/whoami")

    assert rv.status_code == 200
    assert rv.get_json() == LOCAL_USER
    with client.session_transaction() as sess:
        assert sess["user"] == LOCAL_USER


def test_session_user_missing_from_cache_is_resolved_from_database(client):
    # A user created since the last cache refresh must not be logged out.
    with patch("mdvtools.auth.authutils.user_cache", {}), \
         patch("mdvtools.dbutils.dbmodels.User") as mock_user:
        mock_user.query.filter_by.return_value.first.return_value = SimpleNamespace(**LOCAL_USER)
        rv = client.get("/whoami")

    assert rv.status_code == 200
    assert rv.get_json() == LOCAL_USER
    mock_user.query.filter_by.assert_called_once_with(auth_id=LOCAL_USER["auth_id"])


def test_session_user_unknown_to_this_deployment_is_logged_out(client):
    with patch("mdvtools.auth.authutils.user_cache", {}), \
         patch("mdvtools.dbutils.dbmodels.User") as mock_user:
        mock_user.query.filter_by.return_value.first.return_value = None
        rv = client.get("/whoami")

    assert rv.status_code == 302
    assert rv.headers["Location"].endswith("/login_dev")
    with client.session_transaction() as sess:
        assert "user" not in sess


def test_dummy_auth_session_user_is_left_alone(app):
    # The dummy provider's user has no database row, so resolving it would log developers out.
    app.config["DEFAULT_AUTH_METHOD"] = "dummy"
    dummy_user = DummyAuthProvider(app).get_user()
    client = app.test_client()
    with client.session_transaction() as sess:
        sess["user"] = dummy_user

    with patch("mdvtools.auth.authutils.user_cache", {}), \
         patch("mdvtools.dbutils.dbmodels.User") as mock_user:
        mock_user.query.filter_by.return_value.first.return_value = None
        rv = client.get("/whoami")

    assert rv.status_code == 200
    assert rv.get_json() == dummy_user
