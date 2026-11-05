"""What a person signed in to Auth0 sees when this deployment gives them no access.

ADRs 0019 and 0020: someone with no User, or a Deactivated user, gets a page
with the email they used, the site contact and a way to sign out. Before this,
they got JSON or a redirect that sent them back through Auth0 in a loop.

Requests go through the real auth hook and Auth0Provider against an in-memory
database. Only the calls that would reach Auth0 are stubbed.
"""

import re
from pathlib import Path
from unittest.mock import MagicMock, patch
from urllib.parse import urljoin, urlsplit

import pytest
from flask import Flask, jsonify, session

import mdvtools
from mdvtools.auth.auth0_provider import Auth0Provider
from mdvtools.auth.authutils import register_before_request_auth
from mdvtools.dbutils.dbmodels import User, db

CONTACT = "MDV support, mdv-support@example.org"


@pytest.fixture()
def app():
    app = Flask(__name__, template_folder=str(Path(mdvtools.__file__).parent / "templates"))
    app.secret_key = "test"
    app.config.update(
        ENABLE_AUTH=True,
        DEFAULT_AUTH_METHOD="auth0",
        LOGIN_REDIRECT_URL="/login_dev",
        MDV_ACCESS_CONTACT=CONTACT,
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    db.init_app(app)
    register_before_request_auth(app)

    @app.route("/whoami")
    def whoami():
        # Stands in for every route behind the auth hook.
        return jsonify(session["user"])

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


def signed_in_to_auth0(provider, *, sub, email, connection="Username-Password-Authentication"):
    """Make Auth0 report a valid sign-in for this account."""
    sign_in = {
        "sub": sub,
        "email": email,
        "email_verified": True,
        "connection": connection,
        "connection_strategy": "auth0",
    }
    return patch.multiple(
        provider,
        get_token=MagicMock(return_value="token"),
        is_token_valid=MagicMock(return_value=True),
        get_user=MagicMock(return_value=sign_in),
    )


def get_as(app, provider, path, *, sub, email, **sign_in):
    with signed_in_to_auth0(provider, sub=sub, email=email, **sign_in), patch(
        "mdvtools.auth.authutils.get_auth_provider", return_value=provider
    ):
        return app.test_client().get(path)


def test_person_with_no_user_sees_the_no_access_page(app, provider):
    response = get_as(app, provider, "/whoami", sub="auth0|stranger", email="stranger@example.org")

    assert response.status_code == 403
    page = response.get_data(as_text=True)
    assert "stranger@example.org" in page
    assert CONTACT in page
    assert sign_out_target(page, "/whoami") == "/logout"


def test_no_access_page_shows_the_sign_in_method_used(app, provider):
    # ADR 0019: someone who used the wrong Sign-in method can see that they did.
    response = get_as(
        app, provider, "/whoami", sub="waad|stranger", email="stranger@acme.example", connection="acme-entra"
    )

    assert response.status_code == 403
    assert "acme-entra" in response.get_data(as_text=True)


def sign_out_target(page, requested_path):
    """The path the page's sign-out link leads to, resolved as a browser would."""
    href = re.search(r'<a href="([^"]*)">Sign out', page)
    assert href, "the page has no sign-out link"
    return urlsplit(urljoin(f"http://mdv.test{requested_path}", href.group(1))).path


def test_deactivated_user_sees_the_no_access_page_and_stays_deactivated(app, provider):
    db.session.add(
        User(email="former@example.org", auth_id="auth0|former", is_active=False, password="")
    )
    db.session.commit()

    response = get_as(app, provider, "/whoami", sub="auth0|former", email="former@example.org")

    assert response.status_code == 403
    page = response.get_data(as_text=True)
    assert "former@example.org" in page
    assert CONTACT in page
    assert sign_out_target(page, "/whoami") == "/logout"
    db.session.expire_all()
    assert User.query.filter_by(auth_id="auth0|former").one().is_active is False


def test_active_user_goes_through(app, provider):
    db.session.add(
        User(email="member@example.org", auth_id="auth0|member", is_active=True, password="")
    )
    db.session.commit()

    response = get_as(app, provider, "/whoami", sub="auth0|member", email="member@example.org")

    assert response.status_code == 200
    assert response.get_json()["email"] == "member@example.org"


def test_open_session_is_refused_once_the_user_is_deactivated(app, provider):
    # ADR 0020: deactivating someone ends access they already have, not only
    # their next sign-in. MDV-Admin writes is_active straight to the database.
    db.session.add(
        User(email="leaver@example.org", auth_id="auth0|leaver", is_active=True, password="")
    )
    db.session.commit()
    client = app.test_client()

    with signed_in_to_auth0(provider, sub="auth0|leaver", email="leaver@example.org"), patch(
        "mdvtools.auth.authutils.get_auth_provider", return_value=provider
    ):
        assert client.get("/whoami").status_code == 200

        User.query.filter_by(auth_id="auth0|leaver").one().is_active = False
        db.session.commit()

        response = client.get("/whoami")

    assert response.status_code == 403
    page = response.get_data(as_text=True)
    assert "leaver@example.org" in page
    assert CONTACT in page
    assert sign_out_target(page, "/whoami") == "/logout"
