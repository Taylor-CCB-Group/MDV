"""What a person signed in to Auth0 sees when this deployment gives them no access.

ADRs 0019 and 0020: someone with no User, or a Deactivated user, is sent to the
login page with the no_account reason. The page names the email and Sign-in
method they used and the site contact, and offers a way to sign out. A
Deactivated user gets the same page as someone with no User, so it does not show
that their User still exists.

Requests go through the real auth hook, auth routes, login template and
Auth0Provider against an in-memory database. Only the calls that would reach
Auth0 are stubbed.
"""

import re
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from flask import Flask, jsonify, session

import mdvtools
from mdvtools.auth.auth0_provider import Auth0Provider
from mdvtools.auth.authutils import register_before_request_auth
from mdvtools.auth.register_auth_routes import register_auth_routes
from mdvtools.dbutils.dbmodels import User, db

CONTACT = "MDV support, mdv-support@example.org"
PASSWORD = "Username-Password-Authentication"
REFUSED = "/login?denied=no_account"
AUTH0_LOGOUT = "https://test.auth0.com/v2/logout"


@pytest.fixture()
def app():
    app = Flask(__name__, template_folder=str(Path(mdvtools.__file__).parent / "templates"))
    app.secret_key = "test"
    app.config.update(
        ENABLE_AUTH=True,
        DEFAULT_AUTH_METHOD="auth0",
        LOGIN_REDIRECT_URL="/login",
        AUTH0_DOMAIN="test.auth0.com",
        AUTH0_CLIENT_ID="client",
        MDV_ACCESS_CONTACT=CONTACT,
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    # The login page's script and styles do not matter here, only its meta tags.
    app.jinja_env.globals["entry_assets"] = lambda _name: {"css": [], "js": "login.js"}
    db.init_app(app)
    register_before_request_auth(app)

    @app.route("/whoami")
    def whoami():
        # Stands in for every route behind the auth hook.
        return jsonify(session["user"])

    with app.app_context():
        db.create_all()
        # An established deployment. With no users at all, a refused sign-in is
        # told the deployment is not set up yet.
        add_user("admin@example.org", "auth0|admin", is_active=True)
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
        provider = Auth0Provider(
            app,
            oauth=MagicMock(),
            client_id="client",
            client_secret="secret",
            domain="test.auth0.com",
        )
    # The auth routes look up the provider when they are registered.
    with patch("mdvtools.auth.authutils.get_auth_provider", return_value=provider):
        register_auth_routes(app)
    return provider


def signed_in_to_auth0(provider, *, sub, email, connection=PASSWORD):
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


def browser(app):
    """A browser that started at /login, which records the sign-in method."""
    client = app.test_client()
    with client.session_transaction() as browser_session:
        browser_session["auth_method"] = "auth0"
    return client


def refusal_page(client, response):
    """Follow a refusal to the login page and return the details it shows."""
    assert response.status_code == 302, response.get_data(as_text=True)
    assert response.headers["Location"] == REFUSED
    page = client.get(REFUSED)
    assert page.status_code == 200
    html = page.get_data(as_text=True)
    return {
        name: re.search(rf'<meta name="mdv-{name}" content="([^"]*)">', html).group(1)
        for name in ("denied-account", "denied-sign-in-method", "access-contact")
    }


def add_user(email, auth_id, *, is_active):
    db.session.add(User(email=email, auth_id=auth_id, is_active=is_active, password=""))
    db.session.commit()


def test_person_with_no_user_is_told_they_have_no_access(app, provider):
    client = browser(app)
    with signed_in_to_auth0(provider, sub="auth0|stranger", email="stranger@example.org"), patch(
        "mdvtools.auth.authutils.get_auth_provider", return_value=provider
    ):
        shown = refusal_page(client, client.get("/whoami"))

    assert shown == {
        "denied-account": "stranger@example.org",
        "denied-sign-in-method": PASSWORD,
        "access-contact": CONTACT,
    }


def test_refusal_shows_the_sign_in_method_used(app, provider):
    # ADR 0019: someone who used the wrong Sign-in method can see that they did.
    client = browser(app)
    with signed_in_to_auth0(
        provider, sub="waad|stranger", email="stranger@acme.example", connection="acme-entra"
    ), patch("mdvtools.auth.authutils.get_auth_provider", return_value=provider):
        shown = refusal_page(client, client.get("/whoami"))

    assert shown["denied-sign-in-method"] == "acme-entra"


def test_deactivated_user_is_told_they_have_no_access_and_stays_deactivated(app, provider):
    add_user("former@example.org", "auth0|former", is_active=False)
    client = browser(app)

    with signed_in_to_auth0(provider, sub="auth0|former", email="former@example.org"), patch(
        "mdvtools.auth.authutils.get_auth_provider", return_value=provider
    ):
        shown = refusal_page(client, client.get("/whoami"))

    assert shown["denied-account"] == "former@example.org"
    db.session.expire_all()
    assert User.query.filter_by(auth_id="auth0|former").one().is_active is False


def test_active_user_goes_through(app, provider):
    add_user("member@example.org", "auth0|member", is_active=True)

    with signed_in_to_auth0(provider, sub="auth0|member", email="member@example.org"), patch(
        "mdvtools.auth.authutils.get_auth_provider", return_value=provider
    ):
        response = browser(app).get("/whoami")

    assert response.status_code == 200
    assert response.get_json()["email"] == "member@example.org"


def test_open_session_is_refused_once_the_user_is_deactivated(app, provider):
    # ADR 0020: deactivating someone ends access they already have, not only
    # their next sign-in. MDV-Admin writes is_active straight to the database.
    add_user("leaver@example.org", "auth0|leaver", is_active=True)
    client = browser(app)

    with signed_in_to_auth0(provider, sub="auth0|leaver", email="leaver@example.org"), patch(
        "mdvtools.auth.authutils.get_auth_provider", return_value=provider
    ):
        assert client.get("/whoami").status_code == 200

        User.query.filter_by(auth_id="auth0|leaver").one().is_active = False
        db.session.commit()

        shown = refusal_page(client, client.get("/whoami"))

    assert shown["denied-account"] == "leaver@example.org"
    assert shown["access-contact"] == CONTACT


def test_person_with_no_user_can_sign_out_of_auth0(app, provider):
    client = browser(app)

    with signed_in_to_auth0(provider, sub="auth0|stranger", email="stranger@example.org"), patch(
        "mdvtools.auth.authutils.get_auth_provider", return_value=provider
    ):
        refusal_page(client, client.get("/whoami"))
        response = client.get("/logout")

    assert response.status_code == 302
    assert response.headers["Location"].startswith(AUTH0_LOGOUT)


def test_deactivated_user_can_sign_out_of_auth0_after_their_session_is_cleared(app, provider):
    add_user("leaver@example.org", "auth0|leaver", is_active=True)
    client = browser(app)

    with signed_in_to_auth0(provider, sub="auth0|leaver", email="leaver@example.org"), patch(
        "mdvtools.auth.authutils.get_auth_provider", return_value=provider
    ):
        assert client.get("/whoami").status_code == 200

        User.query.filter_by(auth_id="auth0|leaver").one().is_active = False
        db.session.commit()

        refusal_page(client, client.get("/whoami"))
        response = client.get("/logout")

    assert response.status_code == 302
    assert response.headers["Location"].startswith(AUTH0_LOGOUT)
