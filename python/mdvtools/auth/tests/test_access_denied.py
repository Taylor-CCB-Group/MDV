"""Signed in at Auth0 but not let in: a page that says why, not a login loop.

Redirecting such a request to the login page cannot help - Auth0 still holds a
session, completes the login at once, and sends the person straight back. These
tests pin the reasons, the page, and the startup check that reports missing
first-administrator setup before anyone tries to sign in.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from flask import Flask, render_template, session

from mdvtools.auth.access_denied import (
    DENIED_ACCOUNT_SESSION_KEY,
    access_denied_reason,
    deny_access,
    login_url_with_reason,
)
from mdvtools.auth.auth0_provider import AUTH0_TIMEOUT_SECONDS, Auth0Provider, check_first_admin_setup
from mdvtools.auth.authutils import register_before_request_auth
from mdvtools.auth.register_auth_routes import register_auth_routes


# --- Reasons recorded by validate_user --------------------------------------


@pytest.fixture
def bootstrap_app():
    app = Flask("test_access_denied")
    app.config.update(
        TESTING=True,
        AUTH0_DOMAIN="test.auth0.com",
        AUTH0_CLIENT_ID="test_client_id",
        AUTH0_CLIENT_SECRET="test_client_secret",
        MDV_BOOTSTRAP_ADMIN_EMAIL="admin@example.com",
    )
    app.secret_key = "test_secret_key"
    return app


@pytest.fixture
def provider(bootstrap_app):
    with patch("mdvtools.auth.auth0_provider.requests.get") as mock_get:
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = {"jwks_uri": "https://test.auth0.com/.well-known/jwks.json"}
        return Auth0Provider(
            bootstrap_app,
            oauth=MagicMock(),
            client_id="test_client_id",
            client_secret="test_client_secret",
            domain="test.auth0.com",
        )


def _refusal_after_validate(bootstrap_app, provider, user_info, user_count=0, stale_account=None):
    """Run validate_user for an identity with no local row and no Invitation;
    return the reason it denied and the account it kept in the session for the
    login page."""
    with patch.object(provider, "get_token", return_value="dummy-token"), \
         patch.object(provider, "is_token_valid", return_value=True), \
         patch.object(provider, "get_user", return_value=user_info), \
         patch.object(provider, "_accept_invitation", return_value=None), \
         patch("mdvtools.dbutils.dbmodels.User") as mock_user_class, \
         patch("mdvtools.dbutils.dbmodels.db") as mock_db:
        mock_db.session = MagicMock()
        mock_user_class.query.filter_by.return_value.first.return_value = None
        mock_user_class.query.count.return_value = user_count

        with bootstrap_app.test_request_context("/"):
            if stale_account is not None:
                session[DENIED_ACCOUNT_SESSION_KEY] = stale_account
            result, error = provider.validate_user()
            reason = access_denied_reason()
            account = session.get(DENIED_ACCOUNT_SESSION_KEY)

    assert result is None
    assert error is not None
    return reason, account


def _denial_after_validate(bootstrap_app, provider, user_info, user_count=0):
    return _refusal_after_validate(bootstrap_app, provider, user_info, user_count)[0]


def test_refused_account_is_kept_for_the_login_page(bootstrap_app, provider):
    """The page names the account, so "use a different account" has a referent.
    It travels in the session, not the URL."""
    reason, account = _refusal_after_validate(
        bootstrap_app, provider,
        {"sub": "auth0|stranger", "email": "someone@else.com", "email_verified": True},
        user_count=4,
    )

    assert reason == "no_account"
    assert account == "someone@else.com"


def test_refusal_without_a_known_account_drops_a_stale_one(bootstrap_app, provider):
    reason, account = _refusal_after_validate(
        bootstrap_app, provider, None, stale_account="earlier@example.com"
    )

    assert reason == "unavailable"
    assert account is None


def test_login_page_names_the_refused_account_escaped():
    """The real template: the email is written into the page, HTML-escaped."""
    templates = Path(__file__).resolve().parents[2] / "templates"
    app = Flask("test_login_template", template_folder=str(templates))
    app.secret_key = "test_secret"
    app.jinja_env.globals["entry_assets"] = lambda _name: {"css": [], "js": "login.js"}

    with app.test_request_context("/login?denied=no_account"):
        session[DENIED_ACCOUNT_SESSION_KEY] = 'jane"><script>@example.com'
        html = render_template("login.html")

    assert 'content="jane&#34;&gt;&lt;script&gt;@example.com"' in html
    assert "<script>@example.com" not in html


def test_unverified_bootstrap_address_is_told_to_verify(bootstrap_app, provider):
    reason = _denial_after_validate(
        bootstrap_app, provider,
        {"sub": "auth0|new", "email": "admin@example.com", "email_verified": False},
    )

    assert reason == "email_unverified"


def test_other_unverified_address_is_not_told_to_verify(bootstrap_app, provider):
    """Verifying cannot help someone who is not the bootstrap address, so the
    address is compared first and they never see that advice."""
    reason = _denial_after_validate(
        bootstrap_app, provider,
        {"sub": "auth0|stranger", "email": "someone@else.com", "email_verified": False},
    )

    assert reason == "not_set_up"


def test_missing_bootstrap_setting_is_logged_and_explained(bootstrap_app, provider, caplog):
    """Used to return silently: no log line and a login loop."""
    bootstrap_app.config["MDV_BOOTSTRAP_ADMIN_EMAIL"] = None

    reason = _denial_after_validate(
        bootstrap_app, provider,
        {"sub": "auth0|new", "email": "admin@example.com", "email_verified": True},
    )

    assert reason == "not_set_up"
    assert "MDV_BOOTSTRAP_ADMIN_EMAIL is not set" in caplog.text


def test_unknown_user_on_established_deployment_is_logged_and_explained(bootstrap_app, provider, caplog):
    """The other refusal that used to be silent."""
    with caplog.at_level("INFO"):
        reason = _denial_after_validate(
            bootstrap_app, provider,
            {"sub": "auth0|stranger", "email": "someone@else.com", "email_verified": True},
            user_count=4,
        )

    assert reason == "no_account"
    assert "auth0|stranger" in caplog.text
    assert "has no account in this deployment" in caplog.text


def test_failed_userinfo_is_explained(bootstrap_app, provider):
    """get_user returns None when Auth0's userinfo fails or times out; without a
    reason this was a plain redirect to sign-in, which loops while Auth0 is down."""
    reason = _denial_after_validate(bootstrap_app, provider, None)

    assert reason == "unavailable"


def test_failed_bootstrap_is_explained(bootstrap_app, provider):
    with patch.object(provider, "_create_first_admin", side_effect=RuntimeError("database is locked")):
        reason = _denial_after_validate(
            bootstrap_app, provider,
            {"sub": "auth0|new", "email": "admin@example.com", "email_verified": True},
        )

    assert reason == "setup_failed"


# --- Where a denied request is sent ------------------------------------------


@pytest.fixture
def hooked_app():
    app = Flask("test_access_denied_hook")
    app.config.update(
        TESTING=True,
        ENABLE_AUTH=True,
        DEFAULT_AUTH_METHOD="auth0",
        LOGIN_REDIRECT_URL="/login_dev",
    )
    app.secret_key = "test_secret"
    register_before_request_auth(app)

    @app.route("/")
    def index():
        return "the app"

    @app.route("/logout")
    def logout():
        return "signed out"

    return app


def _provider_that(reason=None):
    """A provider whose validation fails, denying with `reason` if one is given."""
    provider = MagicMock()

    def validate_user():
        if reason is not None:
            deny_access(reason)
        return None, ("error", 404)

    provider.validate_user.side_effect = validate_user
    return provider


def test_denied_browser_request_goes_to_login_with_the_reason(hooked_app):
    """The login page shows the reason instead of a sign-in button, which is what
    breaks the loop: signing in again is what sent people straight back."""
    with patch("mdvtools.auth.authutils.get_auth_provider", return_value=_provider_that("email_unverified")):
        rv = hooked_app.test_client().get("/", headers={"Accept": "text/html,*/*;q=0.8"})

    assert rv.status_code == 302
    assert rv.headers["Location"] == "/login_dev?denied=email_unverified"


def test_reason_is_added_to_a_login_url_that_already_has_a_query():
    assert login_url_with_reason("/login_dev", "no_account") == "/login_dev?denied=no_account"
    assert (
        login_url_with_reason("https://example.org/lucas-test/login_dev?callback=1", "no_account")
        == "https://example.org/lucas-test/login_dev?callback=1&denied=no_account"
    )


def test_denied_api_request_is_redirected_too(hooked_app):
    """The same as a request with no session: one behaviour, not two."""
    with patch("mdvtools.auth.authutils.get_auth_provider", return_value=_provider_that("no_account")):
        rv = hooked_app.test_client().get("/", headers={"Accept": "application/json"})

    assert rv.status_code == 302
    assert rv.headers["Location"] == "/login_dev?denied=no_account"


def test_not_signed_in_still_redirects_to_plain_login(hooked_app):
    """Only a denial carries a reason; an ordinary missing session is unchanged."""
    with patch("mdvtools.auth.authutils.get_auth_provider", return_value=_provider_that()):
        rv = hooked_app.test_client().get("/", headers={"Accept": "text/html"})

    assert rv.status_code == 302
    assert rv.headers["Location"] == "/login_dev"


def test_logout_is_reachable_without_a_session(hooked_app):
    """Someone turned away has to be able to end their Auth0 session to switch
    account - before this, /logout itself redirected to login."""
    with patch("mdvtools.auth.authutils.get_auth_provider",
               return_value=_provider_that("no_account")) as mock_provider:
        rv = hooked_app.test_client().get("/logout")

    assert rv.status_code == 200
    assert rv.get_data(as_text=True) == "signed out"
    mock_provider.assert_not_called()


def _app_with_login_route(tmp_path):
    """LOGIN_REDIRECT_URL=/login, as deploy.sh sets it for Auth0, with the real routes.

    Build it while get_auth_provider is patched: the routes bind it at registration.
    """
    (tmp_path / "login.html").write_text("the login page")
    app = Flask("test_access_denied_login_route", template_folder=str(tmp_path))
    app.config.update(
        TESTING=True,
        ENABLE_AUTH=True,
        DEFAULT_AUTH_METHOD="auth0",
        LOGIN_REDIRECT_URL="/login",
    )
    app.secret_key = "test_secret"
    register_before_request_auth(app)
    register_auth_routes(app)

    @app.route("/")
    def index():
        return "the app"

    return app


def test_denial_at_login_route_shows_the_page_instead_of_signing_in_again(tmp_path):
    """/login used to ignore ?denied= and start Auth0, which signed the same person
    straight back in - the loop, on every deployment set up by deploy.sh."""
    provider = _provider_that("no_account")
    with patch("mdvtools.auth.authutils.get_auth_provider", return_value=provider):
        client = _app_with_login_route(tmp_path).test_client()
        denied = client.get("/", headers={"Accept": "text/html"})
        page = client.get(denied.headers["Location"])

    assert denied.headers["Location"] == "/login?denied=no_account"
    assert page.status_code == 200
    assert page.get_data(as_text=True) == "the login page"
    provider.login.assert_not_called()


def test_logout_without_a_session_still_ends_the_auth0_session(tmp_path):
    """With no local session, /logout used to clear nothing and stop, leaving Auth0
    signed in - so the refusal page's Sign out could not switch account."""
    provider = _provider_that()
    provider.logout.return_value = "to Auth0 logout"
    with patch("mdvtools.auth.authutils.get_auth_provider", return_value=provider):
        rv = _app_with_login_route(tmp_path).test_client().get("/logout")

    assert rv.get_data(as_text=True) == "to Auth0 logout"
    provider.logout.assert_called_once()


def test_logout_without_a_session_signs_out_of_the_default_provider_for_other_auth_methods(tmp_path):
    provider = _provider_that()
    provider.logout.return_value = "to the provider's logout"
    with patch("mdvtools.auth.authutils.get_auth_provider", return_value=provider):
        app = _app_with_login_route(tmp_path)
        app.config["DEFAULT_AUTH_METHOD"] = "shibboleth"
        rv = app.test_client().get("/logout")

    assert rv.get_data(as_text=True) == "to the provider's logout"
    provider.logout.assert_called_once()


def test_login_route_without_a_reason_still_starts_sign_in(tmp_path):
    provider = _provider_that()
    provider.login.return_value = "to Auth0"
    with patch("mdvtools.auth.authutils.get_auth_provider", return_value=provider):
        rv = _app_with_login_route(tmp_path).test_client().get("/login")

    assert rv.get_data(as_text=True) == "to Auth0"
    provider.login.assert_called_once()


# --- Startup check ------------------------------------------------------------
# Admin rights live only in MDV (ADR 0019), so the first Admin needs only
# MDV_BOOTSTRAP_ADMIN_EMAIL: the check does not call Auth0.


def _check(app, user_count=0):
    with patch("mdvtools.dbutils.dbmodels.User") as mock_user_class:
        mock_user_class.query.count.return_value = user_count
        return check_first_admin_setup(app)


@pytest.fixture
def startup_app(bootstrap_app):
    bootstrap_app.config["DEFAULT_AUTH_METHOD"] = "auth0"
    return bootstrap_app


def test_startup_check_passes_a_complete_setup(startup_app):
    assert _check(startup_app) == []


def test_startup_check_skips_a_deployment_with_users(startup_app):
    startup_app.config["MDV_BOOTSTRAP_ADMIN_EMAIL"] = ""

    assert _check(startup_app, user_count=3) == []


def test_startup_check_skips_other_auth_methods(startup_app):
    startup_app.config["DEFAULT_AUTH_METHOD"] = "shibboleth"
    startup_app.config["MDV_BOOTSTRAP_ADMIN_EMAIL"] = ""

    assert _check(startup_app) == []


def test_startup_check_reports_a_missing_bootstrap_address(startup_app):
    startup_app.config["MDV_BOOTSTRAP_ADMIN_EMAIL"] = ""

    problems = _check(startup_app)

    assert len(problems) == 1
    assert "MDV_BOOTSTRAP_ADMIN_EMAIL is not set" in problems[0]


# --- Timeouts -----------------------------------------------------------------


def test_user_info_request_has_a_timeout(provider):
    """It had none, so a slow Auth0 held the request - 21 seconds on one login,
    and indefinitely if Auth0 never answered."""
    with patch("mdvtools.auth.auth0_provider.requests.get") as mock_get:
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = {"sub": "auth0|x", "email": "x@example.com"}
        provider.get_user({"access_token": "token"})

    assert mock_get.call_args.kwargs["timeout"] == AUTH0_TIMEOUT_SECONDS


def test_oauth_client_is_registered_with_a_timeout(provider):
    """Covers Authlib's own calls, including the token exchange after login."""
    client_kwargs = provider.oauth.register.call_args.kwargs["client_kwargs"]
    assert client_kwargs["default_timeout"] == AUTH0_TIMEOUT_SECONDS
