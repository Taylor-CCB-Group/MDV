"""Signed in at Auth0 but not let in: a page that says why, not a login loop.

Redirecting such a request to the login page cannot help - Auth0 still holds a
session, completes the login at once, and sends the person straight back. These
tests pin the reasons, the page, and the startup check that reports missing
first-administrator setup before anyone tries to sign in.
"""

from unittest.mock import MagicMock, patch

import pytest
from flask import Flask

from mdvtools.auth.access_denied import access_denied_reason, deny_access, login_url_with_reason
from mdvtools.auth.auth0_provider import AUTH0_TIMEOUT_SECONDS, Auth0Provider, check_first_admin_setup
from mdvtools.auth.authutils import register_before_request_auth


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


def _denial_after_validate(bootstrap_app, provider, user_info, user_count=0):
    """Run validate_user for an identity with no local row; return what it denied."""
    with patch.object(provider, "get_token", return_value="dummy-token"), \
         patch.object(provider, "is_token_valid", return_value=True), \
         patch.object(provider, "get_user", return_value=user_info), \
         patch("mdvtools.dbutils.dbmodels.User") as mock_user_class, \
         patch("mdvtools.dbutils.dbmodels.db") as mock_db:
        mock_db.session = MagicMock()
        mock_user_class.query.filter_by.return_value.first.return_value = None
        mock_user_class.query.count.return_value = user_count

        with bootstrap_app.test_request_context("/"):
            result, error = provider.validate_user()
            reason = access_denied_reason()

    assert result is None
    assert error is not None
    return reason


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


def test_failed_bootstrap_is_explained(bootstrap_app, provider):
    with patch("mdvtools.auth.auth0_provider.GetToken") as mock_get_token, \
         patch("mdvtools.auth.auth0_provider.Auth0") as mock_auth0_class:
        mock_get_token.return_value.client_credentials.return_value = {"access_token": "mgmt-token"}
        mock_auth0_class.return_value.roles.list.return_value = {"roles": []}

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


# --- Startup check ------------------------------------------------------------


def _check(app, user_count=0, scope="read:roles update:users read:users", roles=None, token_error=None):
    with patch("mdvtools.dbutils.dbmodels.User") as mock_user_class, \
         patch("mdvtools.auth.auth0_provider.GetToken") as mock_get_token, \
         patch("mdvtools.auth.auth0_provider.Auth0") as mock_auth0_class:
        mock_user_class.query.count.return_value = user_count
        token_response = {"access_token": "mgmt-token"}
        if scope is not None:
            token_response["scope"] = scope
        if token_error is not None:
            mock_get_token.return_value.client_credentials.side_effect = token_error
        else:
            mock_get_token.return_value.client_credentials.return_value = token_response
        mock_auth0 = mock_auth0_class.return_value
        mock_auth0.roles.list.return_value = {
            "roles": roles if roles is not None else [{"id": "rol_admin", "name": "admin"}]
        }

        problems = check_first_admin_setup(app)

    return problems, mock_get_token, mock_auth0


@pytest.fixture
def startup_app(bootstrap_app):
    bootstrap_app.config["DEFAULT_AUTH_METHOD"] = "auth0"
    return bootstrap_app


def test_startup_check_passes_a_complete_setup(startup_app):
    problems, _get_token, _auth0 = _check(startup_app)

    assert problems == []


def test_startup_check_skips_a_deployment_with_users(startup_app):
    problems, mock_get_token, _auth0 = _check(startup_app, user_count=3)

    assert problems == []
    mock_get_token.assert_not_called()


def test_startup_check_skips_other_auth_methods(startup_app):
    startup_app.config["DEFAULT_AUTH_METHOD"] = "shibboleth"

    problems, mock_get_token, _auth0 = _check(startup_app)

    assert problems == []
    mock_get_token.assert_not_called()


def test_startup_check_reports_a_missing_bootstrap_address(startup_app):
    startup_app.config["MDV_BOOTSTRAP_ADMIN_EMAIL"] = ""

    problems, _get_token, _auth0 = _check(startup_app)

    assert len(problems) == 1
    assert "MDV_BOOTSTRAP_ADMIN_EMAIL is not set" in problems[0]


def test_startup_check_reports_a_blocked_management_grant(startup_app):
    """The Token Endpoint Authentication Method = None case, which looked like a
    session problem on a real server."""
    problems, _get_token, _auth0 = _check(
        startup_app, token_error=Exception("Grant type 'client_credentials' not allowed for the client.")
    )

    assert len(problems) == 1
    assert "Client Secret (Post)" in problems[0]
    assert "not allowed for the client" in problems[0]


def test_startup_check_reports_missing_scopes(startup_app):
    problems, _get_token, mock_auth0 = _check(startup_app, scope="read:roles read:users")

    assert len(problems) == 1
    assert "update:users" in problems[0]
    # read:roles is granted, so the role is still checked.
    mock_auth0.roles.list.assert_called_once()


def test_startup_check_does_not_read_roles_without_permission(startup_app):
    problems, _get_token, mock_auth0 = _check(startup_app, scope="update:users")

    assert len(problems) == 1
    assert "read:roles" in problems[0]
    mock_auth0.roles.list.assert_not_called()


def test_startup_check_reports_a_tenant_with_no_admin_role(startup_app):
    problems, _get_token, _auth0 = _check(startup_app, roles=[{"id": "rol_viewer", "name": "viewer"}])

    assert len(problems) == 1
    assert "no role named exactly 'admin'" in problems[0]


def test_startup_check_trusts_a_token_that_does_not_report_scopes(startup_app):
    """No scope field means unknown, not none - reporting every scope missing
    would send an operator chasing a grant that is already there."""
    problems, _get_token, _auth0 = _check(startup_app, scope=None)

    assert problems == []


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


def test_management_api_calls_have_a_timeout(startup_app):
    """The bootstrap's role assignment failed at the library's default of 5
    seconds on a network where every TLS handshake took 5.06."""
    _problems, mock_get_token, _auth0 = _check(startup_app)

    assert mock_get_token.call_args.kwargs["timeout"] == AUTH0_TIMEOUT_SECONDS


def test_management_client_is_built_with_a_timeout(startup_app):
    with patch("mdvtools.dbutils.dbmodels.User") as mock_user_class, \
         patch("mdvtools.auth.auth0_provider.GetToken") as mock_get_token, \
         patch("mdvtools.auth.auth0_provider.Auth0") as mock_auth0_class:
        mock_user_class.query.count.return_value = 0
        mock_get_token.return_value.client_credentials.return_value = {"access_token": "mgmt-token"}
        mock_auth0_class.return_value.roles.list.return_value = {"roles": [{"id": "rol_admin", "name": "admin"}]}

        check_first_admin_setup(startup_app)

    assert mock_auth0_class.call_args.kwargs["rest_options"].timeout == AUTH0_TIMEOUT_SECONDS
