import pytest
from flask import Flask, session
from unittest.mock import patch, MagicMock
from mdvtools.auth.authutils import get_auth_provider
from mdvtools.auth.dummy_provider import DummyAuthProvider
from mdvtools.auth.shibboleth_provider import ShibbolethProvider

@pytest.fixture
def app():
    """Create and configure a new app instance for each test."""
    app = Flask("test_app")
    app.config["TESTING"] = True
    app.secret_key = 'supersecretkey'
    yield app

def test_get_auth_provider_fails_without_config(app):
    """
    Tests that get_auth_provider raises ValueError if DEFAULT_AUTH_METHOD is not set.
    """
    with app.app_context():
        with pytest.raises(ValueError, match="must be explicitly configured"):
            get_auth_provider()

def test_get_auth_provider_auth0_from_config(app):
    """
    Tests that get_auth_provider attempts to create the Auth0Provider
    when configured, without needing a live connection.
    """
    with patch('mdvtools.auth.auth0_provider.Auth0Provider') as mock_auth0_provider, \
         patch('mdvtools.dbutils.mdv_server_app.oauth', MagicMock()):

        app.config["DEFAULT_AUTH_METHOD"] = "auth0"
        app.config["AUTH0_CLIENT_ID"] = "dummy_id"
        app.config["AUTH0_CLIENT_SECRET"] = "dummy_secret"
        app.config["AUTH0_DOMAIN"] = "dummy.domain"

        with app.app_context():
            provider = get_auth_provider()
            mock_auth0_provider.assert_called_once()
            assert provider == mock_auth0_provider.return_value

def test_get_auth_provider_dummy_override(app):
    """
    Tests that if the default method is 'dummy', it forces a dummy provider
    even if the session requests 'shibboleth'. This is the developer override.
    """
    app.config["DEFAULT_AUTH_METHOD"] = "dummy"
    with app.test_request_context('/'):
        session['auth_method'] = 'shibboleth'
        provider = get_auth_provider()
        assert isinstance(provider, DummyAuthProvider)

def test_get_auth_provider_session_overrides_non_dummy_default(app):
    """
    Tests that the session 'auth_method' overrides the default config
    when the default is NOT 'dummy'.
    """
    app.config["DEFAULT_AUTH_METHOD"] = "auth0"
    # The Auth0Provider requires this config, even if we expect Shibboleth
    app.config["AUTH0_CLIENT_ID"] = "dummy_id"
    app.config["AUTH0_CLIENT_SECRET"] = "dummy_secret"
    app.config["AUTH0_DOMAIN"] = "dummy.domain"

    with app.test_request_context('/'):
        session['auth_method'] = 'shibboleth'
        provider = get_auth_provider()
        assert isinstance(provider, ShibbolethProvider) 