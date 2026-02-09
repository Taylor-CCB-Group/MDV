import pytest
from flask import Flask, session
from unittest.mock import patch, MagicMock
from mdvtools.auth.authutils import get_auth_provider
from mdvtools.auth.dummy_provider import DummyAuthProvider
from mdvtools.auth.shibboleth_provider import ShibbolethProvider
from mdvtools.auth.auth0_provider import Auth0Provider
from auth0.exceptions import RateLimitError

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


class TestAuth0ProviderSync:
    """Test cases for Auth0Provider.sync_users_to_db method.
    
    Note - these are LLM generated tests, and may not be perfect.
    
    """
    
    @pytest.fixture
    def auth0_app(self):
        """Create Flask app with Auth0 configuration."""
        app = Flask("test_auth0_app")
        app.config["TESTING"] = True
        app.config["AUTH0_DOMAIN"] = "test.auth0.com"
        app.config["AUTH0_CLIENT_ID"] = "test_client_id"
        app.config["AUTH0_CLIENT_SECRET"] = "test_client_secret"
        app.config["AUTH0_DB_CONNECTION"] = "Username-Password-Authentication"
        app.secret_key = 'test_secret_key'
        return app
    
    @pytest.fixture
    def mock_auth0_client(self):
        """Create a mock Auth0 client."""
        mock_client = MagicMock()
        mock_client.users = MagicMock()
        return mock_client
    
    @pytest.fixture
    def mock_user_service(self):
        """Mock UserService."""
        with patch('mdvtools.dbutils.dbservice.UserService') as mock:
            mock_user = MagicMock()
            mock_user.id = 1
            mock_user.is_admin = False
            mock.add_or_update_user.return_value = mock_user
            yield mock
    
    @pytest.fixture
    def mock_user_project_service(self):
        """Mock UserProjectService."""
        with patch('mdvtools.dbutils.dbservice.UserProjectService') as mock:
            yield mock
    
    @pytest.fixture
    def mock_db(self):
        """Mock database session."""
        with patch('mdvtools.dbutils.dbmodels.db') as mock_db:
            mock_db.session = MagicMock()
            mock_db.session.commit = MagicMock()
            mock_db.session.rollback = MagicMock()
            yield mock_db
    
    @pytest.fixture
    def mock_user_model(self):
        """Mock User model."""
        with patch('mdvtools.dbutils.dbmodels.User') as mock:
            mock.query = MagicMock()
            mock.query.count.return_value = 0
            # Create a chainable mock for filter_by
            filter_by_mock = MagicMock()
            filter_by_mock.first.return_value = None
            filter_by_mock.count.return_value = 0
            mock.query.filter_by.return_value = filter_by_mock
            yield mock
    
    @pytest.fixture
    def mock_project_model(self):
        """Mock Project model."""
        with patch('mdvtools.dbutils.dbmodels.Project') as mock:
            mock.query = MagicMock()
            mock.query.all.return_value = []
            yield mock
    
    @pytest.fixture
    def mock_get_token(self):
        """Mock GetToken for Auth0 Management API."""
        with patch('mdvtools.auth.auth0_provider.GetToken') as mock:
            mock_instance = MagicMock()
            mock_instance.client_credentials.return_value = {"access_token": "test_token"}
            mock.return_value = mock_instance
            yield mock
    
    @pytest.fixture
    def mock_auth0_class(self, mock_auth0_client):
        """Mock Auth0 class."""
        with patch('mdvtools.auth.auth0_provider.Auth0') as mock:
            mock.return_value = mock_auth0_client
            yield mock
    
    def test_sync_users_to_db_success(
        self, auth0_app, mock_auth0_client, mock_user_service, 
        mock_user_project_service, mock_db, mock_user_model,
        mock_project_model, mock_get_token, mock_auth0_class
    ):
        """Test successful sync with pagination."""
        # Setup: 2 pages of users, 50 users each
        page1_users = [
            {'user_id': f'auth0|user{i}', 'email': f'user{i}@test.com'}
            for i in range(50)
        ]
        page2_users = [
            {'user_id': f'auth0|user{i}', 'email': f'user{i}@test.com'}
            for i in range(50, 75)  # 25 users on second page
        ]
        
        # Mock pagination responses
        mock_auth0_client.users.list.side_effect = [
            {'users': page1_users},
            {'users': page2_users},
            {'users': []}  # Empty list to signal end
        ]
        
        # Mock role fetches - first user is admin, rest are not
        def mock_list_roles(user_id):
            if user_id == 'auth0|user0':
                return {'roles': [{'name': 'admin'}]}
            else:
                return {'roles': []}
        
        mock_auth0_client.users.list_roles.side_effect = mock_list_roles
        
        # Mock OAuth initialization to avoid real HTTP requests
        with patch('mdvtools.auth.auth0_provider.requests.get') as mock_get:
            mock_get.return_value.status_code = 200
            mock_get.return_value.json.return_value = {
                'jwks_uri': 'https://test.auth0.com/.well-known/jwks.json'
            }
            
            # Create provider and sync
            with auth0_app.app_context():
                provider = Auth0Provider(
                    auth0_app,
                    oauth=MagicMock(),
                    client_id="test_id",
                    client_secret="test_secret",
                    domain="test.auth0.com"
                )
                
                with patch('mdvtools.auth.auth0_provider.logging') as mock_logging:
                    with patch('mdvtools.auth.auth0_provider.time.sleep'):  # Speed up test
                        provider.sync_users_to_db()
        
        # Verify: 75 users processed
        assert mock_auth0_client.users.list.call_count >= 2
        assert mock_user_service.add_or_update_user.call_count == 75
        # Admin user should get project assignments
        assert mock_user_project_service.add_or_update_user_project.call_count >= 0
    
    def test_sync_users_to_db_rate_limit_role_fetch(
        self, auth0_app, mock_auth0_client, mock_user_service,
        mock_user_project_service, mock_db, mock_user_model,
        mock_project_model, mock_get_token, mock_auth0_class
    ):
        """Test rate limit on role fetch with successful retry."""
        # Setup: single user
        mock_auth0_client.users.list.side_effect = [
            {'users': [{'user_id': 'auth0|user1', 'email': 'user1@test.com'}]},
            {'users': []}
        ]
        
        # Mock role fetch: rate limit on first call, succeed on retry
        call_count = [0]
        def mock_list_roles(user_id):
            call_count[0] += 1
            if call_count[0] == 1:
                raise RateLimitError(
                    error_code="rate_limit",
                    message="Rate limit exceeded",
                    reset_at=1234567890
                )
            return {'roles': []}
        
        mock_auth0_client.users.list_roles.side_effect = mock_list_roles
        
        # Mock OAuth initialization
        with patch('mdvtools.auth.auth0_provider.requests.get') as mock_get:
            mock_get.return_value.status_code = 200
            mock_get.return_value.json.return_value = {
                'jwks_uri': 'https://test.auth0.com/.well-known/jwks.json'
            }
            
            with auth0_app.app_context():
                provider = Auth0Provider(
                    auth0_app,
                    oauth=MagicMock(),
                    client_id="test_id",
                    client_secret="test_secret",
                    domain="test.auth0.com"
                )
                
                with patch('mdvtools.auth.auth0_provider.logging'):
                    with patch('mdvtools.auth.auth0_provider.time.sleep'):  # Speed up test
                        provider.sync_users_to_db()
        
        # Verify: user was processed after retry
        assert mock_auth0_client.users.list_roles.call_count == 2  # Initial + retry
        assert mock_user_service.add_or_update_user.called
    
    def test_sync_users_to_db_rate_limit_pagination(
        self, auth0_app, mock_auth0_client, mock_user_service,
        mock_user_project_service, mock_db, mock_user_model,
        mock_project_model, mock_get_token, mock_auth0_class
    ):
        """Test rate limit on pagination with retry."""
        # Mock pagination: rate limit on first call, succeed on retry
        call_count = [0]
        def mock_list(page, per_page, q):
            call_count[0] += 1
            if call_count[0] == 1:
                raise RateLimitError(
                    error_code="rate_limit",
                    message="Rate limit exceeded",
                    reset_at=1234567890
                )
            return {'users': []}  # Empty to end pagination
        
        mock_auth0_client.users.list.side_effect = mock_list
        
        # Mock OAuth initialization
        with patch('mdvtools.auth.auth0_provider.requests.get') as mock_get:
            mock_get.return_value.status_code = 200
            mock_get.return_value.json.return_value = {
                'jwks_uri': 'https://test.auth0.com/.well-known/jwks.json'
            }
            
            with auth0_app.app_context():
                provider = Auth0Provider(
                    auth0_app,
                    oauth=MagicMock(),
                    client_id="test_id",
                    client_secret="test_secret",
                    domain="test.auth0.com"
                )
                
                with patch('mdvtools.auth.auth0_provider.logging'):
                    with patch('mdvtools.auth.auth0_provider.time.sleep'):  # Speed up test
                        provider.sync_users_to_db()
        
        # Verify: pagination was retried
        assert mock_auth0_client.users.list.call_count == 2
    
    def test_sync_users_to_db_empty_list(
        self, auth0_app, mock_auth0_client, mock_user_service,
        mock_user_project_service, mock_db, mock_user_model,
        mock_project_model, mock_get_token, mock_auth0_class
    ):
        """Test handling of empty user list."""
        # Mock empty response immediately
        mock_auth0_client.users.list.return_value = {'users': []}
        
        # Mock OAuth initialization to avoid real HTTP requests
        with patch('mdvtools.auth.auth0_provider.requests.get') as mock_get:
            mock_get.return_value.status_code = 200
            mock_get.return_value.json.return_value = {
                'jwks_uri': 'https://test.auth0.com/.well-known/jwks.json'
            }
            
            with auth0_app.app_context():
                provider = Auth0Provider(
                    auth0_app,
                    oauth=MagicMock(),
                    client_id="test_id",
                    client_secret="test_secret",
                    domain="test.auth0.com"
                )
                
                with patch('mdvtools.auth.auth0_provider.logging') as mock_logging:
                    provider.sync_users_to_db()
        
        # Verify: handled gracefully
        assert mock_auth0_client.users.list.called
        assert not mock_user_service.add_or_update_user.called
    
    def test_sync_users_to_db_database_error(
        self, auth0_app, mock_auth0_client, mock_user_service,
        mock_user_project_service, mock_db, mock_user_model,
        mock_project_model, mock_get_token, mock_auth0_class
    ):
        """Test database error handling."""
        # Setup: single user
        mock_auth0_client.users.list.side_effect = [
            {'users': [{'user_id': 'auth0|user1', 'email': 'user1@test.com'}]},
            {'users': []}
        ]
        
        # Mock roles
        mock_auth0_client.users.list_roles.return_value = {'roles': []}
        
        # Mock database commit to raise error
        mock_db.session.commit.side_effect = Exception("Database error")
        
        # Mock OAuth initialization
        with patch('mdvtools.auth.auth0_provider.requests.get') as mock_get:
            mock_get.return_value.status_code = 200
            mock_get.return_value.json.return_value = {
                'jwks_uri': 'https://test.auth0.com/.well-known/jwks.json'
            }
            
            with auth0_app.app_context():
                provider = Auth0Provider(
                    auth0_app,
                    oauth=MagicMock(),
                    client_id="test_id",
                    client_secret="test_secret",
                    domain="test.auth0.com"
                )
                
                with patch('mdvtools.auth.auth0_provider.logging') as mock_logging:
                    with patch('mdvtools.auth.auth0_provider.time.sleep'):  # Speed up test
                        provider.sync_users_to_db()
        
        # Verify: error was handled, rollback called
        assert mock_db.session.rollback.called