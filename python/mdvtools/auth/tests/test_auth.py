import pytest
from flask import Flask, session
from unittest.mock import patch, MagicMock, ANY
from sqlalchemy.exc import IntegrityError
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
        # `>= 0` is vacuous and doesn't test anything
        # assert mock_user_project_service.add_or_update_user_project.call_count >= 0
    
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


class TestValidateUserBootstrapAndActivation:
    """Tests for Auth0Provider.validate_user's first-admin bootstrap and
    pending-user activation behavior."""

    @pytest.fixture
    def bootstrap_app(self):
        app = Flask("test_bootstrap_app")
        app.config["TESTING"] = True
        app.config["AUTH0_DOMAIN"] = "test.auth0.com"
        app.config["AUTH0_CLIENT_ID"] = "test_client_id"
        app.config["AUTH0_CLIENT_SECRET"] = "test_client_secret"
        app.config["MDV_BOOTSTRAP_ADMIN_EMAIL"] = "admin@example.com"
        app.secret_key = "test_secret_key"
        return app

    @pytest.fixture
    def provider(self, bootstrap_app):
        with patch('mdvtools.auth.auth0_provider.requests.get') as mock_get:
            mock_get.return_value.status_code = 200
            mock_get.return_value.json.return_value = {
                'jwks_uri': 'https://test.auth0.com/.well-known/jwks.json'
            }
            return Auth0Provider(
                bootstrap_app,
                oauth=MagicMock(),
                client_id="test_client_id",
                client_secret="test_client_secret",
                domain="test.auth0.com",
            )

    def _run_validate_user(self, bootstrap_app, provider, user_info, existing_user_side_effect):
        """Drive validate_user with get_token/is_token_valid/get_user stubbed, and
        User.query wired so `.filter_by(...).first()` returns each item of
        `existing_user_side_effect` in order across however many times it's called."""
        with patch.object(provider, 'get_token', return_value='dummy-token'), \
             patch.object(provider, 'is_token_valid', return_value=True), \
             patch.object(provider, 'get_user', return_value=user_info), \
             patch('mdvtools.dbutils.dbmodels.User') as mock_user_class, \
             patch('mdvtools.dbutils.dbmodels.db') as mock_db:

            mock_db.session = MagicMock()
            filter_by_mock = MagicMock()
            filter_by_mock.first.side_effect = existing_user_side_effect
            mock_user_class.query.filter_by.return_value = filter_by_mock
            mock_user_class.query.count.return_value = 0

            new_user = MagicMock()
            new_user.id = 99
            new_user.auth_id = user_info.get("sub")
            new_user.email = user_info.get("email")
            new_user.is_admin = True
            mock_user_class.return_value = new_user

            with bootstrap_app.test_request_context('/'):
                result, error = provider.validate_user()

        return result, error, mock_user_class, mock_db, new_user

    def test_bootstrap_creates_first_admin_and_assigns_role(self, bootstrap_app, provider):
        user_info = {"sub": "auth0|new-admin", "email": "admin@example.com", "email_verified": True}

        with patch('mdvtools.auth.auth0_provider.GetToken') as mock_get_token, \
             patch('mdvtools.auth.auth0_provider.Auth0') as mock_auth0_class:
            mock_get_token.return_value.client_credentials.return_value = {"access_token": "mgmt-token"}
            mock_auth0 = MagicMock()
            mock_auth0.roles.list.return_value = {"roles": [{"id": "role_admin", "name": "admin"}]}
            mock_auth0_class.return_value = mock_auth0

            result, error, mock_user_class, mock_db, new_user = self._run_validate_user(
                bootstrap_app, provider, user_info, existing_user_side_effect=[None]
            )

        assert error is None
        assert result == {"id": 99, "auth_id": "auth0|new-admin", "email": "admin@example.com", "is_admin": True}
        mock_user_class.assert_called_once_with(
            email="admin@example.com",
            auth_id="auth0|new-admin",
            confirmed_at=ANY,
            is_active=True,
            administrator=True,
            is_admin=True,
            password="",
        )
        mock_db.session.add.assert_called_once_with(new_user)
        mock_auth0.users.add_roles.assert_called_once_with("auth0|new-admin", ["role_admin"])

    def test_bootstrap_rejects_nonmatching_email(self, bootstrap_app, provider):
        user_info = {"sub": "auth0|someone-else", "email": "someone@else.com", "email_verified": True}

        result, error, mock_user_class, mock_db, _new_user = self._run_validate_user(
            bootstrap_app, provider, user_info, existing_user_side_effect=[None]
        )

        assert result is None
        assert error[1] == 404
        mock_db.session.add.assert_not_called()

    def test_bootstrap_rejects_unverified_email(self, bootstrap_app, provider):
        user_info = {"sub": "auth0|new-admin", "email": "admin@example.com", "email_verified": False}

        result, error, mock_user_class, mock_db, _new_user = self._run_validate_user(
            bootstrap_app, provider, user_info, existing_user_side_effect=[None]
        )

        assert result is None
        assert error[1] == 404
        mock_db.session.add.assert_not_called()

    def test_bootstrap_noop_when_users_already_exist(self, bootstrap_app, provider):
        user_info = {"sub": "auth0|new-admin", "email": "admin@example.com", "email_verified": True}

        with patch.object(provider, 'get_token', return_value='dummy-token'), \
             patch.object(provider, 'is_token_valid', return_value=True), \
             patch.object(provider, 'get_user', return_value=user_info), \
             patch('mdvtools.dbutils.dbmodels.User') as mock_user_class, \
             patch('mdvtools.dbutils.dbmodels.db') as mock_db:

            mock_db.session = MagicMock()
            filter_by_mock = MagicMock()
            filter_by_mock.first.return_value = None
            mock_user_class.query.filter_by.return_value = filter_by_mock
            mock_user_class.query.count.return_value = 5  # an admin already exists

            with bootstrap_app.test_request_context('/'):
                result, error = provider.validate_user()

        assert result is None
        assert error[1] == 404
        mock_db.session.add.assert_not_called()

    def test_bootstrap_rolls_back_when_role_assignment_fails(self, bootstrap_app, provider):
        user_info = {"sub": "auth0|new-admin", "email": "admin@example.com", "email_verified": True}

        with patch('mdvtools.auth.auth0_provider.GetToken') as mock_get_token, \
             patch('mdvtools.auth.auth0_provider.Auth0') as mock_auth0_class:
            mock_get_token.return_value.client_credentials.return_value = {"access_token": "mgmt-token"}
            mock_auth0 = MagicMock()
            mock_auth0.roles.list.return_value = {"roles": []}  # no 'admin' role in this tenant
            mock_auth0_class.return_value = mock_auth0

            result, error, mock_user_class, mock_db, new_user = self._run_validate_user(
                bootstrap_app, provider, user_info, existing_user_side_effect=[None]
            )

        assert result is None
        assert error[1] == 500
        mock_db.session.delete.assert_called_once_with(new_user)
        assert mock_db.session.commit.call_count == 2  # once for the create, once for the rollback delete

    def test_bootstrap_concurrent_race_logs_in_the_winner(self, bootstrap_app, provider):
        """Two requests for the same bootstrap identity race; the loser's commit hits
        the unique constraint on email/auth_id and should just log in as the winner
        rather than surfacing an error to the legitimate administrator."""
        user_info = {"sub": "auth0|new-admin", "email": "admin@example.com", "email_verified": True}
        winner = MagicMock(id=1, auth_id="auth0|new-admin", email="admin@example.com", is_admin=True)

        with patch.object(provider, 'get_token', return_value='dummy-token'), \
             patch.object(provider, 'is_token_valid', return_value=True), \
             patch.object(provider, 'get_user', return_value=user_info), \
             patch('mdvtools.dbutils.dbmodels.User') as mock_user_class, \
             patch('mdvtools.dbutils.dbmodels.db') as mock_db:

            mock_db.session = MagicMock()
            mock_db.session.commit.side_effect = IntegrityError("insert", {}, Exception("unique violation"))
            filter_by_mock = MagicMock()
            # First call: initial lookup -> None. Second call: after the IntegrityError
            # rollback, re-query finds the row the concurrent request already created.
            filter_by_mock.first.side_effect = [None, winner]
            mock_user_class.query.filter_by.return_value = filter_by_mock
            mock_user_class.query.count.return_value = 0
            mock_user_class.return_value = MagicMock()

            with bootstrap_app.test_request_context('/'):
                result, error = provider.validate_user()

        assert error is None
        assert result == {"id": 1, "auth_id": "auth0|new-admin", "email": "admin@example.com", "is_admin": True}
        mock_db.session.rollback.assert_called_once()

    def test_pending_user_activated_on_first_login(self, bootstrap_app, provider):
        user_info = {"sub": "auth0|invited-user", "email": "invited@example.com", "email_verified": True}
        pending_user = MagicMock(id=7, auth_id="auth0|invited-user", email="invited@example.com", is_admin=False)
        pending_user.is_active = False

        with patch.object(provider, 'get_token', return_value='dummy-token'), \
             patch.object(provider, 'is_token_valid', return_value=True), \
             patch.object(provider, 'get_user', return_value=user_info), \
             patch('mdvtools.dbutils.dbmodels.User') as mock_user_class, \
             patch('mdvtools.dbutils.dbmodels.db') as mock_db:

            mock_db.session = MagicMock()
            filter_by_mock = MagicMock()
            filter_by_mock.first.return_value = pending_user
            mock_user_class.query.filter_by.return_value = filter_by_mock

            with bootstrap_app.test_request_context('/'):
                result, error = provider.validate_user()

        assert error is None
        assert pending_user.is_active is True
        assert pending_user.confirmed_at is not None
        mock_db.session.commit.assert_called_once()
        assert result == {"id": 7, "auth_id": "auth0|invited-user", "email": "invited@example.com", "is_admin": False}

    def test_active_user_login_does_not_touch_db(self, bootstrap_app, provider):
        user_info = {"sub": "auth0|active-user", "email": "active@example.com", "email_verified": True}
        active_user = MagicMock(id=3, auth_id="auth0|active-user", email="active@example.com", is_admin=False)
        active_user.is_active = True

        with patch.object(provider, 'get_token', return_value='dummy-token'), \
             patch.object(provider, 'is_token_valid', return_value=True), \
             patch.object(provider, 'get_user', return_value=user_info), \
             patch('mdvtools.dbutils.dbmodels.User') as mock_user_class, \
             patch('mdvtools.dbutils.dbmodels.db') as mock_db:

            mock_db.session = MagicMock()
            filter_by_mock = MagicMock()
            filter_by_mock.first.return_value = active_user
            mock_user_class.query.filter_by.return_value = filter_by_mock

            with bootstrap_app.test_request_context('/'):
                result, error = provider.validate_user()

        assert error is None
        mock_db.session.commit.assert_not_called()
        assert result == {"id": 3, "auth_id": "auth0|active-user", "email": "active@example.com", "is_admin": False}