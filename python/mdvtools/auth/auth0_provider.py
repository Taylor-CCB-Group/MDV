import time
import requests
import threading
from authlib.integrations.flask_client import OAuth
from flask import jsonify, session, redirect
from typing import Optional, Dict, List, Any, TYPE_CHECKING

if TYPE_CHECKING:
    from mdvtools.dbutils.dbmodels import Project
# from flask import Response
from flask.typing import ResponseReturnValue
from mdvtools.auth.auth_provider import AuthProvider
import logging
from jose import jwt
from jose.exceptions import ExpiredSignatureError, JWTError, JWTClaimsError
from auth0.management import Auth0
from auth0.authentication import GetToken
from auth0.exceptions import RateLimitError
import random

# Add JWKS cache with thread-safe access
_jwks_cache = {}
_jwks_cache_expiry = None
_jwks_cache_lock = threading.Lock()  # Lock to protect cache access
JWKS_CACHE_DURATION = 3600  # Cache for 1 hour

# Rate limiting and retry parameters
BASE_DELAY = 1  # Base delay in seconds
MAX_DELAY = 8  # Maximum delay in seconds

class Auth0Provider(AuthProvider):
    def __init__(self, app, oauth: OAuth, client_id: str, client_secret: str, domain: str):
        """
        Initializes the Auth0Provider class with application details.

        :param app: Flask app instance
        :param oauth: Authlib OAuth instance
        :param client_id: Auth0 Client ID
        :param client_secret: Auth0 Client Secret
        :param domain: Auth0 Domain
        """
        try:
            if not all([client_id, client_secret, domain]):
                raise ValueError("Missing required Auth0 configuration parameters.")
            
            self.app = app
            self.oauth = oauth
            self.client_id = client_id
            self.client_secret = client_secret
            self.domain = domain

            self._initialize_oauth()
            logging.info("Auth0Provider initialized successfully.")
        except Exception as e:
            logging.critical(f"Failed to initialize Auth0Provider: {e}")
            raise

    def _initialize_oauth(self):
        """
        Registers the Auth0 OAuth provider and validates OpenID Connect metadata.
        """
        try:
            # Construct the server metadata URL for OpenID Connect discovery
            server_metadata_url = f'https://{self.domain}/.well-known/openid-configuration'

            # Attempt to fetch metadata to ensure it's accessible
            response = requests.get(server_metadata_url)
            if response.status_code != 200:
                logging.error(f"Failed to fetch OpenID configuration from {server_metadata_url}: {response.text}")
                raise RuntimeError(f"Unable to fetch OpenID Connect metadata from {server_metadata_url}")

            # Parse and check the existence of jwks_uri in the metadata
            metadata = response.json()

            jwks_uri = metadata.get('jwks_uri')
            if not jwks_uri:
                logging.error(f"The OpenID configuration is missing 'jwks_uri': {metadata}")
                raise RuntimeError("'jwks_uri' is missing in OpenID Connect metadata.")

            # Register the OAuth provider with server_metadata_url for dynamic metadata fetching
            self.oauth.register(
                'auth0',
                client_id=self.client_id,
                client_secret=self.client_secret,
                server_metadata_url=server_metadata_url,
                client_kwargs={'scope': 'openid profile email'},
            )
            logging.info("Auth0 OAuth provider registered successfully with OpenID Connect metadata.")
        except Exception as e:
            logging.error(f"Error while registering OAuth provider: {e}")
            raise RuntimeError("Failed to initialize OAuth.") from e

    def login(self) -> str:
        """
        Initiates the login process by redirecting to Auth0's authorization page.
        """
        try:
            
            logging.info("Initiating login process.")
            #redirect_uri = url_for('callback', _external=True)
            redirect_uri = self.app.config["AUTH0_CALLBACK_URL"]
            audience = self.app.config["AUTH0_AUDIENCE"]  # The API audience for which the token is requested
            

            # Initiate the redirect to Auth0's authorization endpoint with necessary parameters
            assert self.oauth.auth0 is not None, "Auth0 provider is not registered."
            return self.oauth.auth0.authorize_redirect(
                redirect_uri=redirect_uri,
                audience=audience  # The audience for the token (API identifier)                
            )
        
        except Exception as e:
            logging.error(f"Error during login process: {e}")
            raise RuntimeError("Login failed.") from e

    def logout(self) -> ResponseReturnValue:
        """
        Logs the user out by clearing the session and redirecting to Auth0's logout endpoint.
        """
        try:
            logging.info("Logging out user from Auth0.")
            
            # Clear the server-side session to remove any stored tokens and user data
            session.clear()
            
            # Prepare the redirect URL after logout (i.e., where the user is sent after logging out of Auth0)
            redirect_url = self.app.config["LOGIN_REDIRECT_URL"]  # The URL to redirect after logout
            
            # Redirect the user to Auth0's logout URL, which will handle the Auth0-side logout
            # This will log the user out of Auth0 and redirect them to the provided URL
            logout_url = f"https://{self.app.config['AUTH0_DOMAIN']}/v2/logout?returnTo={redirect_url}&client_id={self.app.config['AUTH0_CLIENT_ID']}"
            
            logging.info(f"Redirecting to Auth0 logout URL: {logout_url}")
            # "type 'response' is not assignable to None"
            return redirect(logout_url)

        except Exception as e:
            logging.error(f"Error during logout process: {e}")
            raise RuntimeError("Auth0 logout failed.") from e


    def get_user(self, token: Optional[dict] = None) -> Optional[dict]:
        """
        Retrieves the user information using the provided token.

        :param token: Dictionary containing access token and user details
        :return: User information dictionary or None
        """
        try:
            logging.info("Fetching user information.")

            if token is None:
                logging.error("Token is None.")
                return None
        
            # Extract access token
            access_token = token.get("access_token")
            if not access_token:
                logging.error("Access token is missing.")
                return None

            # Correct Authorization Header
            headers = {"Authorization": f"Bearer {access_token}"}

            user_info_url = f"https://{self.domain}/userinfo"
            response = requests.get(user_info_url, headers=headers)

            if response.status_code == 200:
                logging.debug("User information retrieved successfully.")
                raw_data = response.json()

                # Extract user metadata if present
                user_metadata = raw_data.get("user_metadata", {})

                user_data = {
                    "sub": raw_data.get("sub"),
                    "first_name": user_metadata.get("first_name", "Unknown"),
                    "last_name": user_metadata.get("last_name", "Unknown"),
                    "email": raw_data.get("email", ""),
                    "association": user_metadata.get("association", "Unknown Organization"),
                    "avatarUrl": raw_data.get("picture", ""),
                }
                
                return user_data
            else:
                logging.warning(f"Failed to fetch user information: {response.status_code} {response.text}")
                return None

        except requests.RequestException as e:
            logging.error(f"Error while fetching user information: {e}")
            return None
        
    def get_token(self) -> Optional[str]:
        """
        Retrieves the token from the session.

        :return: Token string or None
        """
        try:
            logging.info("Retrieving token from session.")
            return session.get('token', {}).get('access_token')
        except Exception as e:
            logging.error(f"Error while retrieving token: {e}")
            return None

    def handle_callback(self) -> Optional[str]:
        """
        Handles the Auth0 callback and retrieves the access token.

        :return: Access token string
        """
        try:
            logging.info("Handling callback from Auth0.")
            assert self.oauth.auth0 is not None, "Auth0 provider is not registered."
            token = self.oauth.auth0.authorize_access_token()
            if 'access_token' not in token:
                raise ValueError("Access token not found in the response.")
            # Check if the token is of the expected type (JWT)
            access_token = token['access_token']
            
            # Decode the token's header to inspect its algorithm and other details
            try:
                header = jwt.get_unverified_header(access_token)
                #print("Token Header:", header)  # Check the header of the token

                # Ensure the algorithm is RS256 (not JWE)
                if header.get('alg') != 'RS256':
                    logging.error(f"Expected RS256 algorithm, but found {header.get('alg')}")
                    raise ValueError("The token is not of type RS256.")
            except Exception as e:
                logging.error(f"Error decoding the token header: {e}")
                raise ValueError("Invalid token format.")

            # Store the token in the session for later use
            session['token'] = token
            session["auth_method"] = "auth0"
            session.modified = True
            logging.info("Access token retrieved and stored in session.")
            
            return token['access_token']

        except Exception as e:
            logging.error(f"Error during callback handling: {e}")
            session.clear()  # Clear session in case of failure
            raise RuntimeError("Callback handling failed.") from e
        
    def is_token_valid(self, token):
        """
        Validates the provided token by verifying its signature using Auth0's public keys
        and ensuring it's not expired.
        """
        global _jwks_cache, _jwks_cache_expiry, _jwks_cache_lock
        
        try:
            # Step 1: Decode the token header without verification to extract the 'kid'
            unverified_header = jwt.get_unverified_header(token)
            
            if unverified_header is None:
                logging.error("Invalid token header.")
                return False

            # Step 2: Get the public key from Auth0's JWKS (JSON Web Key Set) endpoint with caching
            rsa_key = {}
            if 'kid' in unverified_header:
                try:
                    # Double-checked locking pattern for thread-safe cache access
                    current_time = time.time()
                    
                    # First check (without lock) - fast path for cache hit
                    if _jwks_cache_expiry is not None and current_time <= _jwks_cache_expiry:
                        # Cache is valid, read it (this is safe for reads in Python due to GIL,
                        # but we still need lock for consistency with writes)
                        with _jwks_cache_lock:
                            # Re-check after acquiring lock (double-checked locking)
                            if _jwks_cache_expiry is not None and current_time <= _jwks_cache_expiry:
                                # Cache is still valid, use it
                                cache_to_use = _jwks_cache
                            else:
                                # Cache expired while waiting for lock, need to refresh
                                cache_to_use = None
                    else:
                        # Cache expired or doesn't exist, need to refresh
                        cache_to_use = None
                    
                    # If cache is invalid or expired, refresh it
                    needs_refresh = False
                    if cache_to_use is None:
                        with _jwks_cache_lock:
                            # Re-check one more time (another thread might have refreshed it)
                            if _jwks_cache_expiry is not None and current_time <= _jwks_cache_expiry:
                                # Another thread refreshed it, use the cached version
                                cache_to_use = _jwks_cache
                            else:
                                # We need to refresh the cache
                                needs_refresh = True
                    if needs_refresh:
                        response = requests.get(
                            url=self.app.config['AUTH0_PUBLIC_KEY_URI'],
                            timeout=10
                        )
                        if response.status_code != 200:
                            logging.error(f"Failed to fetch JWKS: {response.status_code}")
                            return False
                        new_jwks = response.json()
                        with _jwks_cache_lock:
                            _jwks_cache = new_jwks
                            _jwks_cache_expiry = current_time + JWKS_CACHE_DURATION
                            cache_to_use = _jwks_cache
                        logging.info("JWKS cache refreshed")
                    # Find the key in the cached JWKS that matches the 'kid' in the token header
                    assert cache_to_use is not None, "Cache is None"
                    for key in cache_to_use['keys']:
                        if key['kid'] == unverified_header['kid']:
                            rsa_key = {
                                'kty': key['kty'],
                                'kid': key['kid'],
                                'use': key['use'],
                                'n': key['n'],
                                'e': key['e']
                            }
                            break
                except Exception as e:
                    logging.error(f"Error getting public keys from Auth0: {e}")
                    return False
            
            if not rsa_key:
                logging.error("No valid key found in JWKS for token verification.")
                return False

            # Step 3: Verify the JWT token using the public key
            payload = jwt.decode(
                token,
                rsa_key,
                algorithms=["RS256"],
                audience=self.app.config["AUTH0_AUDIENCE"],  # Your API audience
                issuer=f"https://{self.app.config['AUTH0_DOMAIN']}/"
            )
    
            # Step 4: Check the expiration of the token
            if payload['exp'] > time.time():
                return True
            else:
                logging.error("Token is expired.")
                return False

        except ExpiredSignatureError:
            logging.error("Token is expired.")
            return False
        except JWTClaimsError:
            logging.error("Invalid claims in token.")
            return False
        except JWTError as e:
            logging.error(f"Error decoding token: {e}")
            return False
        except Exception as e:
            logging.error(f"Error during token validation: {e}")
            return False

    def validate_user(self):
        """Validate the user using Auth0."""

        from mdvtools.dbutils.dbmodels import User
        
        try:
            # Check if user information is already cached in session
            if 'user' in session:
                return session['user'], None  # Return the user from session cache

            # Retrieve the token from session
            token = self.get_token()
            if not token:
                return None, (jsonify({"error": "Authentication required"}), 401)

            # Validate token using the provider-specific logic
            if not self.is_token_valid(token):
                return None, (jsonify({"error": "Invalid or expired token"}), 401)
            
            # Retrieve user info from Auth0
            user_info = self.get_user({"access_token": token})
            if not user_info:
                return None, (jsonify({"error": "User not found"}), 404)

            # Get Auth0 user ID
            auth0_id = user_info.get("sub")

            # Query the user from the database if not in cache
            user = User.query.filter_by(auth_id=auth0_id).first()
            if not user:
                return None, (jsonify({"error": "User not found"}), 404)

            # Add the user to the in-memory cache
            user_data = {"id": user.id, "auth_id": user.auth_id, "email": user.email, "is_admin": user.is_admin}

            # Cache the user data in session for future use
            session['user'] = user_data
            session.modified = True

            return user_data, None

        except Exception as e:
            logging.exception(f"Error in validate_user: {e}")
            return None, (jsonify({"error": "Internal server error - user not validated"}), 500)
    
    class SyncContext:
        """Context object for user sync operations."""
        def __init__(
            self,
            auth0: Auth0,
            all_projects: List['Project'],
            page: int,
            processed_users: int = 0,
            new_users: int = 0,
            updated_users: int = 0,
            admin_users_synced: int = 0,
            pagination_rate_limit_count: int = 0,
            max_pagination_rate_limit_retries: int = 5,
            pagination_retry_delay: float = 2.0,
            per_page: int = 50,
            max_role_fetch_retries: int = 3
        ):
            self.auth0 = auth0
            self.all_projects = all_projects
            self.page = page
            self.processed_users = processed_users
            self.new_users = new_users
            self.updated_users = updated_users
            self.admin_users_synced = admin_users_synced
            self.pagination_rate_limit_count = pagination_rate_limit_count
            self.max_pagination_rate_limit_retries = max_pagination_rate_limit_retries
            self.pagination_retry_delay = pagination_retry_delay
            self.per_page = per_page
            self.max_role_fetch_retries = max_role_fetch_retries
        
        def to_dict(self) -> Dict[str, int]:
            """Convert context stats to dictionary."""
            return {
                'processed_users': self.processed_users,
                'new_users': self.new_users,
                'updated_users': self.updated_users,
                'admin_users_synced': self.admin_users_synced
            }
        
        def reset_pagination_rate_limit(self) -> None:
            """Reset pagination rate limit counters after successful request."""
            self.pagination_rate_limit_count = 0
            self.pagination_retry_delay = 2.0
    
    def _fetch_user_roles_with_retry(
        self, 
        auth0_id: str,
        context: 'SyncContext'
    ) -> Optional[Dict[str, Any]]:
        """
        Fetch user roles from Auth0 with retry logic for rate limiting.
        
        Args:
            auth0_id: Auth0 user ID
            context: SyncContext object containing auth0 client, page, and stats
            
        Returns:
            dict: User roles dictionary, or None if fetch failed after retries
        """
        roles = None
        retry_count = 0
        success = False
        max_retries = context.max_role_fetch_retries
        
        while retry_count < max_retries and not success:
            try:
                roles = context.auth0.users.list_roles(auth0_id)
                success = True
            except RateLimitError as e:
                retry_count += 1
                if retry_count >= max_retries:
                    # Max retries exceeded
                    error_info = {
                        'user_id': auth0_id,
                        'page': context.page,
                        'processed': context.processed_users,
                        'status': getattr(e, 'status_code', None),
                        'code': getattr(e, 'error_code', None),
                        'message': getattr(e, 'message', str(e)),
                        'reset_at': getattr(e, 'reset_at', None),
                        'remaining': getattr(e, 'remaining', None),
                        'limit': getattr(e, 'limit', None),
                    }
                    error_info = {k: v for k, v in error_info.items() if v is not None}
                    
                    logging.error(
                        f"Rate limit during role fetch for user {auth0_id} after {max_retries} retries: "
                        f"{error_info}. Skipping user."
                    )
                    return None
                else:
                    # Calculate exponential backoff delay with jitter
                    delay = min(BASE_DELAY * (2 ** retry_count) + random.uniform(0, 1), MAX_DELAY)
                    logging.warning(
                        f"Rate limit during role fetch for user {auth0_id}. "
                        f"Retrying in {delay:.2f} seconds... (Attempt {retry_count}/{max_retries})"
                    )
                    time.sleep(delay)
            except Exception as e:
                # Non-rate-limit error
                logging.error(f"Error fetching roles for user {auth0_id}: {str(e)}")
                return None
        
        return roles if success else None
    
    def _process_single_user(
        self, 
        user: Dict[str, Any], 
        context: 'SyncContext'
    ) -> bool:
        """
        Process a single user: sync to DB, fetch roles, update admin status, assign projects.
        
        Args:
            user: User dict from Auth0 API
            context: SyncContext object containing auth0 client, projects, page, and stats
            
        Returns:
            bool: True if user was successfully processed, False otherwise
        """
        from mdvtools.dbutils.dbservice import UserService, UserProjectService
        from mdvtools.dbutils.dbmodels import db, User
        
        email = user.get('email', '')
        auth0_id = user['user_id']
        
        # Check if user already exists to track new vs updated
        existing_user = User.query.filter_by(auth_id=auth0_id).first()
        is_new_user = existing_user is None
        
        # Use UserService to add or update user
        db_user = UserService.add_or_update_user(
            email=email,
            auth_id=auth0_id
        )
                
        # Add delay between role requests to avoid rate limiting
        time.sleep(0.2)  # 200ms delay between requests
        
        # Fetch user's roles with retry mechanism
        roles = self._fetch_user_roles_with_retry(auth0_id, context)
        if roles is None:
            # Failed to fetch roles, skip this user
            return False
        
        is_admin = any(role['name'] == 'admin' for role in roles['roles'])
        
        # Update admin status
        was_admin = db_user.is_admin
        db_user.is_admin = is_admin
        try:
            db.session.commit()
        except Exception as e:
            db.session.rollback()
            logging.error(f"Error updating user {auth0_id} admin status: {e}")
            return False
                
        if is_admin:
            # Assign all projects to this user as owner via UserProjectService
            for project in context.all_projects:
                UserProjectService.add_or_update_user_project(
                    user_id=db_user.id,
                    project_id=project.id,
                    is_owner=True
                )
        # Update context stats
        if is_admin and not was_admin:
            context.admin_users_synced += 1
        if is_new_user:
            context.new_users += 1
        else:
            context.updated_users += 1        
        context.processed_users += 1
        return True
    
    def _handle_pagination_rate_limit(
        self, 
        e: RateLimitError,
        context: 'SyncContext'
    ) -> None:
        """
        Handle rate limit error during pagination. Updates context in place.
        
        Args:
            e: RateLimitError exception
            context: SyncContext object containing pagination state
            
        Raises:
            RuntimeError: If max retries exceeded
        """
        context.pagination_rate_limit_count += 1
        new_count = context.pagination_rate_limit_count
        
        # Collect all available error information
        error_info = {
            'page': context.page,
            'processed': context.processed_users,
            'consecutive_errors': new_count,
            'retry_delay': context.pagination_retry_delay,
            'status': getattr(e, 'status_code', None),
            'code': getattr(e, 'error_code', None),
            'message': getattr(e, 'message', str(e)),
            'reset_at': getattr(e, 'reset_at', None),
            'remaining': getattr(e, 'remaining', None),
            'limit': getattr(e, 'limit', None),
        }
        # Filter out None values for cleaner output
        error_info = {k: v for k, v in error_info.items() if v is not None}
        
        logging.error(f"Rate limit during pagination: {error_info}")
        
        # Check if we've exceeded max retries
        if new_count >= context.max_pagination_rate_limit_retries:
            logging.error(
                f"Aborting: {new_count} consecutive rate limit errors "
                f"(max: {context.max_pagination_rate_limit_retries}). Processed {context.processed_users} users."
            )
            raise RuntimeError(
                f"Too many consecutive rate limit errors. Check: function call frequency, "
                f"concurrent execution, or insufficient delays between requests."
            ) from e
        
        logging.warning(
            f"Rate limit hit ({new_count}/{context.max_pagination_rate_limit_retries}). "
            f"Possible issues: frequent calls, concurrent execution, batch size ({context.per_page}), "
            f"or insufficient delays (1s pages, 0.2s roles). Retrying in {context.pagination_retry_delay}s..."
        )
    
    def _log_sync_statistics(self, initial_user_count: int, initial_admin_count: int) -> None:
        """
        Log initial and final database statistics for user sync.
        
        Args:
            initial_user_count: User count at start of sync
            initial_admin_count: Admin user count at start of sync
        """
        from mdvtools.dbutils.dbmodels import User
        
        final_user_count = User.query.count()
        final_admin_count = User.query.filter_by(is_admin=True).count()
        
        logging.info(
            f"User sync completed. Database: {final_user_count} total users "
            f"({final_user_count - initial_user_count:+d}), "
            f"{final_admin_count} admin users ({final_admin_count - initial_admin_count:+d})"
        )
        
    def sync_users_to_db(self) -> None:
        """
        Syncs users from Auth0 to the application's database using UserService and UserProjectService.
        Implements rate limiting and retry logic for Auth0 API calls.
        
        WARNING: This function makes many Auth0 Management API calls and should be called sparingly.
        Known call sites:
        - On application startup (acceptable)
        - From manage_project_permissions.py script (acceptable)
        - From project_manager_extension when showing share dialog (PROBLEMATIC - may cause rate limiting)
        
        Design concerns:
        - Auth0 Management API has strict rate limits (typically 2 req/sec for free tier)
        - This function makes 1 + N*2 API calls where N = number of users (list users + list_roles per user)
        - If called frequently (e.g., on every share dialog open), will hit rate limits quickly
        - Consider: caching, background jobs, or incremental sync instead of full sync on-demand
        """
        from mdvtools.dbutils.dbmodels import Project
        
        try:
            # Load Auth0 config from app
            auth0_domain = self.app.config['AUTH0_DOMAIN']
            client_id = self.app.config['AUTH0_CLIENT_ID']
            client_secret = self.app.config['AUTH0_CLIENT_SECRET']
            auth0_db_connection = self.app.config['AUTH0_DB_CONNECTION']
            audience = f"https://{auth0_domain}/api/v2/"

            # Get Auth0 Management API token
            get_token = GetToken(domain=auth0_domain, client_id=client_id, client_secret=client_secret)
            mgmt_api_token = get_token.client_credentials(audience=audience)["access_token"]
            auth0 = Auth0(auth0_domain, mgmt_api_token)

            # Get initial database statistics
            from mdvtools.dbutils.dbmodels import User
            initial_user_count = User.query.count()
            initial_admin_count = User.query.filter_by(is_admin=True).count()
            logging.info(
                f"Starting user sync. Database stats: {initial_user_count} total users, "
                f"{initial_admin_count} admin users"
            )
            
            # Initialize sync context
            all_projects = Project.query.all()
            sync_context = self.SyncContext(
                auth0=auth0,
                all_projects=all_projects,
                page=0
            )
            
            # Fetch users from Auth0 connection with pagination
            while True:
                try:
                    users = auth0.users.list(
                        q=f'identities.connection:"{auth0_db_connection}"',
                        page=sync_context.page,
                        per_page=sync_context.per_page
                    )
                    
                    # Reset rate limit counter on successful request
                    sync_context.reset_pagination_rate_limit()
                    
                    # Check if we've reached the end of pagination
                    # Auth0 returns an empty list when there are no more users to fetch
                    user_list = users.get('users', [])
                    if not user_list:
                        logging.info(f"Reached end of pagination at page {sync_context.page} (empty user list returned)")
                        break
                    
                    # Process users on this page
                    for user in user_list:
                        success = False
                        try:
                            success = self._process_single_user(user, sync_context)
                        except Exception as e:
                            logging.error(f"Error processing user {user['user_id']}: {e}")
                        
                        if success:
                            # Log progress every 10 users
                            if sync_context.processed_users % 10 == 0:
                                logging.info(
                                    f"Progress: {sync_context.processed_users} processed "
                                    f"({sync_context.new_users} new, {sync_context.updated_users} updated, "
                                    f"{sync_context.admin_users_synced} admins)"
                                )
                    
                    # Check if we got fewer users than requested (indicates last page)
                    if len(user_list) < sync_context.per_page:
                        logging.info(
                            f"Reached end of pagination at page {sync_context.page} "
                            f"(got {len(user_list)} users, less than requested {sync_context.per_page})"
                        )
                        break
                    
                    # Continue to next page
                    sync_context.page += 1
                    time.sleep(1)  # Add delay between pagination requests
                    
                except RateLimitError as e:
                    self._handle_pagination_rate_limit(e, sync_context)
                    time.sleep(sync_context.pagination_retry_delay)
                    # Update retry delay with exponential backoff, max 60s
                    sync_context.pagination_retry_delay = min(sync_context.pagination_retry_delay * 2, 60.0)
                    continue
                
            # Log final statistics
            logging.info(
                f"User sync completed. Stats: {sync_context.processed_users} processed "
                f"({sync_context.new_users} new, {sync_context.updated_users} updated, "
                f"{sync_context.admin_users_synced} admins synced)."
            )
            self._log_sync_statistics(initial_user_count, initial_admin_count)

        except Exception as e:
            logging.exception(f"In sync_users_to_db: An unexpected error occurred: {e}")
            raise