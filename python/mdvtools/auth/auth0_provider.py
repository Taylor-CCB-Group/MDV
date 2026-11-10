from mdvtools._optional import require_extra
require_extra("app", "authlib")

import time
import requests
import threading
from datetime import datetime
from authlib.integrations.flask_client import OAuth
from flask import jsonify, session, redirect
from typing import Optional, Dict, Any
# from flask import Response
from flask.typing import ResponseReturnValue
from mdvtools.auth.auth_provider import AuthProvider, NoAccess
from mdvtools.logging_config import get_logger
from jose import jwt
from jose.exceptions import ExpiredSignatureError, JWTError, JWTClaimsError

logger = get_logger(__name__)

# Add JWKS cache with thread-safe access
_jwks_cache = {}
_jwks_cache_expiry = None
_jwks_cache_lock = threading.Lock()  # Lock to protect cache access
JWKS_CACHE_DURATION = 3600  # Cache for 1 hour

# Prefix of the claims the post-login Action adds (docs/auth0/add-sign-in-claims.js).
CLAIM_NAMESPACE = "https://mdv/"

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
            logger.info("Auth0Provider initialized successfully.")
        except Exception as e:
            logger.critical(f"Failed to initialize Auth0Provider: {e}")
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
                logger.error(f"Failed to fetch OpenID configuration from {server_metadata_url}: {response.text}")
                raise RuntimeError(f"Unable to fetch OpenID Connect metadata from {server_metadata_url}")

            # Parse and check the existence of jwks_uri in the metadata
            metadata = response.json()

            jwks_uri = metadata.get('jwks_uri')
            if not jwks_uri:
                logger.error(f"The OpenID configuration is missing 'jwks_uri': {metadata}")
                raise RuntimeError("'jwks_uri' is missing in OpenID Connect metadata.")

            # Register the OAuth provider with server_metadata_url for dynamic metadata fetching
            self.oauth.register(
                'auth0',
                client_id=self.client_id,
                client_secret=self.client_secret,
                server_metadata_url=server_metadata_url,
                client_kwargs={'scope': 'openid profile email'},
            )
            logger.info("Auth0 OAuth provider registered successfully with OpenID Connect metadata.")
        except Exception as e:
            logger.error(f"Error while registering OAuth provider: {e}")
            raise RuntimeError("Failed to initialize OAuth.") from e

    def login(self) -> str:
        """
        Initiates the login process by redirecting to Auth0's authorization page.
        """
        try:
            
            logger.info("Initiating login process.")
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
            logger.error(f"Error during login process: {e}")
            raise RuntimeError("Login failed.") from e

    def logout(self) -> ResponseReturnValue:
        """
        Logs the user out by clearing the session and redirecting to Auth0's logout endpoint.
        """
        try:
            logger.info("Logging out user from Auth0.")
            
            # Clear the server-side session to remove any stored tokens and user data
            session.clear()
            
            # Prepare the redirect URL after logout (i.e., where the user is sent after logging out of Auth0)
            redirect_url = self.app.config["LOGIN_REDIRECT_URL"]  # The URL to redirect after logout
            
            # Redirect the user to Auth0's logout URL, which will handle the Auth0-side logout
            # This will log the user out of Auth0 and redirect them to the provided URL
            logout_url = f"https://{self.app.config['AUTH0_DOMAIN']}/v2/logout?returnTo={redirect_url}&client_id={self.app.config['AUTH0_CLIENT_ID']}"
            
            logger.info(f"Redirecting to Auth0 logout URL: {logout_url}")
            # "type 'response' is not assignable to None"
            return redirect(logout_url)

        except Exception as e:
            logger.error(f"Error during logout process: {e}")
            raise RuntimeError("Auth0 logout failed.") from e


    def get_user(self, token: Optional[dict] = None) -> Optional[dict]:
        """
        Retrieves the user information using the provided token.

        :param token: Dictionary containing access token and user details
        :return: User information dictionary or None
        """
        try:
            logger.info("Fetching user information.")

            if token is None:
                logger.error("Token is None.")
                return None
        
            # Extract access token
            access_token = token.get("access_token")
            if not access_token:
                logger.error("Access token is missing.")
                return None

            # Correct Authorization Header
            headers = {"Authorization": f"Bearer {access_token}"}

            user_info_url = f"https://{self.domain}/userinfo"
            response = requests.get(user_info_url, headers=headers)

            if response.status_code == 200:
                logger.debug("User information retrieved successfully.")
                raw_data = response.json()

                # Extract user metadata if present
                user_metadata = raw_data.get("user_metadata", {})

                user_data = {
                    "sub": raw_data.get("sub"),
                    "first_name": user_metadata.get("first_name", "Unknown"),
                    "last_name": user_metadata.get("last_name", "Unknown"),
                    "email": raw_data.get("email", ""),
                    "email_verified": bool(raw_data.get("email_verified", False)),
                    # Added by the post-login Action, for accepting an Invitation (ADR 0021).
                    "connection": raw_data.get(f"{CLAIM_NAMESPACE}connection"),
                    "connection_strategy": raw_data.get(f"{CLAIM_NAMESPACE}connection_strategy"),
                    "xms_edov": raw_data.get(f"{CLAIM_NAMESPACE}xms_edov") is True,
                    "association": user_metadata.get("association", "Unknown Organization"),
                    "avatarUrl": raw_data.get("picture", ""),
                }
                
                return user_data
            else:
                logger.warning(f"Failed to fetch user information: {response.status_code} {response.text}")
                return None

        except requests.RequestException as e:
            logger.error(f"Error while fetching user information: {e}")
            return None
        
    def get_token(self) -> Optional[str]:
        """
        Retrieves the token from the session.

        :return: Token string or None
        """
        try:
            logger.info("Retrieving token from session.")
            return session.get('token', {}).get('access_token')
        except Exception as e:
            logger.error(f"Error while retrieving token: {e}")
            return None

    def handle_callback(self) -> Optional[str]:
        """
        Handles the Auth0 callback and retrieves the access token.

        :return: Access token string
        """
        try:
            logger.info("Handling callback from Auth0.")
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
                    logger.error(f"Expected RS256 algorithm, but found {header.get('alg')}")
                    raise ValueError("The token is not of type RS256.")
            except Exception as e:
                logger.error(f"Error decoding the token header: {e}")
                raise ValueError("Invalid token format.")

            # Store the token in the session for later use
            session['token'] = token
            session["auth_method"] = "auth0"
            session.modified = True
            logger.info("Access token retrieved and stored in session.")
            
            return token['access_token']

        except Exception as e:
            logger.error(f"Error during callback handling: {e}")
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
                logger.error("Invalid token header.")
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
                            logger.error(f"Failed to fetch JWKS: {response.status_code}")
                            return False
                        new_jwks = response.json()
                        with _jwks_cache_lock:
                            _jwks_cache = new_jwks
                            _jwks_cache_expiry = current_time + JWKS_CACHE_DURATION
                            cache_to_use = _jwks_cache
                        logger.info("JWKS cache refreshed")
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
                    logger.error(f"Error getting public keys from Auth0: {e}")
                    return False
            
            if not rsa_key:
                logger.error("No valid key found in JWKS for token verification.")
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
                logger.error("Token is expired.")
                return False

        except ExpiredSignatureError:
            logger.error("Token is expired.")
            return False
        except JWTClaimsError:
            logger.error("Invalid claims in token.")
            return False
        except JWTError as e:
            logger.error(f"Error decoding token: {e}")
            return False
        except Exception as e:
            logger.error(f"Error during token validation: {e}")
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
            email = (user_info.get("email") or "").strip().lower()

            # Query the user from the database if not in cache
            user = User.query.filter_by(auth_id=auth0_id).first()

            if user is None:
                try:
                    user = self._try_bootstrap_first_admin(user_info, auth0_id, email)
                except Exception as e:
                    logger.exception(f"Bootstrap administrator creation failed: {e}")
                    return None, (jsonify({"error": "Administrator bootstrap failed"}), 500)
                if user is None:
                    user = self._accept_invitation(user_info)
                if user is None:
                    return None, NoAccess(email=email, sign_in_method=user_info.get("connection"))
            elif not user.is_active:
                # A Deactivated user. Signing in reactivates them only by
                # accepting an Invitation an Admin sent them (ADRs 0020, 0021).
                user = self._accept_invitation(user_info)
                if user is None:
                    return None, NoAccess(email=email, sign_in_method=user_info.get("connection"))
            else:
                # An active User can still hold an Invitation that names them, for
                # example when an Admin reactivated them before they signed in. Their
                # account belongs to them, so no other Invitation can match.
                user = self._accept_invitation(user_info) or user

            # Add the user to the in-memory cache
            user_data = {"id": user.id, "auth_id": user.auth_id, "email": user.email, "is_admin": user.is_admin}

            # Cache the user data in session for future use
            session['user'] = user_data
            session.modified = True

            return user_data, None

        except Exception as e:
            logger.exception(f"Error in validate_user: {e}")
            return None, (jsonify({"error": "Internal server error - user not validated"}), 500)

    def _accept_invitation(self, user_info: Dict[str, Any]):
        """The User made or moved by accepting an open Invitation, or None."""
        from mdvtools.dbutils.dbservice import InvitationService
        from mdvtools.auth.authutils import cache_user_projects

        user = InvitationService.accept_for_sign_in(user_info)
        if user is not None:
            # Project routes read grants from the cache, so the new User would
            # otherwise see no projects until something else refreshed it.
            cache_user_projects()
        return user

    def _try_bootstrap_first_admin(self, user_info: Dict[str, Any], auth0_id: Optional[str], email: str):
        """
        Create the first local administrator, but only while the local database has zero
        users and the verified Auth0 identity matches MDV_BOOTSTRAP_ADMIN_EMAIL. Returns
        None (bootstrap not applicable) rather than creating anything if any check fails.

        This path only ever runs once per deployment: as soon as any local user exists,
        User.query.count() > 0 and this method always returns None.
        """
        from mdvtools.dbutils.dbmodels import User, db
        from sqlalchemy.exc import IntegrityError

        if not auth0_id:
            return None

        bootstrap_email = (self.app.config.get("MDV_BOOTSTRAP_ADMIN_EMAIL") or "").strip().lower()
        if not bootstrap_email:
            return None
        if User.query.count() > 0:
            return None
        if not user_info.get("email_verified"):
            logger.warning("Bootstrap admin login rejected: Auth0 email is not verified.")
            return None
        if not email or email != bootstrap_email:
            logger.warning("Bootstrap admin login rejected: email does not match MDV_BOOTSTRAP_ADMIN_EMAIL.")
            return None

        user = User(
            email=email,
            auth_id=auth0_id,
            confirmed_at=datetime.now(),
            is_active=True,
            administrator=True,
            is_admin=True,
            password="",
        )
        db.session.add(user)
        # Admin rights live only in MDV (ADR 0019), so no Auth0 role is assigned.
        try:
            db.session.flush()
        except IntegrityError:
            db.session.rollback()
            # A concurrent request for this same identity already won the race and
            # created the row (e.g. a double-click or two tabs). Log that user in
            # rather than showing the legitimate bootstrap administrator an error.
            return User.query.filter_by(auth_id=auth0_id).first()

        try:
            db.session.commit()
        except IntegrityError:
            db.session.rollback()
            # The same race, seen at commit instead of at flush.
            return User.query.filter_by(auth_id=auth0_id).first()

        # Administrators own every project, but only Admin's Make Admin and a
        # project rescan write those grants, and neither has run on a fresh deployment.
        # Without this the first administrator logs in to an empty project list.
        try:
            from mdvtools.dbutils.dbservice import UserProjectService
            from mdvtools.auth.authutils import cache_user_projects

            granted = UserProjectService.grant_all_projects_to_admins()
            cache_user_projects()
            logger.info(f"Granted the bootstrap administrator access to {granted} project(s).")
        except Exception:
            # To rollback the session if the failure occurs in the querying of grant_all_projects_to_admins()
            db.session.rollback()
            logger.exception(
                "Bootstrapped the first administrator but could not grant project access; "
                "it can still be granted from Admin."
            )

        logger.info(f"Bootstrapped first administrator: {email}")
        return user
