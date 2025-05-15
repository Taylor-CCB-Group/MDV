import time
import requests
from authlib.integrations.flask_client import OAuth
from flask import jsonify, session, redirect
from typing import Optional, Tuple, Union
from flask import Response
from mdvtools.auth.auth_provider import AuthProvider
import logging
from jose import jwt
from jose.exceptions import ExpiredSignatureError, JWTError, JWTClaimsError
from auth0.management import Auth0
from auth0.authentication import GetToken


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
            print(redirect_uri)
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

    def logout(self) -> Union[str, Response, Tuple[Response, int]]:
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
        try:
            # Step 1: Decode the token header without verification to extract the 'kid'
            unverified_header = jwt.get_unverified_header(token)
            
            if unverified_header is None:
    
                logging.error("Invalid token header.")
                return False

            # Step 2: Get the public key from Auth0's JWKS (JSON Web Key Set) endpoint
            rsa_key = {}
            if 'kid' in unverified_header:
                try:
                    # Fetch Auth0 public keys from jwks_uri
                    response = requests.get(self.app.config['AUTH0_PUBLIC_KEY_URI'])
                    if response.status_code != 200:
                        logging.error(f"Failed to fetch JWKS: {response.status_code}")
                        
                        return False
                    jwks = response.json()

                
                    # Find the key in the JWKS that matches the 'kid' in the token header
                    for key in jwks['keys']:
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
        
    def sync_users_to_db(self):
        """
        Syncs users from Auth0 to the application's database using UserService and UserProjectService.
        """
        from mdvtools.dbutils.dbservice import UserService, UserProjectService
        from mdvtools.dbutils.dbmodels import db, Project
        
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

            # Fetch users from Auth0 connection
            users = auth0.users.list(q=f'identities.connection:"{auth0_db_connection}"', per_page=100)

            for user in users['users']:
                email = user.get('email', '')
                auth0_id = user['user_id']

                # Use UserService to add or update user
                db_user = UserService.add_or_update_user(
                    email=email,
                    auth_id=auth0_id
                )

                # Fetch user's roles to determine admin status
                roles = auth0.users.list_roles(auth0_id)
                is_admin = any(role['name'] == 'admin' for role in roles['roles'])

                # Update admin status
                db_user.is_admin = is_admin
                db.session.commit()

                if is_admin:
                    # Assign all projects to this user as owner via UserProjectService
                    for project in Project.query.all():
                        UserProjectService.add_or_update_user_project(
                            user_id=db_user.id,
                            project_id=project.id,
                            is_owner=True
                        )

            logging.info("Successfully synced users from Auth0 to the database.")

        except Exception as e:
            logging.exception(f"In sync_users_to_db: An unexpected error occurred: {e}")
            raise