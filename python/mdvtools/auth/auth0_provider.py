import time
import requests
from authlib.integrations.flask_client import OAuth
from flask import session, redirect, url_for
from typing import Optional
from mdvtools.auth.auth_provider import AuthProvider
import logging
from jose import jwt
from jose.exceptions import ExpiredSignatureError, JWTError, JWTClaimsError



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
            print("@@@@@@")
            #print(metadata)
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
            print("$$$$$$$$$$$$$$$ -login-1")
            logging.info("Initiating login process.")
            #redirect_uri = url_for('callback', _external=True)
            redirect_uri = self.app.config["AUTH0_CALLBACK_URL"]
            print(redirect_uri)
            print(self.oauth.auth0.authorize_redirect(redirect_uri))
            return self.oauth.auth0.authorize_redirect(redirect_uri=redirect_uri)
        
        except Exception as e:
            logging.error(f"Error during login process: {e}")
            raise RuntimeError("Login failed.") from e

    def logout(self) -> None:
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
            return redirect(logout_url)

        except Exception as e:
            logging.error(f"Error during logout process: {e}")
            raise RuntimeError("Auth0 logout failed.") from e


    def get_user(self, token: dict) -> Optional[dict]:
        """
        Retrieves the user information using the provided token.

        :param token: Dictionary containing access token and user details
        :return: User information dictionary or None
        """
        try:
            logging.info("Fetching user information.")

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
                    "first_name": user_metadata.get("first_name", "Unknown"),
                    "last_name": user_metadata.get("last_name", "Unknown"),
                    "email": raw_data.get("email", ""),
                    "association": user_metadata.get("association", "Unknown Organization"),
                    "avatarUrl": raw_data.get("picture", ""),
                }
                print("!!!!!!")
                print(user_data)
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
            token = self.oauth.auth0.authorize_access_token()
            if 'access_token' not in token:
                raise ValueError("Access token not found in the response.")
            session['token'] = token
            session["auth_method"] = "auth0"
            logging.info("Access token retrieved and stored in session.")
            return token['access_token']
        except Exception as e:
            logging.error(f"Error during callback handling: {e}")
            session.clear() # Clear session in case of failure
            raise RuntimeError("Callback handling failed.") from e

    def is_authenticated(self, token: str) -> bool:
        """
        Checks if the user is authenticated by verifying the token.

        :param token: Access token
        :return: True if authenticated, False otherwise
        """
        try:
            logging.info("Checking authentication status.")
            user_info = self.get_user(token)
            return user_info is not None
        except Exception as e:
            logging.error(f"Error while checking authentication: {e}")
            return False
        
    def is_token_valid(self, token):
        """
        Validates the provided token by verifying its signature using Auth0's public keys
        and ensuring it's not expired.
        """
        try:
            # Step 1: Decode the token header without verification to extract the 'kid'
            unverified_header = jwt.get_unverified_header(token)
            
            if unverified_header is None:
                print("++++++++1")
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
                        print("++++++++2")
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
                print("++++++++3")
                return False

            # Step 3: Verify the JWT token using the public key
            payload = jwt.decode(
                token,
                rsa_key,
                algorithms=["RS256"],
                audience=self.app.config["AUTH0_AUDIENCE"],  # Your API audience
                issuer=f"https://{self.app.config['AUTH0_DOMAIN']}/"
            )
            print("++++++++4")
            # Step 4: Check the expiration of the token
            if payload['exp'] > time.time():
                print("++++++++5")
                return True
            else:
                logging.error("Token is expired.")
                print("++++++++6")
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
