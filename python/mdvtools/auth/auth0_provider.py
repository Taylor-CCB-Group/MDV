import requests
from authlib.integrations.flask_client import OAuth
from flask import session, redirect, url_for
from typing import Optional
from mdvtools.auth.auth_provider import AuthProvider
import logging


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
            print(metadata)
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
            redirect_uri = "https://bia.cmd.ox.ac.uk/carroll/callback"
            print(redirect_uri)
            print(self.oauth.auth0.authorize_redirect(redirect_uri))
            return self.oauth.auth0.authorize_redirect(redirect_uri)
        except Exception as e:
            logging.error(f"Error during login process: {e}")
            raise RuntimeError("Login failed.") from e

    def logout(self) -> None:
        """
        Logs the user out by clearing the session.
        """
        try:
            logging.info("Logging out user.")
            session.clear()
        except Exception as e:
            logging.error(f"Error during logout process: {e}")
            raise RuntimeError("Logout failed.") from e

    def get_user(self, token: str) -> Optional[dict]:
        """
        Retrieves the user information using the provided token.

        :param token: Access token
        :return: User information dictionary or None
        """
        try:
            logging.info("Fetching user information.")
            user_info_url = f'https://{self.domain}/userinfo'
            response = requests.get(user_info_url, headers={'Authorization': f'Bearer {token}'})
            if response.status_code == 200:
                logging.debug("User information retrieved successfully.")
                return response.json()
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
            return session.get('token')
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
            logging.info("Access token retrieved and stored in session.")
            return token['access_token']
        except Exception as e:
            logging.error(f"Error during callback handling: {e}")
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
