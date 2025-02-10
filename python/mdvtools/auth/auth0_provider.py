import requests
import logging
from datetime import datetime, timedelta
from authlib.integrations.flask_client import OAuth
from flask import session, redirect, url_for
from typing import Optional
from mdvtools.auth.auth_provider import AuthProvider

class Auth0Provider(AuthProvider):
    SESSION_TIMEOUT = timedelta(minutes=30)  # Set session timeout (adjust as needed)

    def __init__(self, app, oauth: OAuth, client_id: str, client_secret: str, domain: str):
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

    def login(self) -> str:
        """Initiates the login process and stores session timestamp."""
        try:
            session.clear()  # Clear old session
            redirect_uri = url_for('callback', _external=True)
            response = self.oauth.auth0.authorize_redirect(redirect_uri)
            
            # Store login timestamp for session expiration handling
            session['login_timestamp'] = datetime.utcnow().isoformat()

            return response
        except Exception as e:
            logging.error(f"Error during login process: {e}")
            raise RuntimeError("Login failed.") from e

    def logout(self) -> None:
        """Logs the user out by clearing the session."""
        try:
            logging.info("Logging out user.")
            session.clear()
        except Exception as e:
            logging.error(f"Error during logout process: {e}")
            raise RuntimeError("Logout failed.") from e

    def handle_callback(self) -> Optional[str]:
        """Handles the Auth0 callback and retrieves the access token."""
        try:
            logging.info("Handling callback from Auth0.")
            token = self.oauth.auth0.authorize_access_token()

            if 'access_token' not in token:
                raise ValueError("Access token not found in response.")

            # Store access token and login timestamp
            session['token'] = token
            session['login_timestamp'] = datetime.utcnow().isoformat()
            
            logging.info("Access token retrieved and stored in session.")
            return token['access_token']
        except Exception as e:
            logging.error(f"Error during callback handling: {e}")
            raise RuntimeError("Callback handling failed.") from e

    @staticmethod
    def is_authenticated() -> bool:
        """Checks if the user is authenticated and session is valid."""
        token = session.get('token')
        login_time_str = session.get('login_timestamp')

        if not token or not login_time_str:
            return False  # No valid session

        try:
            login_time = datetime.fromisoformat(login_time_str)
            if datetime.utcnow() - login_time > Auth0Provider.SESSION_TIMEOUT:
                logging.info("Session expired, forcing logout.")
                session.clear()
                return False
        except ValueError:
            logging.error("Invalid session timestamp format, forcing logout.")
            session.clear()
            return False

        return True
