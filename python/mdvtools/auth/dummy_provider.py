from flask import session, redirect
from typing import Optional, Tuple, Union
from flask import Response
from mdvtools.auth.auth_provider import AuthProvider
import logging

logger = logging.getLogger(__name__)

class DummyAuthProvider(AuthProvider):
    def __init__(self, app):
        self.app = app

    def login(self) -> str:
        """No-op login for dummy provider."""

        session["auth_method"] = "dummy"
        session.modified = True
        logger.info("Dummy login - no real login redirect.")
        return redirect("/")

    def logout(self) -> Union[str, Response, Tuple[Response, int]]:
        """Clear dummy session."""
        logger.info("Dummy logout - clearing session.")
        session.clear()
        redirect_uri = self.app.config.get("LOGIN_REDIRECT_URL", "/login_dev")
        return redirect(redirect_uri)

    def get_user(self, token: Optional[dict] = None) -> Optional[dict]:
        """Return dummy user data."""
        logger.debug("Returning dummy user data.")
        return {
            "id": 1,
            "auth_id": "dummy|localuser",
            "email": "dev@example.com",
            "is_admin": True,
            "first_name": "Dev",
            "last_name": "User"
        }

    def get_token(self) -> Optional[str]:
        """Dummy token (not used)."""
        return "dummy_token"

    def handle_callback(self) -> Optional[str]:
        """Simulate callback."""
        logger.info("Dummy callback - no real token exchange.")
        return "dummy_token"

    def validate_user(self) -> Tuple[Optional[dict], Optional[Tuple]]:
        """Simulate user validation and cache dummy user."""
        user_data = self.get_user()
        session['user'] = user_data
        session['auth_method'] = "dummy"
        session.modified = True
        logger.info("Dummy user validated and stored in session.")
        return user_data, None

    def sync_users_to_db(self) -> None:
        """No-op sync in dummy provider."""
        logger.info("Dummy sync - skipping user sync.")
