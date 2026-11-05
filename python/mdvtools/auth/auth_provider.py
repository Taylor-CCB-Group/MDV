from abc import ABC, abstractmethod
# from flask import Response
from flask.typing import ResponseReturnValue
from dataclasses import dataclass
from typing import Optional, Tuple, Union


@dataclass(frozen=True)
class NoAccess:
    """A person the provider signed in who has no active User on this deployment.

    validate_user returns it in place of an error response, so the auth hook can
    show the no-access page instead of sending them back to sign in, which would
    loop through the provider (ADR 0019).
    """

    email: str
    # The Auth0 connection they signed in through, when the provider knows it.
    sign_in_method: Optional[str] = None


class AuthProvider(ABC):
    @abstractmethod
    def login(self) -> ResponseReturnValue:
        """Redirects to the login page of the provider."""
        pass

    @abstractmethod
    def logout(self) -> ResponseReturnValue:
        """Logs out the user."""
        pass

    @abstractmethod
    def get_user(self, token: Optional[dict] = None) -> Optional[dict]:
        """Fetches the user profile using the provided token."""
        pass

    @abstractmethod
    def get_token(self) -> Optional[str]:
        """Gets the current user's access token."""
        pass

    @abstractmethod
    def handle_callback(self) -> Optional[str]:
        """
        Handles the callback and exchanges the code for a token (implemented by the provider).
        
        :return: Access token if successfully retrieved, else None
        """
        pass

    @abstractmethod
    def validate_user(self) -> Tuple[Optional[dict], Optional[Union[Tuple, NoAccess]]]:
        """Validates and returns user data."""
        pass

    @abstractmethod
    def sync_users_to_db(self):
        """ Syncs users from the authentication provider (e.g., Auth0) to the application's database."""
        pass
