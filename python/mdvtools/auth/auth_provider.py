from abc import ABC, abstractmethod
from flask import Response
from typing import Optional, Tuple, Union

class AuthProvider(ABC):
    @abstractmethod
    def login(self) -> Union[str, Response, Tuple[Response, int]]:
        """Redirects to the login page of the provider."""
        pass

    @abstractmethod
    def logout(self) -> Union[str, Response, Tuple[Response, int]]:
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
    def validate_user(self) -> Tuple[Optional[dict], Optional[Tuple]]:
        """Validates and returns user data."""
        pass

    @abstractmethod
    def sync_users_to_db(self):
        """ Syncs users from the authentication provider (e.g., Auth0) to the application's database."""
        pass
