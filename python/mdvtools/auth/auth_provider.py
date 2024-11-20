from abc import ABC, abstractmethod
from typing import Optional

class AuthProvider(ABC):
    @abstractmethod
    def login(self) -> str:
        """Redirects to the login page of the provider."""
        pass

    @abstractmethod
    def logout(self) -> None:
        """Logs out the user."""
        pass

    @abstractmethod
    def get_user(self, token: str) -> Optional[dict]:
        """Fetches the user profile using the provided token."""
        pass

    @abstractmethod
    def get_token(self) -> Optional[str]:
        """Gets the current user's access token."""
        pass

    @abstractmethod
    def handle_callback(self, code: str, redirect_uri: str) -> Optional[str]:
        """Handles the callback and exchanges the code for a token."""
        pass

    @abstractmethod
    def is_authenticated(self, token: str) -> bool:
        """Verifies if the user is authenticated using the token."""
        pass
