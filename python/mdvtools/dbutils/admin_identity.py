from __future__ import annotations

import secrets
import string
from dataclasses import dataclass
from typing import Any, Mapping, Protocol

from mdvtools.dbutils.admin_contracts import AdminExternalServiceError, AdminInputError


@dataclass(frozen=True)
class AdminIdentityInput:
    email: str
    first_name: str
    last_name: str


@dataclass(frozen=True)
class AdminIdentityResult:
    email: str
    auth_id: str
    created: bool


class AdminIdentityProvider(Protocol):
    def create_or_resolve_user(self, data: AdminIdentityInput) -> AdminIdentityResult:
        ...

    def rollback_created_user(self, auth_id: str) -> None:
        ...


class LocalAdminIdentityProvider:
    def create_or_resolve_user(self, data: AdminIdentityInput) -> AdminIdentityResult:
        email = _normalise_email(data.email)
        return AdminIdentityResult(
            email=email,
            auth_id=f"local:{email}",
            created=False,
        )

    def rollback_created_user(self, auth_id: str) -> None:
        return None


class ConfiguredAuth0AdminIdentityProvider:
    def __init__(self, config: Mapping[str, Any]):
        self.config = config

    def create_or_resolve_user(self, data: AdminIdentityInput) -> AdminIdentityResult:
        return Auth0AdminIdentityProvider.from_config(self.config).create_or_resolve_user(data)

    def rollback_created_user(self, auth_id: str) -> None:
        return Auth0AdminIdentityProvider.from_config(self.config).rollback_created_user(auth_id)


class Auth0AdminIdentityProvider:
    def __init__(
        self,
        domain: str,
        client_id: str,
        client_secret: str,
        connection: str,
        password_length: int = 16,
    ):
        self.domain = _require_config("AUTH0_DOMAIN", domain)
        self.client_id = _require_config("AUTH0_CLIENT_ID", client_id)
        self.client_secret = _require_config("AUTH0_CLIENT_SECRET", client_secret)
        self.connection = _require_config("AUTH0_DB_CONNECTION", connection)
        self.password_length = password_length

    @classmethod
    def from_config(cls, config: Mapping[str, Any]) -> Auth0AdminIdentityProvider:
        return cls(
            domain=str(config.get("AUTH0_DOMAIN") or ""),
            client_id=str(config.get("AUTH0_CLIENT_ID") or ""),
            client_secret=str(config.get("AUTH0_CLIENT_SECRET") or ""),
            connection=str(config.get("AUTH0_DB_CONNECTION") or ""),
        )

    def create_or_resolve_user(self, data: AdminIdentityInput) -> AdminIdentityResult:
        email = _normalise_email(data.email)
        token = self._get_management_api_token()
        existing_users = self._get_users_by_email(token, email)
        existing_auth_id = self._find_user_id_in_connection(existing_users)
        if existing_auth_id:
            return AdminIdentityResult(email=email, auth_id=existing_auth_id, created=False)

        created_user = self._create_user(token, email)
        auth_id = created_user.get("user_id")
        if not isinstance(auth_id, str) or not auth_id:
            raise AdminExternalServiceError("Auth0 did not return a user ID.")

        return AdminIdentityResult(email=email, auth_id=auth_id, created=True)

    def rollback_created_user(self, auth_id: str) -> None:
        if not auth_id:
            return None
        token = self._get_management_api_token()
        self._delete_user(token, auth_id)

    def _get_management_api_token(self) -> str:
        response = self._request(
            "post",
            f"https://{self.domain}/oauth/token",
            json={
                "grant_type": "client_credentials",
                "client_id": self.client_id,
                "client_secret": self.client_secret,
                "audience": f"https://{self.domain}/api/v2/",
            },
            headers={"Content-Type": "application/json"},
        )
        token = response.get("access_token")
        if not isinstance(token, str) or not token:
            raise AdminExternalServiceError("Auth0 did not return a Management API token.")
        return token

    def _get_users_by_email(self, token: str, email: str) -> list[dict[str, Any]]:
        response = self._request(
            "get",
            f"https://{self.domain}/api/v2/users-by-email",
            headers={
                "Authorization": f"Bearer {token}",
                "Content-Type": "application/json",
            },
            params={"email": email},
        )
        if not isinstance(response, list):
            raise AdminExternalServiceError("Auth0 users-by-email response was not a list.")
        return [user for user in response if isinstance(user, dict)]

    def _create_user(self, token: str, email: str) -> dict[str, Any]:
        response = self._request(
            "post",
            f"https://{self.domain}/api/v2/users",
            headers={
                "Authorization": f"Bearer {token}",
                "Content-Type": "application/json",
            },
            json={
                "connection": self.connection,
                "email": email,
                "password": _generate_random_password(self.password_length),
                "email_verified": True,
            },
        )
        if not isinstance(response, dict):
            raise AdminExternalServiceError("Auth0 create user response was not an object.")
        return response

    def _delete_user(self, token: str, auth_id: str) -> None:
        self._request(
            "delete",
            f"https://{self.domain}/api/v2/users/{auth_id}",
            headers={
                "Authorization": f"Bearer {token}",
                "Content-Type": "application/json",
            },
        )

    def _find_user_id_in_connection(self, users: list[dict[str, Any]]) -> str | None:
        for user in users:
            identities = user.get("identities")
            if not isinstance(identities, list):
                continue
            for identity in identities:
                if not isinstance(identity, dict):
                    continue
                if identity.get("connection") != self.connection:
                    continue
                user_id = user.get("user_id")
                if isinstance(user_id, str) and user_id:
                    return user_id
        return None

    def _request(self, method: str, url: str, **kwargs: Any) -> Any:
        import requests

        try:
            if method == "post":
                response = requests.post(url, **kwargs)
            elif method == "get":
                response = requests.get(url, **kwargs)
            elif method == "delete":
                response = requests.delete(url, **kwargs)
            else:
                raise AdminExternalServiceError(f"Unsupported Auth0 request method: {method}.")
            response.raise_for_status()
            if method == "delete":
                return {}
            return response.json()
        except requests.RequestException as exc:
            raise AdminExternalServiceError("Auth0 operation failed.") from exc


def _normalise_email(email: str) -> str:
    value = email.strip().lower()
    if not value:
        raise AdminInputError("Email is required.")
    if "@" not in value:
        raise AdminInputError("A valid email address is required.")
    return value


def _require_config(name: str, value: str) -> str:
    stripped = value.strip()
    if not stripped:
        raise AdminExternalServiceError(f"{name} is required for Auth0 admin user creation.")
    return stripped


def _generate_random_password(length: int) -> str:
    alphabet = string.ascii_letters + string.digits + "!@#$%^&*()-_=+"
    return "".join(secrets.choice(alphabet) for _ in range(length))
