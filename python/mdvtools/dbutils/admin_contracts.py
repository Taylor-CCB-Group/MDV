from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol


"""
Admin-specific contracts for the Admin Portal extension.

These contracts are intentionally admin-domain specific. They are not the final
generic MDV plugin SDK. The goal is to give the Admin extension a narrow boundary
now, so the same shape can later sit behind a broader plugin host API.
"""


@dataclass(frozen=True)
class AdminUser:
    id: int
    email: str
    first_name: str
    last_name: str
    is_active: bool
    is_admin: bool

    def to_response(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "email": self.email,
            "firstName": self.first_name,
            "lastName": self.last_name,
            "isActive": self.is_active,
            "isAdmin": self.is_admin,
        }


@dataclass(frozen=True)
class AdminProject:
    id: int
    name: str
    path: str
    access_level: str
    is_public: bool
    is_deleted: bool
    updated_at: str | None

    def to_response(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "path": self.path,
            "accessLevel": self.access_level,
            "isPublic": self.is_public,
            "isDeleted": self.is_deleted,
            "updatedAt": self.updated_at,
        }


class AdminHostServices(Protocol):
    def list_users(self) -> list[AdminUser]:
        ...

    def list_projects(self) -> list[AdminProject]:
        ...
