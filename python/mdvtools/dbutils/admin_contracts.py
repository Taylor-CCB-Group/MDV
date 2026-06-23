from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Protocol


"""
Admin-specific contracts for the Admin Portal extension.

These contracts are intentionally admin-domain specific. They are not the final
generic MDV plugin SDK. The goal is to give the Admin extension a narrow boundary
now, so the same shape can later sit behind a broader plugin host API.
"""


AdminPermission = Literal["view", "edit", "owner"]
REQUIRE_INITIAL_PROJECT_ACCESS = False


class AdminServiceError(Exception):
    pass


class AdminInputError(AdminServiceError):
    pass


class AdminConflictError(AdminServiceError):
    pass


class AdminNotFoundError(AdminServiceError):
    pass


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


@dataclass(frozen=True)
class ProjectAccessInput:
    project_id: int
    permission: AdminPermission


@dataclass(frozen=True)
class CreateAdminUserInput:
    email: str
    first_name: str
    last_name: str
    project_access: list[ProjectAccessInput]


@dataclass(frozen=True)
class AdminProjectMembership:
    project_id: int
    permission: AdminPermission
    can_read: bool
    can_write: bool
    is_owner: bool

    def to_response(self) -> dict[str, Any]:
        return {
            "projectId": self.project_id,
            "permission": self.permission,
            "canRead": self.can_read,
            "canWrite": self.can_write,
            "isOwner": self.is_owner,
        }


@dataclass(frozen=True)
class CreateAdminUserResult:
    user: AdminUser
    project_access: list[AdminProjectMembership]
    created: bool

    def to_response(self) -> dict[str, Any]:
        return {
            "user": self.user.to_response(),
            "projectAccess": [access.to_response() for access in self.project_access],
            "created": self.created,
        }


@dataclass(frozen=True)
class ProjectMemberInput:
    user_id: int
    permission: AdminPermission


@dataclass(frozen=True)
class AdminProjectMember:
    user: AdminUser
    project_access: AdminProjectMembership

    def to_response(self) -> dict[str, Any]:
        return {
            "user": self.user.to_response(),
            "projectAccess": self.project_access.to_response(),
        }


class AdminHostServices(Protocol):
    def list_users(self) -> list[AdminUser]:
        ...

    def list_projects(self) -> list[AdminProject]:
        ...

    def create_local_user_with_project_access(
        self,
        data: CreateAdminUserInput,
    ) -> CreateAdminUserResult:
        ...

    def list_project_members(self, project_id: int) -> list[AdminProjectMember]:
        ...

    def add_project_member(self, project_id: int, data: ProjectMemberInput) -> AdminProjectMember:
        ...

    def update_project_member_permission(
        self,
        project_id: int,
        user_id: int,
        permission: AdminPermission,
    ) -> AdminProjectMember:
        ...

    def remove_project_member(self, project_id: int, user_id: int) -> None:
        ...
