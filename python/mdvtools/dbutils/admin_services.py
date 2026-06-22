from __future__ import annotations

from datetime import datetime
from typing import Any

from mdvtools.dbutils.admin_contracts import AdminProject, AdminUser


def _serialize_datetime(value: Any) -> str | None:
    if isinstance(value, datetime):
        return value.isoformat()
    return None


class MDVAdminServices:
    """
    MDV-backed implementation of the Admin Portal service boundary.

    This class is allowed to depend on current MDV internals such as SQLAlchemy
    models. Admin route handlers should depend on this service boundary instead
    of querying MDV models directly.
    """

    def list_users(self) -> list[AdminUser]:
        _Project, User = self._get_models()
        return [self._to_admin_user(user) for user in User.query.all()]

    def list_projects(self) -> list[AdminProject]:
        Project, _User = self._get_models()
        return [self._to_admin_project(project) for project in Project.query.all()]

    def _get_models(self):
        from mdvtools.dbutils.dbmodels import Project, User

        return Project, User

    def _to_admin_user(self, user: Any) -> AdminUser:
        return AdminUser(
            id=int(user.id),
            email=getattr(user, "email", ""),
            first_name=getattr(user, "first_name", ""),
            last_name=getattr(user, "last_name", ""),
            is_active=bool(getattr(user, "is_active", False)),
            is_admin=bool(getattr(user, "is_admin", False) or getattr(user, "administrator", False)),
        )

    def _to_admin_project(self, project: Any) -> AdminProject:
        return AdminProject(
            id=int(project.id),
            name=getattr(project, "name", f"Project {project.id}"),
            path=getattr(project, "path", ""),
            access_level=getattr(project, "access_level", ""),
            is_public=bool(getattr(project, "is_public", False)),
            is_deleted=bool(getattr(project, "is_deleted", False)),
            updated_at=_serialize_datetime(getattr(project, "update_timestamp", None)),
        )
