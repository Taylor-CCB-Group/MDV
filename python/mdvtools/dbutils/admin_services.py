from __future__ import annotations

from datetime import datetime
from typing import Any

from mdvtools.dbutils.admin_contracts import (
    AdminConflictError,
    AdminInputError,
    AdminNotFoundError,
    AdminPermission,
    AdminProject,
    AdminProjectMember,
    AdminProjectMembership,
    AdminUser,
    CreateAdminUserInput,
    CreateAdminUserResult,
    ProjectMemberInput,
)


VALID_PERMISSIONS: set[AdminPermission] = {"view", "edit", "owner"}


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

    def create_local_user_with_project_access(
        self,
        data: CreateAdminUserInput,
    ) -> CreateAdminUserResult:
        email = data.email.strip().lower()
        if not email:
            raise AdminInputError("Email is required.")
        if "@" not in email:
            raise AdminInputError("A valid email address is required.")
        for access in data.project_access:
            if access.permission not in VALID_PERMISSIONS:
                raise AdminInputError("Permission must be one of: view, edit, owner.")

        Project, User, UserProject, db = self._get_write_models()
        try:
            user = User.query.filter_by(email=email).first()
            created = user is None
            if user is None:
                user = User()
                user.email = email
                user.auth_id = f"local:{email}"
                user.first_name = data.first_name.strip()
                user.last_name = data.last_name.strip()
                user.confirmed_at = datetime.now()
                user.is_active = True
                user.password = ""
                user.administrator = False
                user.is_admin = False
                db.session.add(user)
                db.session.flush()
            else:
                first_name = data.first_name.strip()
                last_name = data.last_name.strip()
                if first_name:
                    user.first_name = first_name
                if last_name:
                    user.last_name = last_name

            project_memberships: list[AdminProjectMembership] = []
            for access in data.project_access:
                project = Project.query.get(access.project_id)
                if project is None:
                    raise AdminNotFoundError("Project not found.")

                permission = self._permission_flags(access.permission)
                user_project = UserProject.query.filter_by(
                    user_id=user.id,
                    project_id=project.id,
                ).first()
                if user_project is None:
                    user_project = UserProject(
                        user_id=user.id,
                        project_id=project.id,
                        is_owner=permission["is_owner"],
                        can_read=permission["can_read"],
                        can_write=permission["can_write"],
                    )
                    db.session.add(user_project)
                else:
                    user_project.can_read = permission["can_read"]
                    user_project.can_write = permission["can_write"]
                    user_project.is_owner = permission["is_owner"]
                project_memberships.append(
                    AdminProjectMembership(
                        project_id=int(project.id),
                        permission=access.permission,
                        can_read=permission["can_read"],
                        can_write=permission["can_write"],
                        is_owner=permission["is_owner"],
                    )
                )

            db.session.commit()
            return CreateAdminUserResult(
                user=self._to_admin_user(user),
                project_access=project_memberships,
                created=created,
            )
        except (AdminInputError, AdminNotFoundError, AdminConflictError):
            db.session.rollback()
            raise
        except Exception:
            db.session.rollback()
            raise

    def list_project_members(self, project_id: int) -> list[AdminProjectMember]:
        Project, User, UserProject = self._get_membership_models()
        project = Project.query.get(project_id)
        if project is None:
            raise AdminNotFoundError("Project not found.")

        memberships = UserProject.query.filter_by(project_id=project_id).all()
        members: list[AdminProjectMember] = []
        for membership in memberships:
            user = User.query.get(membership.user_id)
            if user is None:
                continue
            permission = self._permission_from_membership(membership)
            members.append(
                AdminProjectMember(
                    user=self._to_admin_user(user),
                    project_access=AdminProjectMembership(
                        project_id=project_id,
                        permission=permission,
                        can_read=bool(membership.can_read),
                        can_write=bool(membership.can_write),
                        is_owner=bool(membership.is_owner),
                    ),
                )
            )
        return members

    def add_project_member(self, project_id: int, data: ProjectMemberInput) -> AdminProjectMember:
        if data.permission not in VALID_PERMISSIONS:
            raise AdminInputError("Permission must be one of: view, edit, owner.")

        Project, User, UserProject, db = self._get_write_models()
        try:
            project = Project.query.get(project_id)
            if project is None:
                raise AdminNotFoundError("Project not found.")
            user = User.query.get(data.user_id)
            if user is None:
                raise AdminNotFoundError("User not found.")

            membership = UserProject.query.filter_by(
                user_id=data.user_id,
                project_id=project_id,
            ).first()
            permission = self._permission_flags(data.permission)
            if membership is None:
                membership = UserProject(
                    user_id=data.user_id,
                    project_id=project_id,
                    is_owner=permission["is_owner"],
                    can_read=permission["can_read"],
                    can_write=permission["can_write"],
                )
                db.session.add(membership)
            else:
                membership.can_read = permission["can_read"]
                membership.can_write = permission["can_write"]
                membership.is_owner = permission["is_owner"]

            db.session.commit()
            return self._to_admin_project_member(user, membership)
        except (AdminInputError, AdminNotFoundError, AdminConflictError):
            db.session.rollback()
            raise
        except Exception:
            db.session.rollback()
            raise

    def update_project_member_permission(
        self,
        project_id: int,
        user_id: int,
        permission: AdminPermission,
    ) -> AdminProjectMember:
        if permission not in VALID_PERMISSIONS:
            raise AdminInputError("Permission must be one of: view, edit, owner.")

        Project, User, UserProject, db = self._get_write_models()
        try:
            project = Project.query.get(project_id)
            if project is None:
                raise AdminNotFoundError("Project not found.")
            user = User.query.get(user_id)
            if user is None:
                raise AdminNotFoundError("User not found.")
            membership = UserProject.query.filter_by(
                user_id=user_id,
                project_id=project_id,
            ).first()
            if membership is None:
                raise AdminNotFoundError("Project member not found.")
            if bool(membership.is_owner) and permission != "owner":
                self._ensure_not_last_owner(project_id, excluded_user_id=user_id)

            permission_flags = self._permission_flags(permission)
            membership.can_read = permission_flags["can_read"]
            membership.can_write = permission_flags["can_write"]
            membership.is_owner = permission_flags["is_owner"]

            db.session.commit()
            return self._to_admin_project_member(user, membership)
        except (AdminInputError, AdminNotFoundError, AdminConflictError):
            db.session.rollback()
            raise
        except Exception:
            db.session.rollback()
            raise

    def remove_project_member(self, project_id: int, user_id: int) -> None:
        Project, User, UserProject, db = self._get_write_models()
        try:
            project = Project.query.get(project_id)
            if project is None:
                raise AdminNotFoundError("Project not found.")
            user = User.query.get(user_id)
            if user is None:
                raise AdminNotFoundError("User not found.")
            membership = UserProject.query.filter_by(
                user_id=user_id,
                project_id=project_id,
            ).first()
            if membership is None:
                raise AdminNotFoundError("Project member not found.")
            if bool(membership.is_owner):
                self._ensure_not_last_owner(project_id, excluded_user_id=user_id)

            db.session.delete(membership)
            db.session.commit()
        except (AdminInputError, AdminNotFoundError, AdminConflictError):
            db.session.rollback()
            raise
        except Exception:
            db.session.rollback()
            raise

    def _get_models(self):
        from mdvtools.dbutils.dbmodels import Project, User

        return Project, User

    def _get_write_models(self):
        from mdvtools.dbutils.dbmodels import Project, User, UserProject, db

        return Project, User, UserProject, db

    def _get_membership_models(self):
        from mdvtools.dbutils.dbmodels import Project, User, UserProject

        return Project, User, UserProject

    def _permission_flags(self, permission: AdminPermission) -> dict[str, bool]:
        return {
            "can_read": True,
            "can_write": permission in {"edit", "owner"},
            "is_owner": permission == "owner",
        }

    def _permission_from_membership(self, membership: Any) -> AdminPermission:
        if bool(membership.is_owner):
            return "owner"
        if bool(membership.can_write):
            return "edit"
        return "view"

    def _ensure_not_last_owner(self, project_id: int, excluded_user_id: int) -> None:
        _Project, _User, UserProject = self._get_membership_models()
        owner_count = UserProject.query.filter_by(
            project_id=project_id,
            is_owner=True,
        ).count()
        excluded_owner = UserProject.query.filter_by(
            project_id=project_id,
            user_id=excluded_user_id,
            is_owner=True,
        ).first()
        if excluded_owner is not None and owner_count <= 1:
            raise AdminConflictError("Cannot remove or demote the last project owner.")

    def _to_admin_project_member(self, user: Any, membership: Any) -> AdminProjectMember:
        permission = self._permission_from_membership(membership)
        return AdminProjectMember(
            user=self._to_admin_user(user),
            project_access=AdminProjectMembership(
                project_id=int(membership.project_id),
                permission=permission,
                can_read=bool(membership.can_read),
                can_write=bool(membership.can_write),
                is_owner=bool(membership.is_owner),
            ),
        )

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
