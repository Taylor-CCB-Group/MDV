from __future__ import annotations

from typing import Any

from flask import Flask, Response, jsonify, render_template, request, session

from mdvtools.dbutils.admin_contracts import (
    AdminConflictError,
    AdminHostServices,
    AdminInputError,
    AdminNotFoundError,
    AdminPermission,
    CreateAdminUserInput,
    ProjectAccessInput,
    ProjectMemberInput,
    REQUIRE_INITIAL_PROJECT_ACCESS,
)
from mdvtools.dbutils.admin_services import MDVAdminServices
from mdvtools.logging_config import get_logger
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprintProtocol
from mdvtools.server_extension import MDVProjectServerExtension

logger = get_logger(__name__)

ADMIN_PERMISSIONS = [
    "admin:access",
    "admin:users:manage",
    "admin:projects:manage",
]


def _json_error(message: str, status: int):
    return jsonify({"error": message}), status


def _parse_permission(value: Any) -> AdminPermission | None:
    if value == "view":
        return "view"
    if value == "edit":
        return "edit"
    if value == "owner":
        return "owner"
    return None


def _parse_create_user_input(payload: Any) -> CreateAdminUserInput:
    if not isinstance(payload, dict):
        raise AdminInputError("JSON request body is required.")

    email = payload.get("email")
    if not isinstance(email, str):
        raise AdminInputError("Email is required.")

    first_name = payload.get("firstName", "")
    last_name = payload.get("lastName", "")
    if not isinstance(first_name, str) or not isinstance(last_name, str):
        raise AdminInputError("First name and last name must be strings.")

    project_access_payload = payload.get("projectAccess", [])
    if not isinstance(project_access_payload, list):
        raise AdminInputError("Project access must be a list.")
    if REQUIRE_INITIAL_PROJECT_ACCESS and len(project_access_payload) == 0:
        raise AdminInputError("At least one project access assignment is required.")

    seen_project_ids: set[int] = set()
    project_access: list[ProjectAccessInput] = []
    for item in project_access_payload:
        if not isinstance(item, dict):
            raise AdminInputError("Project access entries must be objects.")
        project_id = item.get("projectId")
        if not isinstance(project_id, int):
            raise AdminInputError("Project ID is required.")
        if project_id in seen_project_ids:
            raise AdminInputError("Project access entries must not contain duplicate projects.")
        seen_project_ids.add(project_id)

        permission = _parse_permission(item.get("permission"))
        if permission is None:
            raise AdminInputError("Permission must be one of: view, edit, owner.")
        project_access.append(
            ProjectAccessInput(
                project_id=project_id,
                permission=permission,
            )
        )

    return CreateAdminUserInput(
        email=email,
        first_name=first_name,
        last_name=last_name,
        project_access=project_access,
    )


def _parse_project_member_input(payload: Any) -> ProjectMemberInput:
    if not isinstance(payload, dict):
        raise AdminInputError("JSON request body is required.")

    user_id = payload.get("userId")
    if not isinstance(user_id, int):
        raise AdminInputError("User ID is required.")

    permission = _parse_permission(payload.get("permission"))
    if permission is None:
        raise AdminInputError("Permission must be one of: view, edit, owner.")

    return ProjectMemberInput(user_id=user_id, permission=permission)


def _parse_permission_input(payload: Any) -> AdminPermission:
    if not isinstance(payload, dict):
        raise AdminInputError("JSON request body is required.")

    permission = _parse_permission(payload.get("permission"))
    if permission is None:
        raise AdminInputError("Permission must be one of: view, edit, owner.")
    return permission


def _safe_session_user(user: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": user.get("id"),
        "email": user.get("email"),
        "is_admin": bool(user.get("is_admin")),
        "synthetic": bool(user.get("synthetic", False)),
    }


class AdminExtension(MDVProjectServerExtension):
    extension_id = "admin"

    def __init__(self, services: AdminHostServices | None = None):
        self.services = services or MDVAdminServices()

    def get_session_config(self) -> dict[str, Any]:
        return {
            "extensionId": self.extension_id,
            "permissions": ADMIN_PERMISSIONS,
        }

    def register_global_routes(self, app: Flask, config: dict):
        enable_auth = bool(app.config.get("ENABLE_AUTH", False))

        def require_admin() -> tuple[dict[str, Any], None] | tuple[None, tuple[Response, int]]:
            if not enable_auth:
                return {
                    "id": 0,
                    "email": "local-admin@localhost",
                    "is_admin": True,
                    "synthetic": True,
                }, None

            user = session.get("user")
            if not user:
                return None, _json_error("Admin session required", 401)
            if not bool(user.get("is_admin")):
                return None, _json_error("Admin access required", 403)
            return user, None

        @app.route("/admin", methods=["GET"])
        @app.route("/admin/", methods=["GET"])
        def admin_page():
            return render_template("admin_build.html")

        @app.route("/admin/api/session", methods=["GET"])
        def admin_session():
            user, error = require_admin()
            if error is not None:
                return error
            return jsonify({
                "authEnabled": enable_auth,
                "user": _safe_session_user(user),
                "isAdmin": True,
                "permissions": ADMIN_PERMISSIONS,
            })

        @app.route("/admin/api/users", methods=["GET"])
        def admin_users():
            _user, error = require_admin()
            if error is not None:
                return error
            users = self.services.list_users()
            return jsonify({"users": [user.to_response() for user in users]})

        @app.route("/admin/api/users", methods=["POST"])
        def admin_create_user():
            _user, error = require_admin()
            if error is not None:
                return error
            try:
                data = _parse_create_user_input(request.get_json(silent=True))
                result = self.services.create_local_user_with_project_access(data)
                return jsonify(result.to_response()), 201 if result.created else 200
            except AdminInputError as exc:
                return _json_error(str(exc), 400)
            except AdminConflictError as exc:
                return _json_error(str(exc), 409)
            except AdminNotFoundError as exc:
                return _json_error(str(exc), 404)

        @app.route("/admin/api/projects", methods=["GET"])
        def admin_projects():
            _user, error = require_admin()
            if error is not None:
                return error
            projects = self.services.list_projects()
            return jsonify({"projects": [project.to_response() for project in projects]})

        @app.route("/admin/api/projects/<int:project_id>/users", methods=["GET"])
        def admin_project_users(project_id: int):
            _user, error = require_admin()
            if error is not None:
                return error
            try:
                members = self.services.list_project_members(project_id)
                return jsonify({"members": [member.to_response() for member in members]})
            except AdminNotFoundError as exc:
                return _json_error(str(exc), 404)

        @app.route("/admin/api/projects/<int:project_id>/users", methods=["POST"])
        def admin_add_project_user(project_id: int):
            _user, error = require_admin()
            if error is not None:
                return error
            try:
                data = _parse_project_member_input(request.get_json(silent=True))
                member = self.services.add_project_member(project_id, data)
                return jsonify(member.to_response()), 201
            except AdminInputError as exc:
                return _json_error(str(exc), 400)
            except AdminConflictError as exc:
                return _json_error(str(exc), 409)
            except AdminNotFoundError as exc:
                return _json_error(str(exc), 404)

        @app.route("/admin/api/projects/<int:project_id>/users/<int:user_id>", methods=["PATCH"])
        def admin_update_project_user(project_id: int, user_id: int):
            _user, error = require_admin()
            if error is not None:
                return error
            try:
                permission = _parse_permission_input(request.get_json(silent=True))
                member = self.services.update_project_member_permission(
                    project_id,
                    user_id,
                    permission,
                )
                return jsonify(member.to_response())
            except AdminInputError as exc:
                return _json_error(str(exc), 400)
            except AdminConflictError as exc:
                return _json_error(str(exc), 409)
            except AdminNotFoundError as exc:
                return _json_error(str(exc), 404)

        @app.route("/admin/api/projects/<int:project_id>/users/<int:user_id>", methods=["DELETE"])
        def admin_delete_project_user(project_id: int, user_id: int):
            _user, error = require_admin()
            if error is not None:
                return error
            try:
                self.services.remove_project_member(project_id, user_id)
                return "", 204
            except AdminInputError as exc:
                return _json_error(str(exc), 400)
            except AdminConflictError as exc:
                return _json_error(str(exc), 409)
            except AdminNotFoundError as exc:
                return _json_error(str(exc), 404)

        logger.info("AdminExtension registered /admin and /admin/api routes")

    def register_routes(self, project: MDVProject, project_bp: ProjectBlueprintProtocol):
        return None
