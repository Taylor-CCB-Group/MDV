from __future__ import annotations

from datetime import datetime
from typing import Any

from flask import Flask, Response, jsonify, render_template, session

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


def _serialize_datetime(value: Any) -> str | None:
    if isinstance(value, datetime):
        return value.isoformat()
    return None


def _get_models():
    from mdvtools.dbutils.dbmodels import Project, User

    return Project, User


def _query_all(model: Any):
    return model.query.all()


def _serialize_user(user: Any) -> dict[str, Any]:
    return {
        "id": int(user.id),
        "email": getattr(user, "email", ""),
        "firstName": getattr(user, "first_name", ""),
        "lastName": getattr(user, "last_name", ""),
        "isActive": bool(getattr(user, "is_active", False)),
        "isAdmin": bool(getattr(user, "is_admin", False) or getattr(user, "administrator", False)),
    }


def _serialize_project(project: Any) -> dict[str, Any]:
    return {
        "id": int(project.id),
        "name": getattr(project, "name", f"Project {project.id}"),
        "path": getattr(project, "path", ""),
        "accessLevel": getattr(project, "access_level", ""),
        "isPublic": bool(getattr(project, "is_public", False)),
        "isDeleted": bool(getattr(project, "is_deleted", False)),
        "updatedAt": _serialize_datetime(getattr(project, "update_timestamp", None)),
    }


class AdminExtension(MDVProjectServerExtension):
    extension_id = "admin"

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
                "user": user,
                "isAdmin": True,
                "permissions": ADMIN_PERMISSIONS,
            })

        @app.route("/admin/api/users", methods=["GET"])
        def admin_users():
            _user, error = require_admin()
            if error is not None:
                return error
            _Project, User = _get_models()
            users = _query_all(User)
            return jsonify({"users": [_serialize_user(user) for user in users]})

        @app.route("/admin/api/projects", methods=["GET"])
        def admin_projects():
            _user, error = require_admin()
            if error is not None:
                return error
            Project, _User = _get_models()
            projects = _query_all(Project)
            return jsonify({"projects": [_serialize_project(project) for project in projects]})

        logger.info("AdminExtension registered /admin and /admin/api routes")

    def register_routes(self, project: MDVProject, project_bp: ProjectBlueprintProtocol):
        return None
