from __future__ import annotations

from typing import Any

from flask import Flask, Response, jsonify, render_template, session

from mdvtools.dbutils.admin_contracts import AdminHostServices
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
                "user": user,
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

        @app.route("/admin/api/projects", methods=["GET"])
        def admin_projects():
            _user, error = require_admin()
            if error is not None:
                return error
            projects = self.services.list_projects()
            return jsonify({"projects": [project.to_response() for project in projects]})

        logger.info("AdminExtension registered /admin and /admin/api routes")

    def register_routes(self, project: MDVProject, project_bp: ProjectBlueprintProtocol):
        return None
