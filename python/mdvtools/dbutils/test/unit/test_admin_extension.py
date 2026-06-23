from __future__ import annotations

import sys
from dataclasses import dataclass
from types import ModuleType
from typing import Any

import pytest
from flask import Flask

from mdvtools.dbutils.admin_contracts import (
    AdminConflictError,
    AdminInputError,
    AdminNotFoundError,
    AdminProject,
    AdminProjectMember,
    AdminProjectMembership,
    AdminUser,
    CreateAdminUserInput,
    CreateAdminUserResult,
    ProjectMemberInput,
)
from mdvtools.dbutils.admin_extension import AdminExtension


@dataclass
class FakeUser:
    id: int
    email: str
    first_name: str = ""
    last_name: str = ""
    is_active: bool = True
    is_admin: bool = False
    administrator: bool = False


@dataclass
class FakeProject:
    id: int
    name: str
    path: str
    access_level: str = "editable"
    is_public: bool = False
    is_deleted: bool = False
    update_timestamp: Any = None


class FakeQuery:
    def __init__(self, rows: list[Any]):
        self.rows = rows

    def all(self):
        return list(self.rows)


class FakeAdminServices:
    def __init__(self):
        self.create_user_input: CreateAdminUserInput | None = None
        self.add_member_input: ProjectMemberInput | None = None
        self.updated_member: tuple[int, int, str] | None = None
        self.removed_member: tuple[int, int] | None = None

    def list_users(self):
        return []

    def list_projects(self):
        return [
            AdminProject(
                id=10,
                name="Alpha",
                path="/projects/10",
                access_level="editable",
                is_public=False,
                is_deleted=False,
                updated_at=None,
            )
        ]

    def create_local_user_with_project_access(self, data: CreateAdminUserInput):
        self.create_user_input = data
        return CreateAdminUserResult(
            user=AdminUser(
                id=3,
                email=data.email.strip().lower(),
                first_name=data.first_name,
                last_name=data.last_name,
                is_active=True,
                is_admin=False,
            ),
            project_access=[
                AdminProjectMembership(
                    project_id=access.project_id,
                    permission=access.permission,
                    can_read=True,
                    can_write=access.permission in {"edit", "owner"},
                    is_owner=access.permission == "owner",
                )
                for access in data.project_access
            ],
            created=True,
        )

    def list_project_members(self, project_id: int):
        if project_id != 10:
            raise AdminNotFoundError("Project not found.")
        return [
            AdminProjectMember(
                user=AdminUser(
                    id=3,
                    email="member@example.com",
                    first_name="Project",
                    last_name="Member",
                    is_active=True,
                    is_admin=False,
                ),
                project_access=AdminProjectMembership(
                    project_id=project_id,
                    permission="edit",
                    can_read=True,
                    can_write=True,
                    is_owner=False,
                ),
            )
        ]

    def add_project_member(self, project_id: int, data: ProjectMemberInput):
        self.add_member_input = data
        return AdminProjectMember(
            user=AdminUser(
                id=data.user_id,
                email="added@example.com",
                first_name="Added",
                last_name="Member",
                is_active=True,
                is_admin=False,
            ),
            project_access=AdminProjectMembership(
                project_id=project_id,
                permission=data.permission,
                can_read=True,
                can_write=data.permission in {"edit", "owner"},
                is_owner=data.permission == "owner",
            ),
        )

    def update_project_member_permission(
        self,
        project_id: int,
        user_id: int,
        permission: str,
    ):
        self.updated_member = (project_id, user_id, permission)
        return AdminProjectMember(
            user=AdminUser(
                id=user_id,
                email="updated@example.com",
                first_name="Updated",
                last_name="Member",
                is_active=True,
                is_admin=False,
            ),
            project_access=AdminProjectMembership(
                project_id=project_id,
                permission="owner",
                can_read=True,
                can_write=True,
                is_owner=True,
            ),
        )

    def remove_project_member(self, project_id: int, user_id: int):
        self.removed_member = (project_id, user_id)


@pytest.fixture
def fake_db_modules(monkeypatch):
    users = [
        FakeUser(id=1, email="owner@example.com", first_name="Project", last_name="Owner"),
        FakeUser(id=2, email="viewer@example.com"),
    ]
    projects = [
        FakeProject(id=10, name="Alpha", path="/projects/10"),
        FakeProject(id=20, name="Beta", path="/projects/20"),
    ]

    FakeUser.query = FakeQuery(users)
    FakeProject.query = FakeQuery(projects)

    dbmodels = ModuleType("mdvtools.dbutils.dbmodels")
    dbmodels.User = FakeUser
    dbmodels.Project = FakeProject

    monkeypatch.setitem(sys.modules, "mdvtools.dbutils.dbmodels", dbmodels)

    return {
        "users": users,
        "projects": projects,
    }


@pytest.fixture
def app(fake_db_modules):
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension().register_global_routes(flask_app, {})
    return flask_app


def test_local_dev_can_access_session(app):
    response = app.test_client().get("/admin/api/session")

    assert response.status_code == 200
    body = response.get_json()
    assert body["authEnabled"] is False
    assert body["isAdmin"] is True
    assert body["user"]["synthetic"] is True


def test_auth_enabled_rejects_missing_session(fake_db_modules):
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = True
    AdminExtension().register_global_routes(flask_app, {})

    response = flask_app.test_client().get("/admin/api/session")

    assert response.status_code == 401


def test_auth_enabled_rejects_non_admin_session(fake_db_modules):
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = True
    AdminExtension().register_global_routes(flask_app, {})
    client = flask_app.test_client()
    with client.session_transaction() as session:
        session["user"] = {"id": 1, "email": "user@example.com", "is_admin": False}

    response = client.get("/admin/api/session")

    assert response.status_code == 403


def test_admin_session_can_list_users_and_projects(fake_db_modules):
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = True
    AdminExtension().register_global_routes(flask_app, {})
    client = flask_app.test_client()
    with client.session_transaction() as session:
        session["user"] = {"id": 1, "email": "admin@example.com", "is_admin": True}

    users_response = client.get("/admin/api/users")
    projects_response = client.get("/admin/api/projects")

    assert users_response.status_code == 200
    assert projects_response.status_code == 200
    assert users_response.get_json()["users"][0]["email"] == "owner@example.com"
    assert projects_response.get_json()["projects"][0]["name"] == "Alpha"


def test_users_endpoint_can_return_empty_list(monkeypatch):
    FakeUser.query = FakeQuery([])
    FakeProject.query = FakeQuery([])

    dbmodels = ModuleType("mdvtools.dbutils.dbmodels")
    dbmodels.User = FakeUser
    dbmodels.Project = FakeProject
    monkeypatch.setitem(sys.modules, "mdvtools.dbutils.dbmodels", dbmodels)

    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension().register_global_routes(flask_app, {})

    response = flask_app.test_client().get("/admin/api/users")

    assert response.status_code == 200
    assert response.get_json()["users"] == []


def test_projects_endpoint_can_return_empty_list(monkeypatch):
    FakeUser.query = FakeQuery([])
    FakeProject.query = FakeQuery([])

    dbmodels = ModuleType("mdvtools.dbutils.dbmodels")
    dbmodels.User = FakeUser
    dbmodels.Project = FakeProject
    monkeypatch.setitem(sys.modules, "mdvtools.dbutils.dbmodels", dbmodels)

    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension().register_global_routes(flask_app, {})

    response = flask_app.test_client().get("/admin/api/projects")

    assert response.status_code == 200
    assert response.get_json()["projects"] == []


def test_create_local_user_accepts_project_access_array():
    services = FakeAdminServices()
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=services).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/users",
        json={
            "email": "New.User@Example.com",
            "firstName": "New",
            "lastName": "User",
            "projectAccess": [
                {
                    "projectId": 10,
                    "permission": "edit",
                },
                {
                    "projectId": 20,
                    "permission": "view",
                },
            ],
        },
    )

    assert response.status_code == 201
    body = response.get_json()
    assert body["user"]["email"] == "new.user@example.com"
    assert body["projectAccess"][0]["projectId"] == 10
    assert body["projectAccess"][0]["permission"] == "edit"
    assert body["projectAccess"][1]["projectId"] == 20
    assert body["projectAccess"][1]["permission"] == "view"
    assert services.create_user_input is not None
    assert services.create_user_input.project_access[0].project_id == 10
    assert services.create_user_input.project_access[1].project_id == 20


def test_create_local_user_allows_missing_project_access():
    services = FakeAdminServices()
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=services).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/users",
        json={"email": "new.user@example.com"},
    )

    assert response.status_code == 201
    assert response.get_json()["projectAccess"] == []
    assert services.create_user_input is not None
    assert services.create_user_input.project_access == []


def test_create_local_user_allows_empty_project_access():
    services = FakeAdminServices()
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=services).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/users",
        json={"email": "new.user@example.com", "projectAccess": []},
    )

    assert response.status_code == 201
    assert response.get_json()["projectAccess"] == []


def test_create_local_user_rejects_invalid_permission():
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=FakeAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/users",
        json={
            "email": "new.user@example.com",
            "projectAccess": [
                {
                    "projectId": 10,
                    "permission": "admin",
                },
            ],
        },
    )

    assert response.status_code == 400
    assert response.get_json()["error"] == "Permission must be one of: view, edit, owner."


def test_create_local_user_rejects_duplicate_project_access_entries():
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=FakeAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/users",
        json={
            "email": "new.user@example.com",
            "projectAccess": [
                {
                    "projectId": 10,
                    "permission": "view",
                },
                {
                    "projectId": 10,
                    "permission": "edit",
                },
            ],
        },
    )

    assert response.status_code == 400
    assert response.get_json()["error"] == "Project access entries must not contain duplicate projects."


def test_create_local_user_maps_service_validation_errors():
    class RejectingAdminServices(FakeAdminServices):
        def create_local_user_with_project_access(self, data: CreateAdminUserInput):
            raise AdminInputError("A valid email address is required.")

    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=RejectingAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/users",
        json={
            "email": "not-an-email",
            "projectAccess": [
                {
                    "projectId": 10,
                    "permission": "view",
                },
            ],
        },
    )

    assert response.status_code == 400
    assert response.get_json()["error"] == "A valid email address is required."


def test_create_local_user_maps_missing_project_errors():
    class MissingProjectAdminServices(FakeAdminServices):
        def create_local_user_with_project_access(self, data: CreateAdminUserInput):
            raise AdminNotFoundError("Project not found.")

    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=MissingProjectAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/users",
        json={
            "email": "new.user@example.com",
            "projectAccess": [
                {
                    "projectId": 999,
                    "permission": "view",
                },
            ],
        },
    )

    assert response.status_code == 404
    assert response.get_json()["error"] == "Project not found."


def test_project_members_endpoint_lists_assigned_users():
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=FakeAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().get("/admin/api/projects/10/users")

    assert response.status_code == 200
    body = response.get_json()
    assert body["members"][0]["user"]["email"] == "member@example.com"
    assert body["members"][0]["projectAccess"]["permission"] == "edit"
    assert body["members"][0]["projectAccess"]["canWrite"] is True


def test_project_members_endpoint_maps_missing_project_errors():
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=FakeAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().get("/admin/api/projects/999/users")

    assert response.status_code == 404
    assert response.get_json()["error"] == "Project not found."


def test_add_existing_user_to_project():
    services = FakeAdminServices()
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=services).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/projects/10/users",
        json={"userId": 3, "permission": "view"},
    )

    assert response.status_code == 201
    assert response.get_json()["projectAccess"]["permission"] == "view"
    assert services.add_member_input is not None
    assert services.add_member_input.user_id == 3


def test_add_existing_user_to_project_rejects_invalid_permission():
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=FakeAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/projects/10/users",
        json={"userId": 3, "permission": "admin"},
    )

    assert response.status_code == 400
    assert response.get_json()["error"] == "Permission must be one of: view, edit, owner."


def test_update_project_member_permission():
    services = FakeAdminServices()
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=services).register_global_routes(flask_app, {})

    response = flask_app.test_client().patch(
        "/admin/api/projects/10/users/3",
        json={"permission": "owner"},
    )

    assert response.status_code == 200
    assert response.get_json()["projectAccess"]["permission"] == "owner"
    assert services.updated_member == (10, 3, "owner")


def test_update_project_member_permission_maps_last_owner_errors():
    class LastOwnerAdminServices(FakeAdminServices):
        def update_project_member_permission(
            self,
            project_id: int,
            user_id: int,
            permission: str,
        ):
            raise AdminConflictError("Cannot remove or demote the last project owner.")

    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=LastOwnerAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().patch(
        "/admin/api/projects/10/users/3",
        json={"permission": "edit"},
    )

    assert response.status_code == 409
    assert response.get_json()["error"] == "Cannot remove or demote the last project owner."


def test_remove_project_member():
    services = FakeAdminServices()
    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=services).register_global_routes(flask_app, {})

    response = flask_app.test_client().delete("/admin/api/projects/10/users/3")

    assert response.status_code == 204
    assert services.removed_member == (10, 3)


def test_remove_project_member_maps_last_owner_errors():
    class LastOwnerAdminServices(FakeAdminServices):
        def remove_project_member(self, project_id: int, user_id: int):
            raise AdminConflictError("Cannot remove or demote the last project owner.")

    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=LastOwnerAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().delete("/admin/api/projects/10/users/3")

    assert response.status_code == 409
    assert response.get_json()["error"] == "Cannot remove or demote the last project owner."
