from __future__ import annotations

import sys
from dataclasses import dataclass
from types import ModuleType
from typing import Any

import pytest
from flask import Flask

from mdvtools.dbutils.admin_contracts import (
    AdminConflictError,
    AdminExternalServiceError,
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
from mdvtools.dbutils.admin_identity import AdminIdentityInput, Auth0AdminIdentityProvider
from mdvtools.dbutils.admin_services import MDVAdminServices


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

    def get(self, row_id: int):
        return next((row for row in self.rows if row.id == row_id), None)

    def filter_by(self, **filters: Any):
        return FakeQuery([
            row
            for row in self.rows
            if all(getattr(row, key) == value for key, value in filters.items())
        ])

    def first(self):
        return self.rows[0] if self.rows else None

    def count(self):
        return len(self.rows)


class FakeSession:
    def __init__(self):
        self.added: list[Any] = []
        self.deleted: list[Any] = []
        self.commits = 0
        self.rollbacks = 0

    def add(self, row: Any):
        self.added.append(row)

    def flush(self):
        return None

    def commit(self):
        self.commits += 1

    def rollback(self):
        self.rollbacks += 1

    def delete(self, row: Any):
        self.deleted.append(row)


class FakeDb:
    def __init__(self):
        self.session = FakeSession()


class FakeUserProject:
    query: FakeQuery

    def __init__(
        self,
        user_id: int,
        project_id: int,
        is_owner: bool,
        can_read: bool,
        can_write: bool,
    ):
        self.user_id = user_id
        self.project_id = project_id
        self.is_owner = is_owner
        self.can_read = can_read
        self.can_write = can_write


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

    def create_user_with_project_access(self, data: CreateAdminUserInput):
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
        FakeProject(id=30, name="Deleted", path="/projects/30", is_deleted=True),
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
    assert set(body["user"].keys()) == {"id", "email", "is_admin", "synthetic"}
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
    project_names = [project["name"] for project in projects_response.get_json()["projects"]]
    assert project_names == ["Alpha", "Beta"]


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


def test_admin_project_list_excludes_deleted_projects(fake_db_modules):
    service = MDVAdminServices()

    projects = service.list_projects()

    assert [project.name for project in projects] == ["Alpha", "Beta"]


def test_admin_write_refreshes_auth_cache_when_auth_enabled(monkeypatch):
    FakeUser.query = FakeQuery([FakeUser(id=3, email="member@example.com")])
    FakeProject.query = FakeQuery([FakeProject(id=10, name="Alpha", path="/projects/10")])
    FakeUserProject.query = FakeQuery([])
    fake_db = FakeDb()

    dbmodels = ModuleType("mdvtools.dbutils.dbmodels")
    dbmodels.User = FakeUser
    dbmodels.Project = FakeProject
    dbmodels.UserProject = FakeUserProject
    dbmodels.db = fake_db
    monkeypatch.setitem(sys.modules, "mdvtools.dbutils.dbmodels", dbmodels)

    refresh_calls: list[bool] = []
    service = MDVAdminServices(
        enable_auth=True,
        refresh_auth_cache=refresh_calls.append,
    )

    service.add_project_member(10, ProjectMemberInput(user_id=3, permission="edit"))

    assert fake_db.session.commits == 1
    assert refresh_calls == [True]


def test_admin_write_skips_auth_cache_refresh_when_auth_disabled(monkeypatch):
    FakeUser.query = FakeQuery([FakeUser(id=3, email="member@example.com")])
    FakeProject.query = FakeQuery([FakeProject(id=10, name="Alpha", path="/projects/10")])
    FakeUserProject.query = FakeQuery([])
    fake_db = FakeDb()

    dbmodels = ModuleType("mdvtools.dbutils.dbmodels")
    dbmodels.User = FakeUser
    dbmodels.Project = FakeProject
    dbmodels.UserProject = FakeUserProject
    dbmodels.db = fake_db
    monkeypatch.setitem(sys.modules, "mdvtools.dbutils.dbmodels", dbmodels)

    refresh_calls: list[bool] = []
    service = MDVAdminServices(
        enable_auth=False,
        refresh_auth_cache=refresh_calls.append,
    )

    service.add_project_member(10, ProjectMemberInput(user_id=3, permission="edit"))

    assert fake_db.session.commits == 1
    assert refresh_calls == []


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
        def create_user_with_project_access(self, data: CreateAdminUserInput):
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
        def create_user_with_project_access(self, data: CreateAdminUserInput):
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


def test_create_user_maps_external_identity_errors():
    class FailingIdentityAdminServices(FakeAdminServices):
        def create_user_with_project_access(self, data: CreateAdminUserInput):
            raise AdminExternalServiceError("Auth0 operation failed.")

    flask_app = Flask(__name__)
    flask_app.secret_key = "test-secret"
    flask_app.config["ENABLE_AUTH"] = False
    AdminExtension(services=FailingIdentityAdminServices()).register_global_routes(flask_app, {})

    response = flask_app.test_client().post(
        "/admin/api/users",
        json={"email": "new.user@example.com"},
    )

    assert response.status_code == 502
    assert response.get_json()["error"] == "Auth0 operation failed."


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


class FakeAuth0Response:
    def __init__(self, payload: Any):
        self.payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self.payload


def test_auth0_identity_provider_reuses_user_in_target_connection(monkeypatch):
    import requests

    def fake_post(url: str, **kwargs: Any):
        if url == "https://example.auth0.com/oauth/token":
            return FakeAuth0Response({"access_token": "management-token"})
        raise AssertionError(f"Unexpected POST {url}")

    def fake_get(url: str, **kwargs: Any):
        assert url == "https://example.auth0.com/api/v2/users-by-email"
        return FakeAuth0Response([
            {
                "user_id": "auth0|existing",
                "identities": [{"connection": "Username-Password-Authentication"}],
            }
        ])

    monkeypatch.setattr(requests, "post", fake_post)
    monkeypatch.setattr(requests, "get", fake_get)

    provider = Auth0AdminIdentityProvider(
        domain="example.auth0.com",
        client_id="client-id",
        client_secret="client-secret",
        connection="Username-Password-Authentication",
    )

    result = provider.create_or_resolve_user(
        AdminIdentityInput(
            email=" Existing.User@Example.com ",
            first_name="Existing",
            last_name="User",
        )
    )

    assert result.email == "existing.user@example.com"
    assert result.auth_id == "auth0|existing"
    assert result.created is False


def test_auth0_identity_provider_creates_user_in_target_connection(monkeypatch):
    import requests

    post_calls: list[tuple[str, dict[str, Any]]] = []

    def fake_post(url: str, **kwargs: Any):
        post_calls.append((url, kwargs))
        if url == "https://example.auth0.com/oauth/token":
            return FakeAuth0Response({"access_token": "management-token"})
        if url == "https://example.auth0.com/api/v2/users":
            return FakeAuth0Response({"user_id": "auth0|created"})
        raise AssertionError(f"Unexpected POST {url}")

    def fake_get(url: str, **kwargs: Any):
        assert url == "https://example.auth0.com/api/v2/users-by-email"
        return FakeAuth0Response([
            {
                "user_id": "google-oauth2|existing",
                "identities": [{"connection": "google-oauth2"}],
            }
        ])

    monkeypatch.setattr(requests, "post", fake_post)
    monkeypatch.setattr(requests, "get", fake_get)

    provider = Auth0AdminIdentityProvider(
        domain="example.auth0.com",
        client_id="client-id",
        client_secret="client-secret",
        connection="Username-Password-Authentication",
    )

    result = provider.create_or_resolve_user(
        AdminIdentityInput(
            email="new.user@example.com",
            first_name="New",
            last_name="User",
        )
    )

    create_payload = post_calls[1][1]["json"]
    assert result.auth_id == "auth0|created"
    assert result.created is True
    assert create_payload["connection"] == "Username-Password-Authentication"
    assert create_payload["email"] == "new.user@example.com"
    assert create_payload["email_verified"] is True
    assert isinstance(create_payload["password"], str)
    assert len(create_payload["password"]) == 16
