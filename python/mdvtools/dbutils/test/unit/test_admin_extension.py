from __future__ import annotations

import sys
from dataclasses import dataclass
from types import ModuleType
from typing import Any

import pytest
from flask import Flask

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
