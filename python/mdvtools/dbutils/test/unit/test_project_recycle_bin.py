import os

import pytest
from flask import Flask

from mdvtools.dbutils.dbmodels import File, Project, User, UserProject, db
from mdvtools.dbutils.dbservice import ProjectService
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension
from mdvtools.project_router import ProjectBlueprint


@pytest.fixture()
def app(tmp_path):
    app = Flask(__name__)
    app.config.update(
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
        ENABLE_AUTH=False,
        projects_base_dir=str(tmp_path),
    )
    db.init_app(app)
    with app.app_context():
        db.create_all()
        yield app
        db.session.remove()
        db.drop_all()


def add_project(tmp_path, project_id, *, deleted=False, access_level="editable"):
    path = tmp_path / str(project_id)
    path.mkdir()
    project = Project(
        id=project_id,
        name=f"project-{project_id}",
        path=str(path),
        is_deleted=deleted,
        access_level=access_level,
    )
    db.session.add(project)
    db.session.commit()
    return project


def test_soft_delete_projects_is_atomic(app, tmp_path):
    with app.app_context():
        project = add_project(tmp_path, 1)
        add_project(tmp_path, 2, access_level="view")

        with pytest.raises(ValueError):
            ProjectService.soft_delete_projects([project.id, 2])

        assert not db.session.get(Project, project.id).is_deleted


def test_get_deleted_projects_filters_by_owner(app, tmp_path):
    with app.app_context():
        owned = add_project(tmp_path, 1, deleted=True)
        add_project(tmp_path, 2, deleted=True)
        user = User(email="owner@example.com", auth_id="owner")
        db.session.add(user)
        db.session.flush()
        db.session.add(UserProject(user.id, owned.id, is_owner=True))
        db.session.commit()

        assert [project.id for project in ProjectService.get_deleted_projects(user.id)] == [owned.id]


def test_purge_deleted_project_removes_dependencies_and_allows_missing_path(app, tmp_path):
    with app.app_context():
        project = add_project(tmp_path, 1, deleted=True)
        user = User(email="owner@example.com", auth_id="owner")
        db.session.add(user)
        db.session.flush()
        db.session.add(UserProject(user.id, project.id, is_owner=True))
        db.session.add(File(name="state", file_path="/tmp/state", project_id=project.id))
        db.session.commit()
        os.rmdir(project.path)

        success, message = ProjectService.purge_deleted_project(project.id)

        assert success
        assert message is None
        assert db.session.get(Project, project.id) is None
        assert File.query.count() == 0
        assert UserProject.query.count() == 0


def test_purge_deleted_project_keeps_record_when_filesystem_delete_fails(app, tmp_path, monkeypatch):
    with app.app_context():
        project = add_project(tmp_path, 1, deleted=True)
        monkeypatch.setattr("mdvtools.dbutils.dbservice.shutil.rmtree", lambda path: (_ for _ in ()).throw(OSError("busy")))

        success, message = ProjectService.purge_deleted_project(project.id)

        assert not success
        assert message == "busy"
        assert db.session.get(Project, project.id) is not None


def test_purge_deleted_project_rejects_dangerous_missing_path(app, tmp_path):
    with app.app_context():
        project = add_project(tmp_path, 1, deleted=True)
        project.path = ""
        db.session.commit()

        success, message = ProjectService.purge_deleted_project(project.id)

        assert not success
        assert message == "Project path is missing."
        assert db.session.get(Project, project.id) is not None


def test_recycle_bin_routes_soft_delete_and_purge(app, tmp_path):
    with app.app_context():
        add_project(tmp_path, 1)
        ProjectBlueprint.blueprints["1"] = object()
        ProjectManagerExtension().register_global_routes(app, app.config)

    client = app.test_client()
    response = client.post("/projects/soft-delete", json={"projectIds": [1, 1]})
    assert response.status_code == 200
    assert response.json == {"deletedProjectIds": [1]}
    assert "1" not in ProjectBlueprint.blueprints

    response = client.get("/projects/recycle-bin")
    assert response.status_code == 200
    assert [project["id"] for project in response.json] == [1]

    response = client.delete("/projects/recycle-bin")
    assert response.status_code == 200
    assert response.json == {"deletedProjectIds": [1], "failures": []}


def test_recycle_bin_route_restores_project(app, tmp_path, monkeypatch):
    class DummyMDVProject:
        def __init__(self, *args, **kwargs):
            pass

        def serve(self, *args, **kwargs):
            return None

    monkeypatch.setattr("mdvtools.dbutils.project_manager_service.MDVProject", DummyMDVProject)
    with app.app_context():
        add_project(tmp_path, 1, deleted=True)
        ProjectManagerExtension().register_global_routes(app, app.config)

    response = app.test_client().post("/projects/recycle-bin/1/restore")
    assert response.status_code == 200
    assert response.json == {"restoredProjectId": 1}
    with app.app_context():
        project = db.session.get(Project, 1)
        assert project is not None
        assert not project.is_deleted
        assert project.deleted_timestamp is None


def test_batch_soft_delete_rejects_malformed_ids(app):
    ProjectManagerExtension().register_global_routes(app, app.config)
    response = app.test_client().post("/projects/soft-delete", json={"projectIds": []})
    assert response.status_code == 400


def test_auth_recycle_bin_routes_reject_missing_session(app, monkeypatch):
    from mdvtools.auth import authutils

    app.config.update(ENABLE_AUTH=True, SECRET_KEY="test")
    monkeypatch.setattr(authutils, "user_project_cache", {})
    monkeypatch.setattr(authutils, "active_projects_cache", [])
    monkeypatch.setattr(authutils, "user_cache", {})
    monkeypatch.setattr(authutils, "all_users_cache", [])
    monkeypatch.setattr(authutils, "cache_user_projects", lambda: None)
    ProjectManagerExtension().register_global_routes(app, app.config)

    client = app.test_client()
    assert client.get("/projects/recycle-bin").status_code == 401
    assert client.delete("/projects/recycle-bin").status_code == 401


def test_auth_routes_require_owner_filter_bin_and_refresh_cache(app, tmp_path, monkeypatch):
    from mdvtools.auth import authutils

    cache_refreshes = []
    app.config.update(ENABLE_AUTH=True, SECRET_KEY="test")
    with app.app_context():
        owned = add_project(tmp_path, 1)
        other = add_project(tmp_path, 2, deleted=True)
        owner = User(email="owner@example.com", auth_id="owner")
        other_owner = User(email="other@example.com", auth_id="other")
        db.session.add_all([owner, other_owner])
        db.session.flush()
        db.session.add(UserProject(owner.id, owned.id, is_owner=True))
        db.session.add(UserProject(other_owner.id, other.id, is_owner=True))
        db.session.commit()
        owner_id = owner.id
        owned_id = owned.id

    monkeypatch.setattr(authutils, "user_project_cache", {
        owner_id: {
            owned_id: {"can_read": True, "can_write": True, "is_owner": True},
        },
    })
    monkeypatch.setattr(authutils, "active_projects_cache", [])
    monkeypatch.setattr(authutils, "user_cache", {})
    monkeypatch.setattr(authutils, "all_users_cache", [])
    monkeypatch.setattr(authutils, "cache_user_projects", lambda: cache_refreshes.append(True))
    ProjectManagerExtension().register_global_routes(app, app.config)

    client = app.test_client()
    with client.session_transaction() as session:
        session["user"] = {"id": owner_id}

    response = client.post("/projects/soft-delete", json={"projectIds": [2]})
    assert response.status_code == 403

    response = client.post("/projects/soft-delete", json={"projectIds": [1]})
    assert response.status_code == 200
    assert cache_refreshes == [True]

    response = client.get("/projects/recycle-bin")
    assert response.status_code == 200
    assert [project["id"] for project in response.json] == [1]
