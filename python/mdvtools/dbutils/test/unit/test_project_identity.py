"""Tests for database-assigned Project IDs (ADR 0013)."""

import io
import json
import os
import zipfile

import pytest
from flask import Flask

from mdvtools.dbutils.dbmodels import Project, db
from mdvtools.dbutils.dbservice import ProjectService
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension
from mdvtools.project_router import ProjectBlueprint
from mdvtools.websocket import mdv_socketio


@pytest.fixture()
def app(tmp_path):
    app = Flask(__name__)
    # MDVProject.serve() initialises the SocketIO upload handlers, which need the
    # SocketIO instance the real app creates at startup.
    mdv_socketio(app)
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


@pytest.fixture(autouse=True)
def clear_route_registry():
    """ProjectBlueprint.blueprints is class-level, so it outlives a test."""
    ProjectBlueprint.blueprints.clear()
    yield
    ProjectBlueprint.blueprints.clear()


def retire_one_project_id(tmp_path):
    """Create and purge a project so the next assigned ID is no longer max(id)+1."""
    path = tmp_path / "retired"
    path.mkdir()
    project = ProjectService.add_new_project(path=str(path), name="retired")
    db.session.commit()
    retired_id = project.id
    ProjectService.soft_delete_projects([retired_id])
    purged, message = ProjectService.purge_deleted_project(retired_id)
    assert purged, message
    return retired_id


def mdv_project_archive():
    """A minimal valid MDV project archive, with the required files at the root."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("datasources.json", json.dumps([]))
        archive.writestr("views.json", json.dumps({}))
        archive.writestr("state.json", json.dumps({"all_views": [], "permission": "edit"}))
    buffer.seek(0)
    return buffer


def test_purged_project_id_is_not_reused(app, tmp_path):
    """A Project ID retired by a purge is never given to a later project.

    SQLite reuses max(rowid) once the highest row is deleted unless the table is
    declared AUTOINCREMENT, which would give the second project below the first
    project's ID and therefore its URL.
    """
    with app.app_context():
        first_path = tmp_path / "first"
        first_path.mkdir()
        first = ProjectService.add_new_project(path=str(first_path), name="first")
        db.session.commit()
        first_id = first.id

        ProjectService.soft_delete_projects([first_id])
        purged, message = ProjectService.purge_deleted_project(first_id)
        assert purged, message

        second_path = tmp_path / "second"
        second_path.mkdir()
        second = ProjectService.add_new_project(path=str(second_path), name="second")
        db.session.commit()

        assert second.id != first_id
        assert db.session.get(Project, first_id) is None


def test_add_new_project_assigns_an_id_without_committing(app, tmp_path):
    """The Project ID is readable as soon as the row is added, but the row is
    only durable once the caller commits, so a failed creation can roll back."""
    with app.app_context():
        path = tmp_path / "rolled-back"
        path.mkdir()

        project = ProjectService.add_new_project(path=str(path), name="rolled back")
        assigned_id = project.id
        assert assigned_id is not None

        db.session.rollback()

        assert db.session.get(Project, assigned_id) is None


def test_create_project_registers_the_route_under_the_assigned_id(app, tmp_path):
    """The route a new project is served on is the ID the database gave it, and
    its directory is not named after that ID."""
    with app.app_context():
        retired_id = retire_one_project_id(tmp_path)
        ProjectManagerExtension().register_global_routes(app, app.config)

    response = app.test_client().post("/create_project")
    assert response.status_code == 200, response.get_data(as_text=True)
    assigned_id = response.json["id"]
    assert assigned_id != retired_id

    with app.app_context():
        project = db.session.get(Project, assigned_id)
        assert project is not None
        assert str(assigned_id) in ProjectBlueprint.blueprints
        assert os.path.basename(project.path) != str(assigned_id)


def test_import_project_registers_the_route_under_the_assigned_id(app, tmp_path):
    """An imported project is served on the ID the database gave it, and its
    files land in a directory that is not named after that ID."""
    with app.app_context():
        retired_id = retire_one_project_id(tmp_path)
        ProjectManagerExtension().register_global_routes(app, app.config)

    response = app.test_client().post(
        "/import_project",
        data={"file": (mdv_project_archive(), "project.zip"), "name": "imported"},
        content_type="multipart/form-data",
    )
    assert response.status_code == 200, response.get_data(as_text=True)
    assigned_id = response.json["id"]
    assert assigned_id != retired_id

    with app.app_context():
        assert str(assigned_id) in ProjectBlueprint.blueprints
        project = db.session.get(Project, assigned_id)
        assert project is not None
        assert os.path.basename(project.path) != str(assigned_id)
        assert os.path.exists(os.path.join(project.path, "datasources.json"))
