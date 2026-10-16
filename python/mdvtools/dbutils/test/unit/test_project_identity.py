"""Tests for database-assigned Project IDs (ADR 0013)."""

import io
import json
import os
import zipfile

import pytest
from flask import Flask

from mdvtools.dbutils.dbmodels import Project, db
from mdvtools.dbutils.dbservice import ProjectService
from mdvtools.dbutils.mdv_server_app import serve_projects_from_db, serve_projects_from_filesystem
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension
from mdvtools.file_processing import ValidationError, mdv_project_processing
from mdvtools.mdvproject import MDVProject
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


def write_project_directory(path, name=None):
    """Write the files a directory needs to be recognised as an MDV project."""
    path.mkdir()
    state = {"all_views": [], "permission": "edit"}
    if name is not None:
        state["name"] = name
    (path / "datasources.json").write_text(json.dumps([]))
    (path / "views.json").write_text(json.dumps({}))
    (path / "state.json").write_text(json.dumps(state))
    return path


def mdv_project_archive(name=None, permission="edit"):
    """A minimal valid MDV project archive, with the required files at the root."""
    state = {"all_views": [], "permission": permission}
    if name is not None:
        state["name"] = name
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("datasources.json", json.dumps([]))
        archive.writestr("views.json", json.dumps({}))
        archive.writestr("state.json", json.dumps(state))
    buffer.seek(0)
    return buffer


def test_set_display_name_writes_the_name_into_state_json(tmp_path):
    """The display name is kept in the project directory so a rescan can rebuild
    a catalog row for a directory that has lost one."""
    project = MDVProject(str(tmp_path / "project"))

    project.set_display_name("Tumour atlas")

    assert project.state["name"] == "Tumour atlas"
    assert MDVProject(str(tmp_path / "project")).state["name"] == "Tumour atlas"


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


def test_uploaded_project_registers_the_route_under_the_assigned_id(app, tmp_path):
    """An uploaded archive is served on the ID the database gave it, and the row
    survives the app context the upload was processed in."""
    archive_path = tmp_path / "upload.zip"
    archive_path.write_bytes(mdv_project_archive().getvalue())

    with app.app_context():
        retired_id = retire_one_project_id(tmp_path)

        result = mdv_project_processing(app, str(tmp_path), str(archive_path), "upload.zip")
        assigned_id = result["project_id"]
        assert assigned_id != retired_id

    with app.app_context():
        assert str(assigned_id) in ProjectBlueprint.blueprints
        project = db.session.get(Project, assigned_id)
        assert project is not None
        assert os.path.basename(project.path) != str(assigned_id)


def test_rescan_serves_discovered_directories_under_their_assigned_ids(app, tmp_path):
    """A rescan gives each discovered directory the ID the database assigns and
    keeps the directory where it found it."""
    with app.app_context():
        retire_one_project_id(tmp_path)
        first = write_project_directory(tmp_path / "first-on-disk")
        second = write_project_directory(tmp_path / "second-on-disk")

        created_ids = serve_projects_from_filesystem(app, str(tmp_path))

    assert len(created_ids) == 2
    assert len(set(created_ids)) == 2

    with app.app_context():
        paths = set()
        for assigned_id in created_ids:
            assert str(assigned_id) in ProjectBlueprint.blueprints
            project = db.session.get(Project, assigned_id)
            assert project is not None
            paths.add(project.path)
        assert paths == {str(first), str(second)}


def test_rescan_creates_nothing_when_there_is_nothing_new(app, tmp_path):
    """An empty project root is not an error and leaves the catalog untouched."""
    with app.app_context():
        assert serve_projects_from_filesystem(app, str(tmp_path)) == []
        assert Project.query.count() == 0


def test_rescan_skips_a_failing_project_without_leaving_a_row(app, tmp_path, monkeypatch):
    """One directory that cannot be served does not stop the scan, and the row
    started for it is discarded rather than committed with the next project."""
    good = write_project_directory(tmp_path / "good-project")
    bad = write_project_directory(tmp_path / "bad-project")

    original_serve = MDVProject.serve

    def serve(self, *args, **kwargs):
        if self.dir == str(bad):
            raise RuntimeError("cannot serve this project")
        return original_serve(self, *args, **kwargs)

    monkeypatch.setattr(MDVProject, "serve", serve)

    with app.app_context():
        created_ids = serve_projects_from_filesystem(app, str(tmp_path))

    with app.app_context():
        assert len(created_ids) == 1
        assert db.session.get(Project, created_ids[0]).path == str(good)
        assert Project.query.count() == 1


def test_rescan_unregisters_a_project_whose_row_could_not_be_committed(app, tmp_path, monkeypatch):
    """Serving registers the route before the row is committed, so a failed commit
    has to take the route back out. Otherwise the URL answers for a project the
    catalog does not hold."""
    write_project_directory(tmp_path / "uncommittable")

    def commit():
        raise RuntimeError("database went away")

    with app.app_context():
        monkeypatch.setattr(db.session, "commit", commit)
        created_ids = serve_projects_from_filesystem(app, str(tmp_path))
        monkeypatch.undo()

    assert created_ids == []
    assert ProjectBlueprint.blueprints == {}
    with app.app_context():
        assert Project.query.count() == 0


def test_import_records_the_archive_permission_in_the_catalog(app, tmp_path):
    """The database decides a project's access level, so a read-only project has
    to arrive read-only in the row. A row left editable would have the next
    startup unlock the project on disk."""
    with app.app_context():
        ProjectManagerExtension().register_global_routes(app, app.config)

    client = app.test_client()
    response = client.post(
        "/import_project",
        data={"file": (mdv_project_archive(permission="view"), "project.zip"), "name": "read only"},
        content_type="multipart/form-data",
    )
    assert response.status_code == 200, response.get_data(as_text=True)

    with app.app_context():
        assert db.session.get(Project, response.json["id"]).access_level == "read-only"


def test_startup_treats_an_unrecognised_access_level_as_editable(app, tmp_path):
    """Only "read-only" locks a project. A row holding anything else keeps the
    project editable rather than silently locking it on disk."""
    with app.app_context():
        directory = write_project_directory(tmp_path / "legacy-project")
        legacy = Project()
        legacy.name = "legacy"
        legacy.path = str(directory)
        legacy.access_level = "private"
        db.session.add(legacy)
        db.session.commit()

        serve_projects_from_db(app)

    with open(directory / "state.json") as state_file:
        assert json.load(state_file).get("permission") == "edit"


def test_rescan_recovers_the_display_name_from_state_json(app, tmp_path):
    """A directory copied in from elsewhere keeps the name it had there. Without
    a name on disk there is nothing to recover, so the directory name is used."""
    with app.app_context():
        copied_in = write_project_directory(tmp_path / "e4f5a6b7c8d9", name="Tumour atlas")
        no_name = write_project_directory(tmp_path / "pilot-cohort")

        created_ids = serve_projects_from_filesystem(app, str(tmp_path))

    with app.app_context():
        names = {}
        for assigned_id in created_ids:
            project = db.session.get(Project, assigned_id)
            names[project.path] = project.name
        assert names[str(copied_in)] == "Tumour atlas"
        assert names[str(no_name)] == "pilot-cohort"


def test_rescan_writes_back_a_name_it_had_to_guess(app, tmp_path):
    """A directory with no name on disk is named after itself once, and that name
    is recorded, so copying it on again does not depend on the folder name."""
    with app.app_context():
        discovered = write_project_directory(tmp_path / "pilot-cohort")

        created_ids = serve_projects_from_filesystem(app, str(tmp_path))

    with open(discovered / "state.json") as state_file:
        state = json.load(state_file)

    with app.app_context():
        assert state.get("name") == db.session.get(Project, created_ids[0]).name


def test_rename_writes_the_new_name_to_disk(app, tmp_path):
    """A rename updates the row and the recovery copy together, so a project
    copied out afterwards carries the name it was renamed to."""
    with app.app_context():
        ProjectManagerExtension().register_global_routes(app, app.config)

    client = app.test_client()
    created = client.post("/create_project").json["id"]

    response = client.put(f"/projects/{created}/rename", data={"name": "Tumour atlas"})
    assert response.status_code == 200, response.get_data(as_text=True)

    with app.app_context():
        project = db.session.get(Project, created)
        assert project.name == "Tumour atlas"
        with open(os.path.join(project.path, "state.json")) as state_file:
            assert json.load(state_file).get("name") == "Tumour atlas"


def test_access_change_writes_the_permission_to_disk(app, tmp_path):
    """Making a project read-only updates the row and the recovery copy together,
    so the permission survives the row being lost."""
    with app.app_context():
        ProjectManagerExtension().register_global_routes(app, app.config)

    client = app.test_client()
    created = client.post("/create_project").json["id"]

    response = client.put(f"/projects/{created}/access", data={"type": "read-only"})
    assert response.status_code == 200, response.get_data(as_text=True)

    with app.app_context():
        project = db.session.get(Project, created)
        assert project.access_level == "read-only"
        with open(os.path.join(project.path, "state.json")) as state_file:
            assert json.load(state_file).get("permission") == "view"


def test_startup_writes_the_display_name_of_every_served_project(app, tmp_path):
    """A project whose directory predates the recovery copy gets its name written
    on the next boot, so copying that directory elsewhere carries the name."""
    with app.app_context():
        directory = write_project_directory(tmp_path / "1")
        existing = Project()
        existing.name = "Tumour atlas"
        existing.path = str(directory)
        db.session.add(existing)
        db.session.commit()

        serve_projects_from_db(app)

    with open(directory / "state.json") as state_file:
        assert json.load(state_file).get("name") == "Tumour atlas"


def test_startup_does_not_let_state_json_overwrite_an_access_level(app, tmp_path):
    """The database decides the access level for a project it holds a row for, so
    a stale permission on disk cannot unlock a project that was made read-only."""
    with app.app_context():
        directory = write_project_directory(tmp_path / "locked-project")
        locked = Project()
        locked.name = "locked"
        locked.path = str(directory)
        locked.access_level = "read-only"
        db.session.add(locked)
        db.session.commit()
        project_id = locked.id

        serve_projects_from_db(app)

    with app.app_context():
        assert db.session.get(Project, project_id).access_level == "read-only"
        with open(directory / "state.json") as state_file:
            assert json.load(state_file).get("permission") == "view"


def unsafe_archive():
    """An archive whose entry escapes the directory it is extracted into."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("../escaped.json", "{}")
    buffer.seek(0)
    return buffer


def test_a_rejected_upload_leaves_no_directory_behind(app, tmp_path):
    """A rejected archive must not leave its directory in the Project root, where
    no catalog row points at it and no later scan can tell it from a project."""
    projects_root = tmp_path / "projects"
    projects_root.mkdir()
    archive_path = tmp_path / "unsafe.zip"
    archive_path.write_bytes(unsafe_archive().getvalue())

    with app.app_context():
        with pytest.raises(ValidationError):
            mdv_project_processing(app, str(projects_root), str(archive_path), "unsafe.zip")

    assert list(projects_root.iterdir()) == []


def test_a_rejected_import_leaves_no_directory_behind(app, tmp_path):
    """The import route rejects the same archives and has the same directory to
    clean up."""
    with app.app_context():
        ProjectManagerExtension().register_global_routes(app, app.config)

    client = app.test_client()
    response = client.post(
        "/import_project",
        data={"file": (unsafe_archive(), "unsafe.zip")},
        content_type="multipart/form-data",
    )

    assert response.status_code == 400
    assert list(tmp_path.iterdir()) == []


def test_creation_paths_record_the_display_name_on_disk(app, tmp_path):
    """Every path that creates a project directory writes the display name into
    it, so the name survives the row being lost or the directory being copied."""
    archive_path = tmp_path / "upload.zip"
    archive_path.write_bytes(mdv_project_archive().getvalue())

    with app.app_context():
        ProjectManagerExtension().register_global_routes(app, app.config)

    client = app.test_client()
    created = client.post("/create_project").json["id"]
    imported = client.post(
        "/import_project",
        data={"file": (mdv_project_archive(), "project.zip"), "name": "imported project"},
        content_type="multipart/form-data",
    ).json["id"]

    with app.app_context():
        uploaded = mdv_project_processing(app, str(tmp_path), str(archive_path), "upload.zip")["project_id"]

    with app.app_context():
        for assigned_id in (created, imported, uploaded):
            project = db.session.get(Project, assigned_id)
            with open(os.path.join(project.path, "state.json")) as state_file:
                state = json.load(state_file)
            assert state.get("name") == project.name


def test_import_and_upload_take_the_name_from_the_archive(app, tmp_path):
    """An archive carrying a display name keeps it when the caller supplies no
    name, instead of the name being replaced by the default."""
    archive_path = tmp_path / "upload.zip"
    archive_path.write_bytes(mdv_project_archive(name="Pilot cohort").getvalue())

    with app.app_context():
        ProjectManagerExtension().register_global_routes(app, app.config)

    response = app.test_client().post(
        "/import_project",
        data={"file": (mdv_project_archive(name="Tumour atlas"), "project.zip")},
        content_type="multipart/form-data",
    )
    assert response.status_code == 200, response.get_data(as_text=True)
    imported = response.json["id"]

    with app.app_context():
        uploaded = mdv_project_processing(app, str(tmp_path), str(archive_path), "upload.zip")["project_id"]

    with app.app_context():
        assert db.session.get(Project, imported).name == "Tumour atlas"
        assert db.session.get(Project, uploaded).name == "Pilot cohort"
