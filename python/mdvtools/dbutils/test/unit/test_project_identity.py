"""Tests for database-assigned Project IDs (ADR 0013)."""

import pytest
from flask import Flask

from mdvtools.dbutils.dbmodels import Project, db
from mdvtools.dbutils.dbservice import ProjectService


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
