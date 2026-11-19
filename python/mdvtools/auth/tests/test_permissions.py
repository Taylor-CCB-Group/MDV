"""Project roles are read from the database each time they are checked.

A role written by another process, such as MDV-Admin on a shared database or a
script run with `docker exec`, applies on the next check without a restart.
The file-backed database lets a plain sqlite3 connection play that process.
"""

import sqlite3

import pytest
from flask import Flask

from mdvtools.auth.permissions import can_view, project_owner_emails, project_role, project_roles
from mdvtools.dbutils.dbmodels import Project, User, UserProject, db


@pytest.fixture()
def db_path(tmp_path):
    return tmp_path / "mdv.sqlite3"


@pytest.fixture()
def app(db_path):
    app = Flask(__name__)
    app.config.update(
        SQLALCHEMY_DATABASE_URI=f"sqlite:///{db_path}",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    db.init_app(app)
    with app.app_context():
        db.create_all()
        db.session.add_all([
            Project(id=1, name="atlas", path="/projects/atlas"),
            Project(id=2, name="census", path="/projects/census"),
            User(id=7, email="reader@example.org", auth_id="reader"),
            User(id=8, email="owner-b@example.org", auth_id="owner-b"),
            User(id=9, email="owner-a@example.org", auth_id="owner-a"),
        ])
        db.session.add_all([
            UserProject(7, 1, can_read=True),
            UserProject(8, 1, is_owner=True),
            UserProject(9, 1, is_owner=True),
        ])
        db.session.commit()
        yield app
        db.session.remove()
        db.drop_all()


def test_project_role_reads_the_users_row(app):
    assert project_role(7, 1) == {"can_read": True, "can_write": False, "is_owner": False}
    assert project_role(7, 2) is None


def test_project_roles_lists_every_project_the_user_holds_a_role_on(app):
    assert project_roles(7) == {1: {"can_read": True, "can_write": False, "is_owner": False}}
    assert project_roles(99) == {}


def test_can_view_needs_read_write_or_owner(app):
    assert not can_view(None)
    assert not can_view({"can_read": False, "can_write": False, "is_owner": False})
    assert can_view({"can_read": True, "can_write": False, "is_owner": False})
    assert can_view({"can_read": False, "can_write": True, "is_owner": False})
    assert can_view({"can_read": False, "can_write": False, "is_owner": True})


def test_project_owner_emails_are_sorted_and_skip_unowned_projects(app):
    assert project_owner_emails([1, 2]) == {1: ["owner-a@example.org", "owner-b@example.org"]}
    assert project_owner_emails([]) == {}


def test_a_role_written_by_another_process_applies_on_the_next_check(app, db_path):
    assert project_role(7, 2) is None

    with sqlite3.connect(db_path) as other_process:
        other_process.execute(
            "INSERT INTO user_projects (user_id, project_id, can_read, can_write, is_owner)"
            " VALUES (7, 2, 1, 0, 0)"
        )
        other_process.execute("DELETE FROM user_projects WHERE user_id = 7 AND project_id = 1")
    # Flask-SQLAlchemy ends the session after each request, so the next check starts fresh.
    db.session.remove()

    assert can_view(project_role(7, 2))
    assert project_role(7, 1) is None
