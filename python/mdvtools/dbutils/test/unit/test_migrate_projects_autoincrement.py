"""Tests for the script that rebuilds the projects table with AUTOINCREMENT."""

import sqlite3

import pytest
from sqlalchemy.dialects import sqlite
from sqlalchemy.schema import CreateTable

from mdvtools.dbutils.dbmodels import Project
from mdvtools.scripts.migrate_projects_autoincrement import migrate_projects_table


# Every NOT NULL column on Project has a Python-side default only, so raw SQL has
# to supply all of them.
INSERT_PROJECT = (
    "INSERT INTO projects (name, path, created_timestamp, is_deleted, update_timestamp,"
    " accessed_timestamp, is_public, access_level)"
    " VALUES (?, ?, '2026-01-01', 0, '2026-01-01', '2026-01-01', 0, 'editable')"
)


def create_pre_migration_database(db_path):
    """Build a projects table as it exists before this migration.

    The model now carries sqlite_autoincrement, so the earlier schema is the
    same DDL with that keyword removed.
    """
    ddl = str(CreateTable(Project.__table__).compile(dialect=sqlite.dialect()))
    ddl = ddl.replace(" AUTOINCREMENT", "")
    assert "AUTOINCREMENT" not in ddl

    connection = sqlite3.connect(db_path)
    try:
        connection.execute(ddl)
        for name in ("Tumour atlas", "Pilot cohort"):
            connection.execute(INSERT_PROJECT, (name, f"/app/mdv/{name}"))
        connection.commit()
    finally:
        connection.close()


def highest_row_then_insert(db_path, name):
    """Delete the highest project, insert another, and return the new ID."""
    connection = sqlite3.connect(db_path)
    try:
        connection.execute("DELETE FROM projects WHERE id = (SELECT max(id) FROM projects)")
        cursor = connection.execute(INSERT_PROJECT, (name, f"/app/mdv/{name}"))
        connection.commit()
        return cursor.lastrowid
    finally:
        connection.close()


@pytest.fixture()
def db_path(tmp_path):
    path = tmp_path / "mdv.sqlite"
    create_pre_migration_database(str(path))
    return str(path)


def test_projects_table_reuses_ids_before_migration(db_path):
    """The state this migration exists to fix."""
    assert highest_row_then_insert(db_path, "Third project") == 2


def test_migration_stops_ids_being_reused(db_path):
    migrate_projects_table(db_path)

    assert highest_row_then_insert(db_path, "Third project") == 3


def test_migration_keeps_every_project(db_path):
    migrate_projects_table(db_path)

    connection = sqlite3.connect(db_path)
    try:
        rows = connection.execute("SELECT id, name, path, access_level FROM projects ORDER BY id").fetchall()
    finally:
        connection.close()

    assert rows == [
        (1, "Tumour atlas", "/app/mdv/Tumour atlas", "editable"),
        (2, "Pilot cohort", "/app/mdv/Pilot cohort", "editable"),
    ]


def test_migration_is_safe_to_run_twice(db_path):
    first = migrate_projects_table(db_path)
    second = migrate_projects_table(db_path)

    assert first["migrated"] is True
    assert second["migrated"] is False
    assert highest_row_then_insert(db_path, "Third project") == 3
