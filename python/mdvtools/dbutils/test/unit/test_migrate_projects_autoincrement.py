"""Tests for the script that rebuilds the projects table with AUTOINCREMENT."""

import sqlite3

import pytest
from sqlalchemy.dialects import sqlite
from sqlalchemy.schema import CreateTable

from mdvtools.dbutils.dbmodels import Project
from mdvtools.scripts.migrate_projects_autoincrement import BACKUP_SUFFIX, migrate_projects_table


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


def test_backup_keeps_rows_still_in_the_write_ahead_log(tmp_path):
    """A database in WAL mode holds committed rows in the -wal file until a
    checkpoint, so a backup that copies only the main file can lose them."""
    path = str(tmp_path / "wal.sqlite")
    create_pre_migration_database(path)

    # Kept open for the whole test, because closing the last connection
    # checkpoints the write-ahead log and hides what this is testing.
    writer = sqlite3.connect(path)
    try:
        writer.execute("PRAGMA journal_mode = WAL")
        writer.execute(INSERT_PROJECT, ("Written in WAL", "/app/mdv/wal"))
        writer.commit()

        migrate_projects_table(path)

        backup = sqlite3.connect(path + BACKUP_SUFFIX)
        try:
            names = [row[0] for row in backup.execute("SELECT name FROM projects ORDER BY id")]
        finally:
            backup.close()
    finally:
        writer.close()

    assert "Written in WAL" in names


def test_migration_keeps_an_index_the_deployment_added(db_path):
    """Rebuilding the table drops every index it had, so one added by hand has to
    be put back rather than lost to the rebuild."""
    connection = sqlite3.connect(db_path)
    try:
        connection.execute("CREATE INDEX idx_projects_name ON projects (name)")
        connection.commit()
    finally:
        connection.close()

    migrate_projects_table(db_path)

    connection = sqlite3.connect(db_path)
    try:
        indexes = {
            row[0] for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'index' AND tbl_name = 'projects'"
            )
        }
    finally:
        connection.close()

    assert "idx_projects_name" in indexes
    assert {"idx_projects_genome", "idx_projects_owner"} <= indexes


def test_an_id_another_table_still_points_at_is_not_handed_out(db_path):
    """max(id) on the projects table is not the highest ID ever assigned. A row
    elsewhere pointing at a project is evidence that its ID was used."""
    connection = sqlite3.connect(db_path)
    try:
        connection.execute(
            "CREATE TABLE user_projects (id INTEGER PRIMARY KEY, user_id INTEGER,"
            " project_id INTEGER NOT NULL REFERENCES projects (id))"
        )
        connection.execute("INSERT INTO user_projects (user_id, project_id) VALUES (1, 7)")
        connection.commit()
    finally:
        connection.close()

    migrate_projects_table(db_path)

    connection = sqlite3.connect(db_path)
    try:
        cursor = connection.execute(INSERT_PROJECT, ("Third project", "/app/mdv/third"))
        connection.commit()
        assigned = cursor.lastrowid
    finally:
        connection.close()

    assert assigned == 8


def test_min_next_id_reserves_ids_that_no_row_remembers(db_path):
    """A project purged before the migration leaves nothing behind to show its ID
    was used, so an operator who knows the highest ID ever served supplies it."""
    migrate_projects_table(db_path, min_next_id=100)

    connection = sqlite3.connect(db_path)
    try:
        cursor = connection.execute(INSERT_PROJECT, ("Third project", "/app/mdv/third"))
        connection.commit()
        assigned = cursor.lastrowid
    finally:
        connection.close()

    assert assigned == 100
