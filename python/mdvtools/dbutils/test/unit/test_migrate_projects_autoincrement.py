"""Tests for the script that rebuilds a SQLite database so Project IDs are never reused."""

import sqlite3
from contextlib import closing

import pytest
from sqlalchemy.dialects import sqlite
from sqlalchemy.schema import CreateIndex, CreateTable

from mdvtools.dbutils.dbmodels import db
from mdvtools.scripts.migrate_projects_autoincrement import main

# Every NOT NULL column in the models has a Python-side default only, so raw SQL
# has to supply all of them.
INSERT_USER = (
    "INSERT INTO users (email, password, is_active, first_name, last_name, administrator,"
    " auth_id, is_admin) VALUES ('admin@example.com', '', 1, '', '', 1, 'admin', 1)"
)
INSERT_PROJECT = (
    "INSERT INTO projects (name, path, created_timestamp, is_deleted, update_timestamp,"
    " accessed_timestamp, is_public, access_level)"
    " VALUES (?, ?, '2026-01-01', 0, '2026-01-01', '2026-01-01', 0, 'editable')"
)
INSERT_OWNER = (
    "INSERT INTO user_projects (user_id, project_id, can_read, can_write, is_owner)"
    " VALUES (1, ?, 1, 1, 1)"
)


def compile_sqlite(element):
    return str(element.compile(dialect=sqlite.dialect()))


def table_rows(path, table):
    with closing(sqlite3.connect(path)) as connection:
        return connection.execute(f'SELECT * FROM "{table}" ORDER BY rowid').fetchall()


@pytest.fixture()
def db_path(tmp_path):
    """A database as the app created it before projects was declared AUTOINCREMENT.

    The models now carry that declaration, so the earlier schema is the models'
    DDL with the keyword removed.
    """
    path = tmp_path / "mdv.sqlite3"
    with closing(sqlite3.connect(path)) as connection:
        for table in db.metadata.sorted_tables:
            connection.execute(compile_sqlite(CreateTable(table)).replace(" AUTOINCREMENT", ""))
            for index in table.indexes:
                connection.execute(compile_sqlite(CreateIndex(index)))
        connection.execute(INSERT_USER)
        for name in ("Tumour atlas", "Pilot cohort"):
            project_id = connection.execute(INSERT_PROJECT, (name, f"/app/mdv/{name}")).lastrowid
            connection.execute(INSERT_OWNER, (project_id,))
        connection.commit()
    return path


def test_a_purged_project_id_is_not_reused_after_the_migration(db_path):
    assert main([str(db_path)]) == 0

    with closing(sqlite3.connect(db_path)) as connection:
        connection.execute("DELETE FROM user_projects WHERE project_id = 2")
        connection.execute("DELETE FROM projects WHERE id = 2")
        new_id = connection.execute(INSERT_PROJECT, ("Third project", "/app/mdv/third")).lastrowid
        connection.commit()

    assert new_id == 3


def test_every_table_keeps_its_rows_and_the_original_is_kept_unchanged(db_path, tmp_path):
    original = db_path.read_bytes()

    assert main([str(db_path)]) == 0

    backup = tmp_path / "mdv.sqlite3.pre-autoincrement"
    assert sorted(path.name for path in tmp_path.iterdir()) == [db_path.name, backup.name]
    assert backup.read_bytes() == original
    for table in db.metadata.sorted_tables:
        assert table_rows(db_path, table.name) == table_rows(backup, table.name), table.name
