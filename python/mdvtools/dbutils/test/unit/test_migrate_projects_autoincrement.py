"""Tests for the script that rebuilds a SQLite database so Project IDs are never reused."""

import sqlite3
from contextlib import closing
from typing import NamedTuple

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


class Extra(NamedTuple):
    """Something a deployment added that the models do not declare."""

    statements: list[str]
    name: str
    # Run against the migrated database to show the extra came across.
    query: str
    expected: list[tuple]


EXTRAS = {
    "table": Extra(
        [
            "CREATE TABLE legacy_notes (id INTEGER PRIMARY KEY, note TEXT)",
            "INSERT INTO legacy_notes (note) VALUES ('kept')",
        ],
        "legacy_notes",
        "SELECT note FROM legacy_notes",
        [("kept",)],
    ),
    "column": Extra(
        [
            "ALTER TABLE projects ADD COLUMN legacy_note TEXT",
            "UPDATE projects SET legacy_note = 'kept' WHERE id = 1",
        ],
        "legacy_note",
        "SELECT legacy_note FROM projects WHERE id = 1",
        [("kept",)],
    ),
    "index": Extra(
        ["CREATE INDEX idx_projects_name ON projects (name)"],
        "idx_projects_name",
        "SELECT type FROM sqlite_master WHERE name = 'idx_projects_name'",
        [("index",)],
    ),
    "trigger": Extra(
        [
            "CREATE TRIGGER touch_project AFTER UPDATE OF name ON projects"
            " BEGIN UPDATE projects SET update_timestamp = '2026-09-14' WHERE id = new.id; END"
        ],
        "touch_project",
        "SELECT type FROM sqlite_master WHERE name = 'touch_project'",
        [("trigger",)],
    ),
    "view": Extra(
        ["CREATE VIEW live_projects AS SELECT name FROM projects WHERE is_deleted = 0"],
        "live_projects",
        "SELECT name FROM live_projects ORDER BY name",
        [("Pilot cohort",), ("Tumour atlas",)],
    ),
}


def compile_sqlite(element):
    return str(element.compile(dialect=sqlite.dialect()))


def table_rows(path, table):
    with closing(sqlite3.connect(path)) as connection:
        return connection.execute(f'SELECT * FROM "{table}" ORDER BY rowid').fetchall()


def directory_contents(directory):
    return {path.name: path.read_bytes() for path in directory.iterdir()}


def run_sql(path, statements):
    with closing(sqlite3.connect(path)) as connection:
        for statement in statements:
            connection.execute(statement)
        connection.commit()


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


def test_an_already_migrated_database_reports_nothing_to_do_and_writes_nothing(db_path, tmp_path, capsys):
    assert main([str(db_path)]) == 0
    capsys.readouterr()
    before = directory_contents(tmp_path)

    assert main([str(db_path)]) == 0

    assert "nothing to do" in capsys.readouterr().out
    assert directory_contents(tmp_path) == before


def test_a_dry_run_writes_nothing_and_reports_the_next_id(db_path, tmp_path, capsys):
    before = directory_contents(tmp_path)

    assert main([str(db_path), "--dry-run"]) == 0

    assert "Next project ID: 3" in capsys.readouterr().out
    assert directory_contents(tmp_path) == before


def test_min_next_id_sets_the_next_id(db_path):
    assert main([str(db_path), "--min-next-id", "100"]) == 0

    with closing(sqlite3.connect(db_path)) as connection:
        new_id = connection.execute(INSERT_PROJECT, ("Third project", "/app/mdv/third")).lastrowid
        connection.commit()

    assert new_id == 100


@pytest.mark.parametrize("extra", EXTRAS.values(), ids=list(EXTRAS))
def test_something_the_models_do_not_declare_stops_the_run(db_path, tmp_path, capsys, extra):
    run_sql(db_path, extra.statements)
    before = directory_contents(tmp_path)

    assert main([str(db_path)]) == 1

    assert extra.name in capsys.readouterr().err
    assert directory_contents(tmp_path) == before


@pytest.mark.parametrize("extra", EXTRAS.values(), ids=list(EXTRAS))
def test_keep_extras_brings_across_something_the_models_do_not_declare(db_path, extra):
    run_sql(db_path, extra.statements)

    assert main([str(db_path), "--keep-extras"]) == 0

    with closing(sqlite3.connect(db_path)) as connection:
        assert connection.execute(extra.query).fetchall() == extra.expected
