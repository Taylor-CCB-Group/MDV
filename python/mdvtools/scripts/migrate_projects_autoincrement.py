#!/usr/bin/env python3
"""Rebuild a SQLite MDV database so Project IDs are never reused.

SQLite gives a deleted project's ID to the next project unless the projects table
is declared AUTOINCREMENT, and it only reads that declaration when the table is
created. A database created before the models declared it keeps reusing IDs until
it is rebuilt. This script builds a new database file from the models, copies every
row into it and swaps it in, keeping the original beside it as
<name>.pre-autoincrement.

Stop the app before migrating, and run the script in a one-off container:

    docker compose stop <service>
    docker compose run --rm --no-deps <service> uv run python mdvtools/scripts/migrate_projects_autoincrement.py
    docker compose start <service>

Do not run the migration with docker compose exec in the running container. The app
holds the database open, so anything it writes after the copy would be lost.

To roll back, stop the service, move <name>.pre-autoincrement back to the live name,
and start the service.
"""

import argparse
import os
import shutil
import sqlite3
import sys
import tempfile
from contextlib import closing
from pathlib import Path

# Add the parent directory to sys.path to import mdvtools modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sqlalchemy.dialects import sqlite as sqlite_dialect
from sqlalchemy.schema import CreateIndex, CreateTable

from mdvtools.dbutils.dbmodels import db

DEFAULT_DATABASE = "/app/mdv/mdv.sqlite3"
BACKUP_SUFFIX = ".pre-autoincrement"


def _compile(element):
    return str(element.compile(dialect=sqlite_dialect.dialect()))


def _fsync(path):
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _build_from_models(new_path, original_path):
    """Create the models' tables in new_path, copy the original's rows into them and
    check the result.

    Returns the rows in each table and the next Project ID.
    """
    with closing(sqlite3.connect(new_path.as_uri(), uri=True, isolation_level=None)) as connection:
        connection.execute("ATTACH DATABASE ? AS original", (original_path.as_uri() + "?mode=ro",))

        connection.execute("BEGIN")
        for table in db.metadata.sorted_tables:
            connection.execute(_compile(CreateTable(table)))
            for index in table.indexes:
                connection.execute(_compile(CreateIndex(index)))
            columns = ", ".join(f'"{column.name}"' for column in table.columns)
            connection.execute(
                f'INSERT INTO main."{table.name}" ({columns})'
                f' SELECT {columns} FROM original."{table.name}"'
            )
        # SQLite does not promise that inserting explicit ids moves the sequence,
        # so it is set here.
        highest_id = connection.execute("SELECT max(id) FROM main.projects").fetchone()[0] or 0
        connection.execute("DELETE FROM main.sqlite_sequence WHERE name = 'projects'")
        connection.execute(
            "INSERT INTO main.sqlite_sequence (name, seq) VALUES ('projects', ?)", (highest_id,)
        )
        connection.execute("COMMIT")

        rows = {}
        for table in db.metadata.sorted_tables:
            copied = connection.execute(f'SELECT count(*) FROM main."{table.name}"').fetchone()[0]
            expected = connection.execute(f'SELECT count(*) FROM original."{table.name}"').fetchone()[0]
            if copied != expected:
                raise RuntimeError(f"Copied {copied} of {expected} rows from {table.name}")
            rows[table.name] = copied

        # integrity_check, unlike quick_check, verifies UNIQUE constraints.
        integrity = connection.execute("PRAGMA main.integrity_check").fetchall()
        if integrity != [("ok",)]:
            raise RuntimeError(f"The new database failed its integrity check: {integrity}")

        sequence = connection.execute(
            "SELECT seq FROM main.sqlite_sequence WHERE name = 'projects'"
        ).fetchall()
        if sequence != [(highest_id,)]:
            raise RuntimeError(
                f"The projects sequence is {sequence} after the copy, expected {highest_id}"
            )

        connection.execute("DETACH DATABASE original")
    return rows, highest_id + 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "database",
        nargs="?",
        default=(os.getenv("SQLITE_DB_PATH") or DEFAULT_DATABASE).strip(),
        help=f"the SQLite database file, by default SQLITE_DB_PATH or {DEFAULT_DATABASE}",
    )
    args = parser.parse_args(argv)

    db_path = Path(args.database).resolve()
    backup_path = db_path.with_name(db_path.name + BACKUP_SUFFIX)
    if backup_path.exists():
        print(f"{backup_path} already exists. Move it away before migrating.", file=sys.stderr)
        return 1

    descriptor, new_name = tempfile.mkstemp(dir=db_path.parent, prefix=f"{db_path.name}.", suffix=".tmp")
    os.close(descriptor)
    new_path = Path(new_name)

    rows, next_id = _build_from_models(new_path, db_path)

    # The new file replaces the original under the same name, so it takes the
    # original's permissions and reaches the disk before either rename.
    shutil.copymode(db_path, new_path)
    _fsync(new_path)
    os.rename(db_path, backup_path)
    os.rename(new_path, db_path)
    _fsync(db_path.parent)

    print("Rows per table:")
    for table_name, count in rows.items():
        print(f"  {table_name}: {count}")
    print(f"Next project ID: {next_id}")
    print(f"The original database is now at {backup_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
