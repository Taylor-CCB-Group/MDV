#!/usr/bin/env python3
"""Rebuild a SQLite MDV database so Project IDs are never reused.

SQLite gives a deleted project's ID to the next project unless the projects table
is declared AUTOINCREMENT, and it only reads that declaration when the table is
created. A database created before the models declared it keeps reusing IDs until
it is rebuilt. This script builds a new database file from the models, copies every
row into it and swaps it in, keeping the original beside it as
<name>.pre-autoincrement.

A table or nullable column that the models declare and the database lacks is
created empty. The run stops and lists anything the database holds that the models
do not declare, which can be a table, column, index, trigger, view or non-zero
user_version. Pass --keep-extras to copy those into the new database. Rows that
already break a foreign key are copied as they are and listed.

The run also stops if a -journal, -wal or -shm file sits beside the database, or if
the database is in WAL mode. If a copied row breaks a constraint, the run stops,
removes the file it was building and leaves the original untouched.

A dry run writes nothing, so it is safe at any time inside the running app
container. It prints the rows in each table and the ID the next project will get:

    docker compose -f <compose file> exec <service> uv run python mdvtools/scripts/migrate_projects_autoincrement.py --dry-run

Stop the app before migrating, and run the script in a one-off container:

    docker compose -f <compose file> stop <service>
    docker compose -f <compose file> run --rm --no-deps <service> uv run python mdvtools/scripts/migrate_projects_autoincrement.py
    docker compose -f <compose file> start <service>

Do not run the migration with docker compose exec in the running container. The app
holds the database open, so anything it writes after the copy would be lost.

Run it as the image's default user, as the command above does, and leave out -u root.
The new database file belongs to whoever runs the script, so a run as root leaves a
database the app cannot write to.

To roll back, stop the service, move <name>.pre-autoincrement back to the live name,
and start the service.
"""

import argparse
import os
import shutil
import sqlite3
import sys
import tempfile
from collections import Counter
from contextlib import closing
from dataclasses import dataclass, field
from pathlib import Path

# Add the parent directory to sys.path to import mdvtools modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sqlalchemy.dialects import sqlite as sqlite_dialect
from sqlalchemy.schema import CreateIndex, CreateTable

from mdvtools.dbutils.dbmodels import db

DEFAULT_DATABASE = "/app/mdv/mdv.sqlite3"
BACKUP_SUFFIX = ".pre-autoincrement"


@dataclass
class _Comparison:
    """How a database differs from the models."""

    # Each model table the database holds, with the columns to copy from it.
    copy_columns: dict[str, list[str]] = field(default_factory=dict)
    missing_tables: list[str] = field(default_factory=list)
    # Missing columns are named as table.column.
    missing_columns: list[str] = field(default_factory=list)
    missing_required_columns: list[str] = field(default_factory=list)
    # Each table maps to its CREATE TABLE statement.
    extra_tables: dict[str, str] = field(default_factory=dict)
    # Each model table maps to the name and definition of every column the models lack.
    extra_columns: dict[str, list[tuple[str, str]]] = field(default_factory=dict)
    # The type, name and SQL of each index, trigger and view, in the order they were created.
    extra_objects: list[tuple[str, str, str]] = field(default_factory=list)
    user_version: int = 0

    def copied_tables(self):
        return [*self.copy_columns, *self.extra_tables]

    def extras(self):
        """One line for each thing the database holds that the models do not declare."""
        lines = [f"table {name}" for name in self.extra_tables]
        lines += [
            f"column {table}.{name}"
            for table, columns in self.extra_columns.items()
            for name, _ in columns
        ]
        lines += [f"{kind} {name}" for kind, name, _ in self.extra_objects]
        if self.user_version:
            lines.append(f"user_version {self.user_version}")
        return lines


def _quote(identifier):
    return '"' + identifier.replace('"', '""') + '"'


def _compile(element):
    return str(element.compile(dialect=sqlite_dialect.dialect()))


def _fsync(path):
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _print_list(heading, lines, file=None):
    print(heading, file=file)
    for line in lines:
        print(f"  {line}", file=file)


def _print_report(comparison, violations, rows, next_id):
    for heading, lines in (
        ("Tables the models declare and the database lacks, created empty:", comparison.missing_tables),
        ("Columns the models declare and the database lacks, left empty:", comparison.missing_columns),
        ("Copied from the database, though the models do not declare them:", comparison.extras()),
        (
            "Rows that already break a foreign key, copied as they are:",
            [f"{table} row {rowid} has no matching row in {parent}" for table, rowid, parent, _ in violations],
        ),
        ("Rows per table:", [f"{name}: {count}" for name, count in rows.items()]),
    ):
        if lines:
            _print_list(heading, lines)
    print(f"Next project ID: {next_id}")


def _compare(connection):
    """Compare the database open on connection with the models."""
    comparison = _Comparison(user_version=connection.execute("PRAGMA user_version").fetchone()[0])
    stored_tables = {
        name: sql
        for name, sql in connection.execute("SELECT name, sql FROM sqlite_master WHERE type = 'table'")
        if not name.startswith("sqlite_")
    }

    for table in db.metadata.sorted_tables:
        if table.name not in stored_tables:
            comparison.missing_tables.append(table.name)
            continue
        stored_columns = connection.execute(f"PRAGMA table_info({_quote(table.name)})").fetchall()
        stored_names = {row[1] for row in stored_columns}
        copy_columns = []
        for column in table.columns:
            if column.name in stored_names:
                copy_columns.append(column.name)
            elif column.nullable:
                comparison.missing_columns.append(f"{table.name}.{column.name}")
            else:
                comparison.missing_required_columns.append(f"{table.name}.{column.name}")

        model_names = {column.name for column in table.columns}
        for _, name, declared_type, not_null, default, _ in stored_columns:
            if name in model_names:
                continue
            definition = f"{_quote(name)} {declared_type}".rstrip()
            if not_null:
                definition += " NOT NULL"
            if default is not None:
                definition += f" DEFAULT {default}"
            comparison.extra_columns.setdefault(table.name, []).append((name, definition))
            copy_columns.append(name)
        comparison.copy_columns[table.name] = copy_columns

    model_tables = {table.name for table in db.metadata.sorted_tables}
    comparison.extra_tables = {
        name: sql for name, sql in stored_tables.items() if name not in model_tables
    }

    # SQLite stores no SQL for the indexes it builds behind UNIQUE and PRIMARY KEY
    # constraints, and those come back with their table.
    model_indexes = {index.name for table in db.metadata.sorted_tables for index in table.indexes}
    comparison.extra_objects = [
        (kind, name, sql)
        for kind, name, sql in connection.execute(
            "SELECT type, name, sql FROM sqlite_master"
            " WHERE type IN ('index', 'trigger', 'view') AND sql IS NOT NULL ORDER BY rowid"
        )
        if name not in model_indexes
    ]
    return comparison


def _build_from_models(new_path, original_path, comparison, sequence, violations):
    """Create the models' tables in new_path, copy the original's rows into them,
    set the projects sequence and check the result.

    violations are the rows PRAGMA foreign_key_check reports in the original.
    Returns the rows in each copied table.
    """
    with closing(sqlite3.connect(new_path.as_uri(), uri=True, isolation_level=None)) as connection:
        connection.execute("ATTACH DATABASE ? AS original", (original_path.as_uri() + "?mode=ro",))
        # Rows that already break a foreign key in the original are copied as they
        # are. The check after the copy compares them with the new file.
        connection.execute("PRAGMA foreign_keys = OFF")

        connection.execute("BEGIN")
        for table in db.metadata.sorted_tables:
            connection.execute(_compile(CreateTable(table)))
            for index in table.indexes:
                connection.execute(_compile(CreateIndex(index)))
        for table_name, columns in comparison.extra_columns.items():
            for _, definition in columns:
                connection.execute(f"ALTER TABLE main.{_quote(table_name)} ADD COLUMN {definition}")
        for sql in comparison.extra_tables.values():
            connection.execute(sql)

        for table_name, columns in comparison.copy_columns.items():
            names = ", ".join(_quote(name) for name in columns)
            connection.execute(
                f"INSERT INTO main.{_quote(table_name)} ({names})"
                f" SELECT {names} FROM original.{_quote(table_name)}"
            )
        for table_name in comparison.extra_tables:
            connection.execute(
                f"INSERT INTO main.{_quote(table_name)} SELECT * FROM original.{_quote(table_name)}"
            )

        # Created after the copy, so no trigger fires on the rows being copied.
        for _, _, sql in comparison.extra_objects:
            connection.execute(sql)
        if comparison.user_version:
            connection.execute(f"PRAGMA main.user_version = {int(comparison.user_version)}")

        # SQLite does not promise that inserting explicit ids moves the sequence,
        # so it is set here.
        connection.execute("DELETE FROM main.sqlite_sequence WHERE name = 'projects'")
        connection.execute(
            "INSERT INTO main.sqlite_sequence (name, seq) VALUES ('projects', ?)", (sequence,)
        )
        connection.execute("COMMIT")

        rows = {}
        for table_name in comparison.copied_tables():
            copied = connection.execute(f"SELECT count(*) FROM main.{_quote(table_name)}").fetchone()[0]
            expected = connection.execute(
                f"SELECT count(*) FROM original.{_quote(table_name)}"
            ).fetchone()[0]
            if copied != expected:
                raise RuntimeError(f"Copied {copied} of {expected} rows from {table_name}")
            rows[table_name] = copied

        # integrity_check, unlike quick_check, verifies UNIQUE constraints.
        integrity = connection.execute("PRAGMA main.integrity_check").fetchall()
        if integrity != [("ok",)]:
            raise RuntimeError(f"The new database failed its integrity check: {integrity}")

        # Counted per table and parent, because a table without an INTEGER PRIMARY
        # KEY can give copied rows new rowids.
        introduced = Counter(
            (table, parent)
            for table, _, parent, _ in connection.execute("PRAGMA main.foreign_key_check")
        ) - Counter((table, parent) for table, _, parent, _ in violations)
        if introduced:
            raise RuntimeError(
                f"The new database breaks foreign keys the original does not: {dict(introduced)}"
            )

        sequence_rows = connection.execute(
            "SELECT seq FROM main.sqlite_sequence WHERE name = 'projects'"
        ).fetchall()
        if sequence_rows != [(sequence,)]:
            raise RuntimeError(
                f"The projects sequence is {sequence_rows} after the copy, expected {sequence}"
            )

        connection.execute("DETACH DATABASE original")
    return rows


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
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the rows in each table and the next Project ID, and write nothing",
    )
    parser.add_argument(
        "--keep-extras",
        action="store_true",
        help="copy tables, columns, indexes, triggers, views and a user_version that the "
        "models do not declare into the new database",
    )
    parser.add_argument(
        "--min-next-id",
        type=int,
        help="the lowest ID the next project may get. Use it when projects with higher IDs "
        "were purged before the migration, because nothing in the database records them",
    )
    args = parser.parse_args(argv)

    db_path = Path(args.database).resolve()
    if not db_path.is_file():
        print(f"There is no database at {db_path}.", file=sys.stderr)
        return 1
    # Each of these files means a connection is open or did not close cleanly. SQLite
    # would read one left beside the new file as belonging to it.
    for suffix in ("-journal", "-wal", "-shm"):
        sidecar = db_path.with_name(db_path.name + suffix)
        if sidecar.exists():
            print(
                f"{sidecar} exists, so something may still have the database open or may not "
                "have closed it cleanly. Stop the app before migrating. If the file is still "
                "there, start and stop the app once so SQLite can recover the database.",
                file=sys.stderr,
            )
            return 1
    # Bytes 18 and 19 of the header are 2 for a database in WAL mode, where even a
    # read-only connection can create -wal and -shm files.
    with open(db_path, "rb") as database_file:
        header = database_file.read(20)
    if 2 in header[18:20]:
        print(
            f"{db_path} is in WAL mode. Switch it back with PRAGMA journal_mode = DELETE "
            "before migrating.",
            file=sys.stderr,
        )
        return 1

    with closing(sqlite3.connect(db_path.as_uri() + "?mode=ro", uri=True)) as original:
        projects = original.execute(
            "SELECT sql FROM sqlite_master WHERE type = 'table' AND name = 'projects'"
        ).fetchone()
        if projects is not None and "AUTOINCREMENT" in projects[0].upper():
            print("The projects table is already declared AUTOINCREMENT, so there is nothing to do.")
            return 0
        comparison = _compare(original)
        rows = {
            table_name: original.execute(f"SELECT count(*) FROM {_quote(table_name)}").fetchone()[0]
            for table_name in comparison.copied_tables()
        }
        highest_id = 0
        if "id" in comparison.copy_columns.get("projects", []):
            highest_id = original.execute("SELECT max(id) FROM projects").fetchone()[0] or 0
        violations = original.execute("PRAGMA foreign_key_check").fetchall()

    backup_path = db_path.with_name(db_path.name + BACKUP_SUFFIX)
    if backup_path.exists():
        print(f"{backup_path} already exists. Move it away before migrating.", file=sys.stderr)
        return 1

    if comparison.missing_required_columns:
        _print_list(
            "The models declare these columns NOT NULL, and the database has no values for them:",
            comparison.missing_required_columns,
            file=sys.stderr,
        )
        return 1
    extras = comparison.extras()
    if extras and not args.keep_extras:
        _print_list(
            "The database holds these, which the models do not declare:", extras, file=sys.stderr
        )
        print("Run with --keep-extras to copy them into the new database.", file=sys.stderr)
        return 1

    sequence = highest_id if args.min_next_id is None else max(highest_id, args.min_next_id - 1)
    next_id = sequence + 1
    if args.dry_run:
        _print_report(comparison, violations, rows, next_id)
        print("This was a dry run, so nothing was written.")
        return 0

    descriptor, new_name = tempfile.mkstemp(dir=db_path.parent, prefix=f"{db_path.name}.", suffix=".tmp")
    os.close(descriptor)
    new_path = Path(new_name)

    try:
        rows = _build_from_models(new_path, db_path, comparison, sequence, violations)
        # The new file replaces the original under the same name, so it takes the
        # original's permissions and reaches the disk before either rename.
        shutil.copymode(db_path, new_path)
        _fsync(new_path)
    except Exception as error:
        new_path.unlink(missing_ok=True)
        print(
            "The migration stopped and the original database is unchanged. "
            f"{type(error).__name__}: {error}",
            file=sys.stderr,
        )
        return 1

    os.rename(db_path, backup_path)
    os.rename(new_path, db_path)
    _fsync(db_path.parent)

    _print_report(comparison, violations, rows, next_id)
    print(f"The original database is now at {backup_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
