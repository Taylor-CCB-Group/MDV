#!/usr/bin/env python3
"""Rebuild an existing SQLite projects table so Project IDs are never reused.

SQLite hands back max(rowid) after the highest row is deleted unless the table is
declared AUTOINCREMENT, and that keyword can only be set at CREATE TABLE. A
database created before the declaration was added therefore keeps reusing IDs
until its projects table is rebuilt, which is what this script does.

The rebuild copies every row, so project names, paths and access levels are
preserved. Run it with --dry-run first to see what it will do.
"""

import argparse
import os
import shutil
import sqlite3
import sys

# Add the parent directory to sys.path to import mdvtools modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sqlalchemy.dialects import sqlite as sqlite_dialect
from sqlalchemy.schema import CreateIndex, CreateTable

from mdvtools.dbutils.dbmodels import Project

TABLE = "projects"
OLD_TABLE = "projects_pre_autoincrement"
BACKUP_SUFFIX = ".pre-autoincrement"


def _existing_table_sql(connection):
    row = connection.execute(
        "SELECT sql FROM sqlite_master WHERE type = 'table' AND name = ?", (TABLE,)
    ).fetchone()
    return row[0] if row else None


def _existing_columns(connection):
    return [row[1] for row in connection.execute(f"PRAGMA table_info({TABLE})")]


def _projects_digest(connection):
    """The fields worth checking survived the rebuild."""
    return connection.execute(f"SELECT id, name, path FROM {TABLE} ORDER BY id").fetchall()


def _compile(element):
    return str(element.compile(dialect=sqlite_dialect.dialect()))


def migrate_projects_table(db_path, dry_run=False, backup=True):
    """Rebuild the projects table with AUTOINCREMENT, keeping every row.

    Returns a summary dict. 'migrated' is False when there was nothing to do,
    which makes the script safe to run more than once.
    """
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"No such database: {db_path}")

    connection = sqlite3.connect(db_path)
    connection.isolation_level = None
    try:
        existing_sql = _existing_table_sql(connection)
        if existing_sql is None:
            raise ValueError(f"{db_path} has no '{TABLE}' table")

        if "AUTOINCREMENT" in existing_sql.upper():
            return {
                "migrated": False,
                "reason": f"'{TABLE}' is already declared AUTOINCREMENT",
                "rows": connection.execute(f"SELECT count(*) FROM {TABLE}").fetchone()[0],
            }

        leftover = connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ?", (OLD_TABLE,)
        ).fetchone()
        if leftover:
            raise ValueError(
                f"{db_path} still has a '{OLD_TABLE}' table from an interrupted run. "
                "Check which of the two holds the projects you want before continuing."
            )

        existing_columns = _existing_columns(connection)
        target_columns = [column.name for column in Project.__table__.columns]
        shared = [name for name in target_columns if name in existing_columns]
        dropped = [name for name in existing_columns if name not in target_columns]
        added = [name for name in target_columns if name not in existing_columns]

        # A column the model requires and the database has no value for would fail
        # the copy. Say so rather than leaving a half-migrated database.
        unfillable = [
            name for name in added
            if not Project.__table__.columns[name].nullable
            and Project.__table__.columns[name].server_default is None
        ]
        if unfillable:
            raise ValueError(
                f"Cannot rebuild '{TABLE}': the model requires {unfillable}, which the "
                "database has no column for and no server default to fill"
            )

        digest_before = _projects_digest(connection)
        summary = {
            "migrated": False,
            "rows": len(digest_before),
            "copied_columns": shared,
            "dropped_columns": dropped,
            "added_columns": added,
            "backup": None,
        }

        if dry_run:
            summary["reason"] = "dry run, nothing written"
            return summary

        if backup:
            backup_path = db_path + BACKUP_SUFFIX
            if os.path.exists(backup_path):
                raise FileExistsError(
                    f"{backup_path} already exists. Move or remove it so this run cannot "
                    "overwrite an earlier backup."
                )
            shutil.copy2(db_path, backup_path)
            summary["backup"] = backup_path

        # Move the old table aside and create the new one under the real name, so
        # the DDL comes straight from the model and its foreign key to genomes
        # still resolves.
        create_table = _compile(CreateTable(Project.__table__))
        columns = ", ".join(f'"{name}"' for name in shared)

        # The rename must not rewrite the foreign keys that other tables declare
        # against this one, and the drop must not cascade.
        connection.execute("PRAGMA legacy_alter_table = ON")
        connection.execute("PRAGMA foreign_keys = OFF")
        connection.execute("BEGIN")
        try:
            connection.execute(f"ALTER TABLE {TABLE} RENAME TO {OLD_TABLE}")
            connection.execute(create_table)
            connection.execute(
                f"INSERT INTO {TABLE} ({columns}) SELECT {columns} FROM {OLD_TABLE}"
            )
            copied = connection.execute(f"SELECT count(*) FROM {TABLE}").fetchone()[0]
            if copied != len(digest_before):
                raise RuntimeError(
                    f"Copied {copied} of {len(digest_before)} projects, rolling back"
                )

            # Dropping the old table takes its indexes with it, which frees their
            # names for the ones the model declares.
            connection.execute(f"DROP TABLE {OLD_TABLE}")
            for index in Project.__table__.indexes:
                connection.execute(_compile(CreateIndex(index)))

            # Start the sequence above every ID ever used, so the rows just copied
            # cannot be handed out again.
            highest = connection.execute(f"SELECT max(id) FROM {TABLE}").fetchone()[0]
            if highest is not None:
                connection.execute(
                    "INSERT OR REPLACE INTO sqlite_sequence (name, seq) VALUES (?, ?)",
                    (TABLE, highest),
                )
            connection.execute("COMMIT")
        except Exception:
            connection.execute("ROLLBACK")
            raise
        finally:
            connection.execute("PRAGMA foreign_keys = ON")
            connection.execute("PRAGMA legacy_alter_table = OFF")

        digest_after = _projects_digest(connection)
        if digest_after != digest_before:
            raise RuntimeError(
                f"Projects differ after the rebuild. The database before it is at "
                f"{summary['backup']}"
            )

        summary["migrated"] = True
        summary["reason"] = f"rebuilt '{TABLE}' with AUTOINCREMENT"
        return summary
    finally:
        connection.close()


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild an existing SQLite projects table so Project IDs are never reused"
    )
    parser.add_argument("database", help="Path to the SQLite database file")
    parser.add_argument(
        "--dry-run", action="store_true", help="Report what would change and write nothing"
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help=f"Skip copying the database to <database>{BACKUP_SUFFIX} first",
    )
    args = parser.parse_args()

    try:
        summary = migrate_projects_table(
            args.database, dry_run=args.dry_run, backup=not args.no_backup
        )
    except Exception as error:
        print(f"Error: {error}")
        sys.exit(1)

    print(f"Projects: {summary['rows']}")
    if summary.get("dropped_columns"):
        print(f"Columns in the database but not the model, not copied: {summary['dropped_columns']}")
    if summary.get("added_columns"):
        print(f"Columns in the model but not the database, left at their default: {summary['added_columns']}")
    if summary.get("backup"):
        print(f"Backup written to {summary['backup']}")
    print(summary["reason"])


if __name__ == '__main__':
    main()
