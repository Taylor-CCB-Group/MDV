# Project Cleanup Script

## Overview

The `cleanup_projects.py` script provides utilities for cleaning up projects in the MDV database. It handles both database records and filesystem directories.

**⚠️ WARNING**: This script is intended primarily for local development environments. Be very cautious about running it in production/server deployments!

## Database Relationships

The script properly handles the following database foreign key relationships:

```
Project (projects table)
├── File (files table)
│   └── project_id -> projects.id (NOT NULL)
└── UserProject (user_projects table)
    └── project_id -> projects.id (NOT NULL)
```

When deleting a project, the script:
1. First deletes all related `File` records
2. Then deletes all related `UserProject` records  
3. Finally deletes the `Project` record itself

This prevents foreign key constraint violations like:
```
psycopg2.errors.NotNullViolation: null value in column "project_id" 
of relation "files" violates not-null constraint
```

## Available Operations

### List Operations (Read-only)

- `--list-deleted`: List all projects marked as deleted (is_deleted=True)
- `--list-empty`: List projects with empty datasources (datasources.json contains `[]`)
- `--list-orphaned`: List projects that exist in filesystem but not in database

### Delete Operations (Destructive)

- `--purge-deleted`: Remove deleted projects from filesystem AND database
  - Deletes project directory
  - Deletes related File records
  - Deletes related UserProject records
  - Deletes Project record

- `--purge-empty-datasources`: Handle projects with empty datasources
  - Default: Soft-delete (sets is_deleted=True, deleted_timestamp=now())
  - With `--hard-delete`: Complete removal (filesystem + database)

- `--purge-orphaned`: Remove orphaned project directories
  - Deletes directories that exist in filesystem but not registered in database

### Safety Options

- `--dry-run`: Preview changes without executing (shows what would be done)
- `--confirm`: Skip confirmation prompts (for automation, use with caution)

## Usage Examples

### Safe Exploration

```bash
# List all deleted projects
python -m mdvtools.dbutils.cleanup_projects --list-deleted

# List projects with empty datasources
python -m mdvtools.dbutils.cleanup_projects --list-empty

# List orphaned filesystem directories
python -m mdvtools.dbutils.cleanup_projects --list-orphaned
```

### Preview Changes (Dry Run)

```bash
# See what would be purged (deleted projects)
python -m mdvtools.dbutils.cleanup_projects --purge-deleted --dry-run

# See what would be purged (empty datasource projects)
python -m mdvtools.dbutils.cleanup_projects --purge-empty-datasources --hard-delete --dry-run

# See what would be purged (orphaned directories)
python -m mdvtools.dbutils.cleanup_projects --purge-orphaned --dry-run
```

### Actual Deletion (with confirmation)

```bash
# Purge deleted projects (will prompt for confirmation)
python -m mdvtools.dbutils.cleanup_projects --purge-deleted

# Soft-delete projects with empty datasources (will prompt)
python -m mdvtools.dbutils.cleanup_projects --purge-empty-datasources

# Hard-delete projects with empty datasources (will prompt)
python -m mdvtools.dbutils.cleanup_projects --purge-empty-datasources --hard-delete

# Remove orphaned directories (will prompt)
python -m mdvtools.dbutils.cleanup_projects --purge-orphaned
```

### Automated Deletion (skip confirmation)

```bash
# For automation - skips confirmation prompts
python -m mdvtools.dbutils.cleanup_projects --purge-deleted --confirm
python -m mdvtools.dbutils.cleanup_projects --purge-orphaned --confirm
```

## Error Handling

The script includes comprehensive error handling:

- **Filesystem errors**: Catches permission errors, missing paths, etc.
- **Database errors**: Includes transaction rollback on failures
- **Partial failures**: Continues processing remaining projects even if one fails
- **Error reporting**: Collects and displays all errors at the end with details

## Example Output

```
================================================================================
DELETED PROJECTS
================================================================================
Found 23 project(s)

1. [✓ EXISTS] ID 25: unnamed_project
   Path: /app/mdv/25
   Deleted: 2026-01-16 18:23:45.123456
   Datasources: 2

...

================================================================================
SUMMARY: PURGE DELETED PROJECTS
================================================================================
Total projects processed: 23
Filesystem removals: 23
Database removals: 23

✓ All operations completed successfully
```

## Implementation Details

### delete_project_from_database()

This function handles the proper deletion order:

```python
1. Query and delete all File records (project_id foreign key)
2. Query and delete all UserProject records (project_id foreign key)  
3. Delete the Project record itself
4. Commit transaction (or rollback on error)
```

### Dry Run Mode

When `--dry-run` is used:
- No filesystem changes are made
- No database changes are made
- Operations are logged as "[DRY RUN]"
- Results show what would have been done

### Confirmation Prompts

Unless `--confirm` flag is used, destructive operations will prompt:

```
⚠️  This will PERMANENTLY DELETE projects from filesystem and database. Continue? (yes/no):
```

User must type "yes" or "y" to proceed.

## Development Notes

- The script uses the Flask app context from `mdvtools.dbutils.mdv_server_app`
- All database operations use SQLAlchemy ORM
- Logging uses the centralized `mdvtools.logging_config` logger
- Project validation uses `is_valid_mdv_project()` to check for required files:
  - `datasources.json`
  - `views.json`
  - `state.json`

