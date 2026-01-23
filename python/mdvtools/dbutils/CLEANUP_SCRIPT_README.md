# Project Cleanup Script

## Overview

The `cleanup_projects.py` script provides utilities for cleaning up projects in the MDV database. It handles both database records and filesystem directories.

**⚠️ WARNING**: This script is intended primarily for local development environments. It is mostly LLM authored, as a quick convenience - but hopefully of some wider use.
Be very cautious and please review carefully before running it in production/server deployments, or if you are not confident that any important data is backed up!

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

## API Structure

The cleanup script uses an orthogonal API design where **actions** and **selectors** are separate, making it easier to compose operations.

### Actions (Mutually Exclusive)

Choose **one** action to perform:

- `--list SELECTOR`: List projects matching selector (read-only)
- `--soft-delete SELECTOR`: Mark projects as deleted (sets is_deleted=True)
- `--purge SELECTOR`: Remove projects completely (filesystem + database)
- `--restore SELECTOR`: Restore deleted projects (sets is_deleted=False)

### Selectors

Specify which projects to target:

- `deleted`: Projects where is_deleted=True
- `empty`: Projects with empty datasources (datasources.json contains `[]`)
- `orphaned`: Projects in filesystem but not in database
- `test-generated`: Projects with provenance metadata (created by test scripts)
- `all`: All projects (⚠️ requires --confirm flag for non-list actions)

### Additional Filters

- `--filter-name PATTERN`: Further narrow selection by name pattern
  - Supports `*` and `%` wildcards
  - Examples: `"pbmc*"`, `"%test%"`, `"*mock*"`

### Safety Options

- `--dry-run`: Preview changes without executing
- `--confirm`: Skip confirmation prompts (for automation, use with caution)

## Usage Examples

### List Operations (Read-Only)

```bash
# List all deleted projects
python -m mdvtools.dbutils.cleanup_projects --list deleted

# List projects with empty datasources
python -m mdvtools.dbutils.cleanup_projects --list empty

# List orphaned filesystem directories
python -m mdvtools.dbutils.cleanup_projects --list orphaned

# List test-generated projects
python -m mdvtools.dbutils.cleanup_projects --list test-generated

# List with name filtering
python -m mdvtools.dbutils.cleanup_projects --list deleted --filter-name "pbmc*"
python -m mdvtools.dbutils.cleanup_projects --list empty --filter-name "%test%"
python -m mdvtools.dbutils.cleanup_projects --list test-generated --filter-name "*mock*"

# List all projects (everything in database)
python -m mdvtools.dbutils.cleanup_projects --list all
```

### Soft-Delete Operations (Mark as Deleted)

```bash
# Soft-delete projects with empty datasources (preview)
python -m mdvtools.dbutils.cleanup_projects --soft-delete empty --dry-run

# Soft-delete projects with empty datasources (execute)
python -m mdvtools.dbutils.cleanup_projects --soft-delete empty

# Soft-delete test-generated projects without confirmation
python -m mdvtools.dbutils.cleanup_projects --soft-delete test-generated --confirm

# Soft-delete with name filtering
python -m mdvtools.dbutils.cleanup_projects --soft-delete test-generated --filter-name "pbmc*" --confirm
```

### Purge Operations (Complete Removal)

```bash
# Preview purge of deleted projects
python -m mdvtools.dbutils.cleanup_projects --purge deleted --dry-run

# Purge deleted projects (will prompt for confirmation)
python -m mdvtools.dbutils.cleanup_projects --purge deleted

# Purge projects with empty datasources
python -m mdvtools.dbutils.cleanup_projects --purge empty --confirm

# Purge orphaned directories
python -m mdvtools.dbutils.cleanup_projects --purge orphaned --confirm

# Purge with name filtering
python -m mdvtools.dbutils.cleanup_projects --purge deleted --filter-name "pbmc*" --confirm
python -m mdvtools.dbutils.cleanup_projects --purge test-generated --filter-name "*mock*" --confirm
```

### Restore Operations (Undelete)

```bash
# Restore deleted projects (preview)
python -m mdvtools.dbutils.cleanup_projects --restore deleted --dry-run

# Restore specific deleted projects by name
python -m mdvtools.dbutils.cleanup_projects --restore deleted --filter-name "important*"

# Restore all deleted projects (requires --confirm)
python -m mdvtools.dbutils.cleanup_projects --restore all --confirm
```

### Common Workflows

```bash
# Workflow 1: Clean up test data
# Step 1: See what test projects exist
python -m mdvtools.dbutils.cleanup_projects --list test-generated

# Step 2: Preview removal
python -m mdvtools.dbutils.cleanup_projects --purge test-generated --dry-run

# Step 3: Execute removal
python -m mdvtools.dbutils.cleanup_projects --purge test-generated --confirm

# Workflow 2: Handle empty datasource projects
# Step 1: List them
python -m mdvtools.dbutils.cleanup_projects --list empty

# Step 2: Soft-delete (mark as deleted)
python -m mdvtools.dbutils.cleanup_projects --soft-delete empty

# Step 3: Later, purge the deleted projects
python -m mdvtools.dbutils.cleanup_projects --purge deleted --confirm

# Workflow 3: Clean up orphaned directories
# Step 1: List orphaned projects
python -m mdvtools.dbutils.cleanup_projects --list orphaned

# Step 2: Remove them
python -m mdvtools.dbutils.cleanup_projects --purge orphaned --confirm
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

## Provenance Tracking

Projects generated by `python/mdvtools/tests/generate_test_data.py` include provenance metadata in `state.json`:

```json
{
  "all_views": [],
  "provenance": {
    "created_by": "generate_test_data.py",
    "created_at": "2026-01-22T12:00:00",
    "source": "mock",
    "parameters": {
      "n_cells": 100,
      "n_genes": 200
    }
  }
}
```

This enables:
- Identification of test-generated projects with `--list test-generated`
- Display of creation details (source, parameters, timestamp)
- Future extension for more comprehensive provenance tracking

## Development Notes

- The script uses the Flask app context from `mdvtools.dbutils.mdv_server_app`
- All database operations use SQLAlchemy ORM
- Logging uses the centralized `mdvtools.logging_config` logger
- Project validation uses `is_valid_mdv_project()` to check for required files:
  - `datasources.json`
  - `views.json`
  - `state.json`
- Name filtering supports SQL LIKE patterns (`%` for multiple chars, `_` for single char)
- Shell-style wildcards (`*`, `?`) are automatically converted to SQL patterns

