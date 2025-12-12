# MDV Database Backend Manual

## Table of Contents
1. [Overview](#overview)
2. [Database Architecture](#database-architecture)
3. [Database Models](#database-models)
4. [Database Services](#database-services)
5. [Database Operations](#database-operations)
6. [Configuration](#configuration)
7. [API Routes and Endpoints](#api-routes-and-endpoints)
8. [Setup and Initialization](#setup-and-initialization)
9. [File Synchronization](#file-synchronization)
10. [Authentication Integration](#authentication-integration)
11. [Best Practices](#best-practices)
12. [Troubleshooting](#troubleshooting)

---

## Overview

The MDV database backend provides a PostgreSQL-based persistence layer for managing projects, users, files, and related metadata. The system uses SQLAlchemy ORM for database operations and integrates with Flask for web services.

### Key Components

- **Database Models** (`dbutils/dbmodels.py`): SQLAlchemy model definitions
- **Database Services** (`dbutils/dbservice.py`): Business logic layer for database operations
- **Server Application** (`dbutils/mdv_server_app.py`): Flask app initialization and database setup
- **Routes** (`dbutils/routes.py`): API endpoints for database-backed operations
- **Project Manager Extension** (`dbutils/project_manager_extension.py`): Extended project management features

### Technology Stack

- **Database**: PostgreSQL 16
- **ORM**: SQLAlchemy (Flask-SQLAlchemy)
- **Web Framework**: Flask
- **Connection Pooling**: psycopg2 with gevent support (psycogreen)

---

## Database Architecture

### Database Connection

The system connects to PostgreSQL using connection strings in the format:
```
postgresql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}/{DB_NAME}
```

### Connection Management

- Database credentials are loaded from environment variables or Docker secrets
- Automatic database creation if the target database doesn't exist
- Connection retry logic with exponential backoff (30 retries, 5-second intervals)
- Gevent-compatible connection pooling for async operations

### Database Initialization Flow

1. **Wait for Database**: `wait_for_database()` ensures PostgreSQL is ready
2. **Create Database**: Automatically creates the database if it doesn't exist
3. **Create Tables**: `db.create_all()` creates all tables defined in models
4. **Sync Users**: If authentication is enabled, syncs users from Auth provider
5. **Cache Setup**: Initializes Redis cache for user-project mappings (if auth enabled)

---

## Database Models

All models are defined in `python/mdvtools/dbutils/dbmodels.py` using SQLAlchemy.

### User Model

**Table**: `users`

Stores user account information and authentication details.

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique user identifier |
| `email` | String(255) | Unique, Not Null | User email address |
| `auth_id` | String(255) | Unique, Not Null | Auth0 User ID |
| `password` | String(255) | Not Null | Password hash (may be empty for OAuth users) |
| `first_name` | String(50) | Not Null | User's first name |
| `last_name` | String(50) | Not Null | User's last name |
| `is_active` | Boolean | Not Null, Default: False | Account active status |
| `is_admin` | Boolean | Not Null, Default: False | Administrator flag |
| `administrator` | Boolean | Not Null, Default: False | Legacy admin flag |
| `institution` | Text | Nullable | User's institution |
| `confirmed_at` | DateTime | Nullable | Account confirmation timestamp |

**Relationships**:
- `projects`: One-to-many with `UserProject`
- `jobs`: One-to-many with `Job`
- `preferences`: One-to-many with `UserPreference`

### Project Model

**Table**: `projects`

Core model for MDV projects, tracking metadata and access control.

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique project identifier |
| `name` | String(255) | Not Null, Default: 'unnamed_project' | Project name |
| `path` | String(1024) | Not Null, Unique | Filesystem path to project |
| `created_timestamp` | DateTime | Not Null, Default: now | Project creation time |
| `update_timestamp` | DateTime | Not Null, Default: now | Last modification time |
| `accessed_timestamp` | DateTime | Not Null, Default: now | Last access time |
| `is_deleted` | Boolean | Not Null, Default: False | Soft delete flag |
| `deleted_timestamp` | DateTime | Nullable | Soft delete timestamp |
| `owner` | Integer | Nullable | Owner user ID (legacy) |
| `type` | Text | Nullable | Project type |
| `data` | JSON | Nullable | Additional project metadata |
| `is_public` | Boolean | Not Null, Default: False | Public visibility flag |
| `date_made_public` | DateTime | Nullable | When project was made public |
| `status` | Text | Nullable | Project status |
| `genome` | String | Foreign Key → `genomes.name` | Associated genome |
| `parent` | Integer | Nullable | Parent project ID |
| `description` | Text | Nullable | Project description |
| `access_level` | String(50) | Not Null, Default: 'editable' | Access level: 'editable' or 'read-only' |

**Indexes**:
- `idx_projects_genome`: Index on `genome` column
- `idx_projects_owner`: Index on `owner` column

**Relationships**:
- `users`: Many-to-many with `User` via `UserProject`
- `files`: One-to-many with `File`
- `genomes`: Many-to-one with `Genome`

### File Model

**Table**: `files`

Tracks files associated with projects.

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique file identifier |
| `name` | String(255) | Not Null | File name |
| `file_path` | String(255) | Not Null, Unique | Full filesystem path |
| `upload_timestamp` | DateTime | Not Null, Default: now | File upload time |
| `update_timestamp` | DateTime | Not Null, Default: now, On Update: now | Last modification time |
| `project_id` | Integer | Foreign Key → `projects.id`, Not Null | Associated project |

**Relationships**:
- `project`: Many-to-one with `Project`

### UserProject Model

**Table**: `user_projects`

Junction table for user-project permissions and access control.

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique relationship identifier |
| `user_id` | Integer | Foreign Key → `users.id`, Not Null | User identifier |
| `project_id` | Integer | Foreign Key → `projects.id`, Not Null | Project identifier |
| `can_read` | Boolean | Not Null, Default: False | Read permission |
| `can_write` | Boolean | Not Null, Default: False | Write permission |
| `is_owner` | Boolean | Not Null, Default: False | Ownership flag |

**Permission Logic**:
- If `is_owner` is True, both `can_read` and `can_write` are automatically True
- If `can_write` is True, `can_read` is automatically True
- Default `can_read` is True when user is added to a project

**Relationships**:
- `user`: Many-to-one with `User`
- `project`: Many-to-one with `Project`

### Genome Model

**Table**: `genomes`

Stores genome reference information.

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique genome identifier |
| `name` | String(50) | Not Null, Unique | Genome name |
| `label` | Text | Nullable | Display label |
| `data` | JSON | Nullable | Additional genome data |
| `database` | Text | Nullable | Database reference |
| `date_added` | DateTime | Not Null, Default: now | Creation timestamp |
| `connections` | Integer | Nullable | Connection count |
| `icon` | Text | Nullable | Icon URL/path |
| `small_icon` | Text | Nullable | Small icon URL/path |
| `is_public` | Boolean | Not Null, Default: True | Public visibility |
| `chrom_sizes` | JSON | Nullable | Chromosome sizes data |

**Relationships**:
- `projects`: One-to-many with `Project`

### Job Model

**Table**: `jobs`

Tracks background jobs and processing tasks.

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique job identifier |
| `user_id` | Integer | Foreign Key → `users.id`, Nullable | Associated user |
| `inputs` | JSON | Nullable | Job input parameters |
| `outputs` | JSON | Nullable | Job output results |
| `sent_on` | DateTime | Not Null, Default: now | Job submission time |
| `finished_on` | DateTime | Nullable | Job completion time |
| `status` | String(200) | Nullable | Job status |
| `class_name` | String(200) | Nullable | Job class identifier |
| `genome` | String(100) | Nullable | Associated genome |
| `type` | String(200) | Nullable | Job type |
| `is_deleted` | Boolean | Not Null, Default: False | Soft delete flag |

**Relationships**:
- `user`: Many-to-one with `User`

### SharedObject Model

**Table**: `shared_objects`

Tracks shared objects between users (legacy/optional feature).

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique share identifier |
| `owner` | Integer | Nullable | Owner user ID |
| `shared_with` | Integer | Nullable | Recipient user ID |
| `object_id` | Integer | Nullable | Shared object ID |
| `date_shared` | DateTime | Not Null, Default: now | Share timestamp |
| `level` | Text | Not Null, Default: 'view' | Share permission level |

### UserPreference Model

**Table**: `user_preferences`

Stores user-specific preferences and settings.

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique preference identifier |
| `user_id` | Integer | Foreign Key → `users.id`, Nullable | Associated user |
| `preference` | Text | Not Null | Preference key |
| `data` | JSON | Nullable | Preference value/data |

**Relationships**:
- `user`: Many-to-one with `User`

### ViewSet Model

**Table**: `view_sets`

Stores view configurations and presets.

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique view set identifier |
| `table_name` | String(100) | Nullable | Associated table name |
| `name` | String(200) | Nullable | View set name |
| `description` | Text | Nullable | Description |
| `date_added` | DateTime | Not Null, Default: now | Creation timestamp |
| `date_modified` | DateTime | Not Null, Default: now, On Update: now | Last modification |
| `owner` | Integer | Default: 0 | Owner user ID |
| `is_public` | Boolean | Default: False | Public visibility |
| `fields` | JSON | Nullable | View field configuration |
| `data` | JSON | Nullable | Additional view data |
| `date_made_public` | DateTime | Nullable | Publication timestamp |
| `status` | Text | Nullable | View set status |
| `is_deleted` | Boolean | Default: False | Soft delete flag |

**Indexes**:
- `idx_views_table_name`: Index on `table_name`
- `idx_views_name`: Index on `name`

### GeneSet Model

**Table**: `gene_sets`

Stores gene set definitions and annotations.

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `id` | Integer | Primary Key, Auto-increment | Unique gene set identifier |
| `name` | String(100) | Not Null | Gene set name |
| `table_name` | String(150) | Nullable | Associated table name |
| `data` | JSON | Nullable | Gene set data |
| `date_added` | DateTime | Not Null, Default: now | Creation timestamp |
| `date_modified` | DateTime | Not Null, Default: now, On Update: now | Last modification |
| `is_deleted` | Boolean | Default: False | Soft delete flag |
| `description` | Text | Nullable | Description |

**Indexes**:
- `idx_genes_name`: Index on `name`

---

## Database Services

Database services (`dbutils/dbservice.py`) provide a business logic layer for database operations. All services use static methods and handle transaction management.

### ProjectService

Manages project-related database operations.

#### `get_active_projects()`
Returns a list of all non-deleted projects with metadata.

**Returns**: List of dictionaries containing:
- `id`: Project ID
- `name`: Project name
- `lastModified`: Formatted update timestamp
- `thumbnail`: Project thumbnail image (if available)
- `readme`: Project README content (if available)

**Filters**: Excludes projects in `failed_projects` list and soft-deleted projects.

#### `get_next_project_id()`
Calculates the next available project ID.

**Returns**: Integer (next project ID)

**Logic**: Gets maximum existing ID and increments by 1, or returns 1 if no projects exist.

#### `add_new_project(path, name='unnamed_project')`
Creates a new project record in the database.

**Parameters**:
- `path` (str): Filesystem path to project directory
- `name` (str, optional): Project name (default: 'unnamed_project')

**Returns**: `Project` model instance

**Transaction**: Commits on success, rolls back on error

#### `get_project_by_id(id)`
Retrieves a project by its ID.

**Parameters**:
- `id` (int): Project ID

**Returns**: `Project` model instance or None

#### `soft_delete_project(id)`
Performs a soft delete on a project.

**Parameters**:
- `id` (int): Project ID

**Returns**: Boolean (True if successful, False if project not found)

**Side Effects**:
- Sets `is_deleted` to True
- Sets `deleted_timestamp` to current time

#### `update_project_name(project_id, new_name)`
Updates a project's name.

**Parameters**:
- `project_id` (int): Project ID
- `new_name` (str): New project name

**Returns**: Boolean (True if successful, False if project not found)

**Side Effects**:
- Updates `name`
- Updates `update_timestamp` and `accessed_timestamp`

#### `change_project_access(project_id, new_access_level)`
Changes a project's access level.

**Parameters**:
- `project_id` (int): Project ID
- `new_access_level` (str): 'editable' or 'read-only'

**Returns**: Tuple `(access_level, message, status_code)`

**Side Effects**:
- Updates `access_level`
- Updates `update_timestamp` and `accessed_timestamp`

#### `set_project_update_timestamp(project_id)`
Updates the project's modification timestamp.

**Parameters**:
- `project_id` (str/int): Project ID

**Side Effects**: Updates `update_timestamp` to current time

#### `set_project_accessed_timestamp(project_id)`
Updates the project's access timestamp.

**Parameters**:
- `project_id` (str/int): Project ID

**Side Effects**: Updates `accessed_timestamp` to current time

**Class Variable**:
- `failed_projects`: List of tuples `(project_id, error)` for projects that failed to serve

### FileService

Manages file-related database operations.

#### `add_or_update_file_in_project(file_name, file_path, project_id)`
Adds a new file record or updates an existing one.

**Parameters**:
- `file_name` (str): File name
- `file_path` (str): Full filesystem path
- `project_id` (int): Associated project ID

**Logic**:
- Checks if file exists by path and project ID
- Updates existing record if found
- Creates new record if not found

**Transaction**: Commits on success, rolls back on error

#### `get_file_by_path_and_project(file_path, project_id)`
Retrieves a file by its path and project ID.

**Parameters**:
- `file_path` (str): Full filesystem path
- `project_id` (int): Project ID

**Returns**: `File` model instance or None

#### `file_exists_in_project(file_path, project_id)`
Checks if a file exists in the database for a project.

**Parameters**:
- `file_path` (str): Full filesystem path
- `project_id` (int): Project ID

**Returns**: Boolean

#### `get_files_by_project(project_id)`
Retrieves all files for a project.

**Parameters**:
- `project_id` (int): Project ID

**Returns**: List of `File` model instances

#### `delete_files_by_project(project_id)`
Deletes all file records for a project.

**Parameters**:
- `project_id` (int): Project ID

**Returns**: Boolean (True if successful)

**Transaction**: Commits on success, rolls back on error

#### `update_file_timestamp(file_id)`
Updates a file's modification timestamp.

**Parameters**:
- `file_id` (int): File ID

**Returns**: Boolean (True if successful, False if file not found)

**Side Effects**: Updates `update_timestamp` to current time

### UserService

Manages user-related database operations.

#### `add_or_update_user(email, auth_id=None, first_name=None, last_name=None, institution=None)`
Adds a new user or updates an existing user.

**Parameters**:
- `email` (str, required): User email address
- `auth_id` (str, optional): Auth0 user ID
- `first_name` (str, optional): User's first name
- `last_name` (str, optional): User's last name
- `institution` (str, optional): User's institution

**Returns**: `User` model instance

**Logic**:
- If user exists (by email), updates provided fields
- If user doesn't exist, creates new user with defaults:
  - `is_active`: True
  - `confirmed_at`: Current timestamp
  - `password`: Empty string
  - `administrator`/`is_admin`: False

**Transaction**: Commits on success, rolls back on error

### UserProjectService

Manages user-project permission relationships.

#### `add_or_update_user_project(user_id, project_id, is_owner=False, can_write=False)`
Creates or updates a user-project permission relationship.

**Parameters**:
- `user_id` (int): User ID
- `project_id` (int): Project ID
- `is_owner` (bool, default: False): Ownership flag
- `can_write` (bool, default: False): Write permission

**Permission Logic**:
- If `is_owner` is True: `can_read` and `can_write` are set to True
- If `can_write` is True: `can_read` is set to True
- Otherwise: `can_read` defaults to True

**Returns**: `UserProject` model instance

**Transaction**: Commits on success, rolls back on error

#### `get_user_project_permissions(user_id, project_id)`
Retrieves permission information for a user-project relationship.

**Parameters**:
- `user_id` (int): User ID
- `project_id` (int): Project ID

**Returns**: Dictionary with keys:
- `can_read` (bool)
- `can_write` (bool)
- `is_owner` (bool)

**Default**: Returns all False if relationship doesn't exist

#### `remove_user_from_project(user_id, project_id)`
Removes a user from a project (deletes the relationship).

**Parameters**:
- `user_id` (int): User ID
- `project_id` (int): Project ID

**Returns**: None

**Transaction**: Commits on success, rolls back on error

---

## Database Operations

### Project Lifecycle Operations

#### Creating Projects

1. **From Filesystem** (`serve_projects_from_filesystem`):
   - Scans filesystem for projects not in database
   - Creates database records for new projects
   - Serves projects via MDVProject
   - Optionally syncs files if `ENABLE_FILE_SYNC` is enabled

2. **Via API** (`/create_project`):
   - Generates next project ID
   - Creates project directory
   - Initializes MDVProject
   - Creates database record
   - Assigns ownership (if auth enabled)
   - Updates cache

#### Serving Projects

**From Database** (`serve_projects_from_db`):
- Queries all non-deleted projects
- Validates project paths exist
- Serves projects via MDVProject
- Syncs access levels from state.json
- Optionally syncs files if `ENABLE_FILE_SYNC` is enabled
- Tracks failed projects in `ProjectService.failed_projects`

#### Updating Projects

- **Name Update**: Updates `name`, `update_timestamp`, `accessed_timestamp`
- **Access Level Update**: Updates `access_level`, timestamps
- **Timestamp Updates**: Automatic on file changes (if file sync enabled)

#### Deleting Projects

- **Soft Delete**: Sets `is_deleted=True`, `deleted_timestamp=now()`
- Projects are filtered out from active project queries
- Filesystem data is preserved

### File Synchronization

When `ENABLE_FILE_SYNC` is enabled, the system synchronizes filesystem files with the database:

1. **During Project Serve**: Walks project directory tree
2. **For Each File**: Calls `FileService.add_or_update_file_in_project()`
3. **Update Logic**: Updates existing records or creates new ones

**Configuration**:
- Environment variable: `ENABLE_FILE_SYNC` (set to "1", "true", or "yes")
- Config file: `config.json` → `enable_file_sync` (boolean)

### Timestamp Management

**Automatic Updates**:
- `update_timestamp`: Updated on project modifications (name, access level, etc.)
- `accessed_timestamp`: Updated on project access
- File `update_timestamp`: Updated on file modifications

**Manual Updates**:
- `ProjectService.set_project_update_timestamp()`
- `ProjectService.set_project_accessed_timestamp()`
- `FileService.update_file_timestamp()`

### Permission Management

**Permission Hierarchy**:
1. **Owner** (`is_owner=True`): Full access (read + write)
2. **Editor** (`can_write=True`): Read + write access
3. **Viewer** (`can_read=True`): Read-only access

**Permission Assignment**:
- New projects: Creator becomes owner
- Imported projects: Importer becomes owner
- Shared projects: Owner can grant permissions via API

**Permission Checks**:
- API endpoints verify permissions before operations
- Cache-based lookups for performance (if auth enabled)

---

## Configuration

### Configuration File

Location: `python/mdvtools/dbutils/config.json`

```json
{
  "database_uri": "postgresql://admin:password@database/mydatabase",
  "db_container": "mdv_db",
  "db_track_modifications": false,
  "upload_folder": "/python",
  "projects_base_dir": "/app/mdv",
  "extensions": ["chat", "project_manager"]
}
```

**Fields**:
- `database_uri`: Default database connection string (overridden by environment)
- `db_container`: Database container name (legacy)
- `db_track_modifications`: SQLAlchemy track modifications flag
- `upload_folder`: File upload directory
- `projects_base_dir`: Base directory for project storage
- `extensions`: List of enabled extensions

### Environment Variables

**Database Configuration**:
- `DB_USER`: Database username
- `DB_PASSWORD`: Database password
- `DB_NAME`: Database name
- `DB_HOST`: Database hostname

**Authentication**:
- `ENABLE_AUTH`: Enable authentication ("1", "true", "yes")
- `DEFAULT_AUTH_METHOD`: Authentication method (e.g., "auth0")
- `AUTH0_DOMAIN`: Auth0 domain
- `AUTH0_CLIENT_ID`: Auth0 client ID
- `AUTH0_CLIENT_SECRET`: Auth0 client secret
- `AUTH0_CALLBACK_URL`: Auth0 callback URL
- `AUTH0_PUBLIC_KEY_URI`: Auth0 public key URI
- `AUTH0_AUDIENCE`: Auth0 audience
- `AUTH0_DB_CONNECTION`: Auth0 database connection name
- `FLASK_SECRET_KEY`: Flask session secret key

**File Synchronization**:
- `ENABLE_FILE_SYNC`: Enable file synchronization ("1", "true", "yes")

**Extensions**:
- `MDV_USER_CONFIG_PATH`: Path to user-provided config file for extensions

**API Configuration**:
- `MDV_API_ROOT`: API root path prefix

### Docker Secrets

The system can read secrets from Docker secrets directory (`/run/secrets/`):

- `db_user`: Database username
- `db_password`: Database password
- `db_name`: Database name
- `flask_secret_key`: Flask session secret
- `auth0_client_secret`: Auth0 client secret

**Priority**: Environment variables override Docker secrets.

---

## API Routes and Endpoints

### Base Routes (`routes.py`)

#### `GET /`
Renders the main index page.

#### `GET /api_root`
Returns the MDV API root path.

**Response**:
```json
{
  "mdv_api_root": "/"
}
```

#### `GET /enable_auth`
Returns authentication status.

**Response**:
```json
{
  "enable_auth": true
}
```

#### `GET /extension_config`
Returns extension configuration flags.

**Response**:
```json
{
  "project_manager": {
    "createProject": true,
    "importProject": true,
    "exportProject": true,
    "deleteProject": true,
    "renameProject": true,
    "changeProjectAccess": true,
    "shareProject": true,
    "editUserPermissions": true,
    "removeUserFromProject": true
  }
}
```

#### `GET /projects`
Returns list of active projects (filtered by user permissions if auth enabled).

**Response**:
```json
[
  {
    "id": 1,
    "name": "My Project",
    "lastModified": "2024-01-01 12:00:00",
    "thumbnail": "path/to/thumbnail.png",
    "permissions": {
      "can_read": true,
      "can_write": true,
      "is_owner": true
    },
    "owner": ["user@example.com"],
    "readme": "# Project README"
  }
]
```

**Authentication**: Required if `ENABLE_AUTH` is true

#### `GET /rescan_projects`
Rescans filesystem for new projects and adds them to database.

**Authentication**: Required (admin only)

**Response**: Redirects to index or `MDV_API_ROOT`

### Project Management Routes (`project_manager_extension.py`)

#### `POST /create_project`
Creates a new project.

**Response**:
```json
{
  "id": 1,
  "name": "unnamed_project",
  "status": "success"
}
```

**Authentication**: Required if `ENABLE_AUTH` is true

**Side Effects**:
- Creates project directory
- Creates database record
- Assigns ownership to creator
- Updates cache

#### `POST /import_project`
Imports a project from a ZIP archive.

**Request**: Multipart form data
- `file`: ZIP file containing project
- `name` (optional): Project name

**Response**:
```json
{
  "id": 1,
  "name": "Imported Project",
  "status": "success"
}
```

**Authentication**: Required if `ENABLE_AUTH` is true

**Validation**:
- Checks for required files: `views.json`, `state.json`, `datasources.json`
- Rejects unsafe paths (absolute paths, "..")

#### `GET /export_project/<project_id>`
Exports a project as a ZIP archive.

**Parameters**:
- `project_id` (int): Project ID

**Response**: ZIP file download

**Authentication**: Required (owner only)

**Access Control**: Project must be editable

#### `DELETE /delete_project/<project_id>`
Soft deletes a project.

**Parameters**:
- `project_id` (int): Project ID

**Response**:
```json
{
  "message": "Project '1' deleted successfully."
}
```

**Authentication**: Required (owner only)

**Access Control**: Project must be editable

#### `PUT /projects/<project_id>/rename`
Renames a project.

**Request**: Multipart form data
- `name`: New project name

**Response**:
```json
{
  "status": "success",
  "id": 1,
  "new_name": "New Name"
}
```

**Authentication**: Required (owner only)

**Access Control**: Project must be editable

#### `PUT /projects/<project_id>/access`
Changes project access level.

**Request**: Multipart form data
- `type`: "editable" or "read-only"

**Response**:
```json
{
  "status": "success",
  "access_level": "read-only"
}
```

**Authentication**: Required (owner only)

#### `GET /projects/<project_id>/share`
Gets list of users with project access.

**Response**:
```json
{
  "shared_users": [
    {
      "id": 1,
      "email": "user@example.com",
      "permission": "Owner"
    }
  ],
  "all_users": [
    {
      "id": 2,
      "email": "other@example.com"
    }
  ]
}
```

**Authentication**: Required (owner only)

#### `POST /projects/<project_id>/share`
Adds a user to a project with specified permissions.

**Request**: JSON
```json
{
  "user_id": 2,
  "permission": "edit"
}
```

**Permissions**: "view", "edit", or "owner"

**Response**:
```json
{
  "message": "User 2 added to project 1 with edit permission."
}
```

**Authentication**: Required (owner only)

#### `POST /projects/<project_id>/share/<user_id>/edit`
Updates user permissions for a project.

**Request**: JSON
```json
{
  "permission": "view"
}
```

**Response**:
```json
{
  "message": "Permissions updated successfully"
}
```

**Authentication**: Required (owner only)

#### `POST /projects/<project_id>/share/<user_id>/delete`
Removes a user from a project.

**Response**:
```json
{
  "message": "User removed successfully"
}
```

**Authentication**: Required (owner only)

---

## Setup and Initialization

### Database Setup

1. **Start PostgreSQL**:
   ```bash
   docker-compose up -d mdv_db
   ```

2. **Configure Environment**:
   Create `.env` file:
   ```
   DB_NAME=mdv_db
   DB_USER=admin
   DB_PASSWORD=password
   DB_HOST=mdv_db
   ```

3. **Initialize Application**:
   The application automatically:
   - Waits for database to be ready
   - Creates database if it doesn't exist
   - Creates all tables
   - Syncs users (if auth enabled)

### Application Initialization Flow

1. **Load Configuration** (`load_config`):
   - Reads `config.json`
   - Loads environment variables
   - Reads Docker secrets (if available)
   - Configures database URI

2. **Initialize Database**:
   - Calls `wait_for_database()` (retries up to 30 times)
   - Creates database if needed
   - Initializes SQLAlchemy
   - Creates tables via `db.create_all()`

3. **Sync Users** (if auth enabled):
   - Syncs users from Auth provider
   - Caches user-project mappings

4. **Serve Projects**:
   - `serve_projects_from_db()`: Serves existing projects
   - `serve_projects_from_filesystem()`: Adds new projects from filesystem

5. **Register Routes**:
   - Base routes (`routes.py`)
   - Project management routes (`project_manager_extension.py`)
   - Authentication routes (if enabled)

### Database Connection Retry Logic

The `wait_for_database()` function implements retry logic:

- **Max Retries**: 30 attempts
- **Delay**: 5 seconds between retries
- **Operations**:
  1. Connect to PostgreSQL database
  2. Check if target database exists
  3. Create database if needed
  4. Test connection to target database

### Table Creation

Tables are created automatically on first run using `db.create_all()`. The system checks if tables exist before creating them to avoid errors on subsequent runs.

**Note**: The system does not use database migrations (Alembic). Schema changes require manual database updates or recreation.

---

## File Synchronization

### Overview

File synchronization maintains a database record of all files within project directories. This enables:
- File tracking and metadata
- Project file listings
- File change detection

### Configuration

**Enable File Sync**:
- Environment variable: `ENABLE_FILE_SYNC=1`
- Config file: `"enable_file_sync": true`

### Synchronization Process

1. **Trigger**: During project serve operations
2. **Scope**: All files in project directory tree
3. **Operation**: For each file:
   - Check if file exists in database (by path and project ID)
   - Update existing record or create new record
   - Update timestamps

### File Service Operations

- `add_or_update_file_in_project()`: Main sync operation
- `get_file_by_path_and_project()`: Lookup by path
- `get_files_by_project()`: List all project files
- `delete_files_by_project()`: Remove all file records
- `update_file_timestamp()`: Update modification time

### Performance Considerations

- File sync can be expensive for large projects
- Disable if not needed: `ENABLE_FILE_SYNC=false`
- Sync runs during project serve, not on every file access

---

## Authentication Integration

### User Synchronization

When authentication is enabled, users are synchronized from the Auth provider (e.g., Auth0) to the database:

1. **On Startup**: `auth_provider.sync_users_to_db()`
2. **After Project Creation**: Users are synced to ensure database consistency
3. **User Service**: `UserService.add_or_update_user()` handles user creation/updates

### Permission Caching

User-project permissions are cached in Redis (if available) for performance:

- **Cache Structure**: `user_project_cache[user_id][project_id] = {can_read, can_write, is_owner}`
- **Initialization**: `cache_user_projects()` on startup
- **Updates**: Cache updated when permissions change

### Access Control

API endpoints check permissions before operations:

1. **Owner Checks**: Verify `is_owner` flag
2. **Write Checks**: Verify `can_write` for modifications
3. **Read Checks**: Verify `can_read` for access

### Session Management

- Flask sessions store user information
- Session cookies: HTTPOnly, Secure, SameSite=Lax
- User data: `session['user']` contains `{id, email, auth_id, is_admin}`

---

## Best Practices

### Database Operations

1. **Always Use Services**: Use service layer methods instead of direct model queries
2. **Transaction Management**: Services handle commits/rollbacks automatically
3. **Error Handling**: Services log errors and handle rollbacks
4. **Soft Deletes**: Use soft deletes for projects to preserve data

### Project Management

1. **Path Uniqueness**: Ensure project paths are unique
2. **ID Generation**: Use `get_next_project_id()` for new projects
3. **Timestamp Updates**: Update timestamps on modifications
4. **Access Level Sync**: Sync access levels from `state.json` when serving projects

### File Synchronization

1. **Disable if Unused**: Set `ENABLE_FILE_SYNC=false` if file tracking isn't needed
2. **Monitor Performance**: File sync can be slow for large projects
3. **Error Handling**: File sync errors are logged but don't block project serving

### Authentication

1. **User Sync**: Ensure users are synced before assigning permissions
2. **Cache Updates**: Update permission cache after permission changes
3. **Permission Checks**: Always verify permissions in API endpoints

### Configuration

1. **Environment Variables**: Prefer environment variables over config file for secrets
2. **Docker Secrets**: Use Docker secrets in production
3. **External Config**: Use `MDV_USER_CONFIG_PATH` for deployment-specific extensions

### Error Handling

1. **Logging**: All errors are logged with context
2. **Rollbacks**: Database operations roll back on errors
3. **Failed Projects**: Track failed projects in `ProjectService.failed_projects`
4. **Graceful Degradation**: Continue serving other projects if one fails

---

## Troubleshooting

### Database Connection Issues

**Problem**: Database connection fails on startup

**Solutions**:
1. Verify environment variables are set correctly
2. Check database container is running: `docker ps`
3. Verify network connectivity between containers
4. Check database logs for errors
5. Increase retry count/delay in `wait_for_database()`

### Table Creation Errors

**Problem**: Tables fail to create

**Solutions**:
1. Check database permissions
2. Verify database exists
3. Check for existing tables with conflicting names
4. Review SQLAlchemy model definitions

### Project Serving Failures

**Problem**: Projects fail to serve

**Solutions**:
1. Check `ProjectService.failed_projects` for error details
2. Verify project paths exist in filesystem
3. Check project directory contains required files
4. Review application logs for specific errors

### Permission Issues

**Problem**: Users can't access projects

**Solutions**:
1. Verify user exists in database
2. Check `UserProject` records exist
3. Verify permission cache is updated
4. Check authentication is properly configured
5. Review API endpoint permission checks

### File Sync Issues

**Problem**: Files not syncing to database

**Solutions**:
1. Verify `ENABLE_FILE_SYNC` is enabled
2. Check file paths are valid
3. Review file service error logs
4. Verify project ID is correct
5. Check database connection is active

### Performance Issues

**Problem**: Slow database operations

**Solutions**:
1. Disable file sync if not needed
2. Check database indexes are created
3. Review query patterns for optimization
4. Consider connection pooling settings
5. Monitor database query performance

### Migration Issues

**Problem**: Schema changes not applied

**Solutions**:
1. The system uses `db.create_all()` (no migrations)
2. Manual schema updates may be required
3. Consider database backup before changes
4. Test schema changes in development first

---

## Appendix

### Database Schema Diagram

```
users
  ├── user_projects (many-to-many)
  │   └── projects
  ├── jobs (one-to-many)
  └── user_preferences (one-to-many)

projects
  ├── files (one-to-many)
  ├── user_projects (many-to-many)
  │   └── users
  └── genomes (many-to-one)

genomes
  └── projects (one-to-many)
```

### Key Files Reference

- **Models**: `python/mdvtools/dbutils/dbmodels.py`
- **Services**: `python/mdvtools/dbutils/dbservice.py`
- **Server App**: `python/mdvtools/dbutils/mdv_server_app.py`
- **Routes**: `python/mdvtools/dbutils/routes.py`
- **Project Manager**: `python/mdvtools/dbutils/project_manager_extension.py`
- **Config**: `python/mdvtools/dbutils/config.json`

### Related Documentation

- Authentication: See `docs/changes/PERFORMANCE_OPTIMIZATION_AUTH.md`
- Project Router: See `docs/extradocs/views.md`
- Server Extensions: See extension-specific documentation

---

**Last Updated**: 2024
**Version**: 1.0

