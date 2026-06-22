# Admin Portal API Contract

Status: Draft for team review

Related architecture doc: [admin-portal-architecture.md](./admin-portal-architecture.md)

## Summary

All Admin Portal APIs live under:

```text
/admin/api/*
```

All responses are JSON.

The admin frontend should use these APIs instead of importing catalog, project-manager, or database internals.

## Shared Types

### Admin User

```ts
type AdminUser = {
  id: number;
  email: string;
  firstName: string;
  lastName: string;
  isActive: boolean;
  isAdmin: boolean;
  authId?: string;
};
```

### Admin Project

```ts
type AdminProject = {
  id: number;
  name: string;
  path: string;
  accessLevel: string;
  isPublic: boolean;
  isDeleted: boolean;
  updatedAt: string | null;
};
```

### Project Permission

```ts
type ProjectPermission = "view" | "edit" | "owner";
```

Permission mapping:

| Permission | canRead | canWrite | isOwner |
| --- | --- | --- | --- |
| `view` | `true` | `false` | `false` |
| `edit` | `true` | `true` | `false` |
| `owner` | `true` | `true` | `true` |

### Project Access Input

```ts
type ProjectAccessInput = {
  projectId: number;
  permission: ProjectPermission;
};
```

### Admin Project User

```ts
type AdminProjectUser = {
  userId: number;
  projectId: number;
  email: string;
  permission: ProjectPermission;
  canRead: boolean;
  canWrite: boolean;
  isOwner: boolean;
};
```

## Error Shape

Standard error response:

```json
{
  "error": "Human-readable message"
}
```

## Authentication Behavior

Production behavior:

| Case | Response |
| --- | --- |
| No session | `401` |
| Logged in but not admin | `403` |
| Admin user | request allowed |

POC/local behavior:

```text
When auth is disabled, the backend may use a synthetic local admin.
```

## GET `/admin/api/session`

Purpose:

```text
Return the current admin/session status.
```

Response:

```json
{
  "authEnabled": true,
  "isAdmin": true,
  "permissions": [
    "admin:access",
    "admin:users:manage",
    "admin:projects:manage"
  ],
  "user": {
    "id": 1,
    "email": "admin@example.com",
    "is_admin": true
  }
}
```

POC local response:

```json
{
  "authEnabled": false,
  "isAdmin": true,
  "permissions": [
    "admin:access",
    "admin:users:manage",
    "admin:projects:manage"
  ],
  "user": {
    "id": 0,
    "email": "local-admin@localhost",
    "is_admin": true,
    "synthetic": true
  }
}
```

Errors:

| Status | Meaning |
| --- | --- |
| `401` | no authenticated session |
| `403` | authenticated user is not an admin |

## GET `/admin/api/users`

Purpose:

```text
List users known to this MDV deployment.
```

Response:

```json
{
  "users": [
    {
      "id": 1,
      "email": "user@example.com",
      "firstName": "Jane",
      "lastName": "User",
      "isActive": true,
      "isAdmin": false,
      "authId": "auth0|abc123"
    }
  ]
}
```

POC note:

```text
Without Auth0, this returns local/dev MDV users.
For demo seeding details, see admin-portal-poc-plan.md.
```

## POST `/admin/api/users`

Purpose:

```text
Create/add a user for this deployment and assign at least one project permission.
```

Request:

```json
{
  "email": "new.user@example.com",
  "firstName": "New",
  "lastName": "User",
  "projectAccess": [
    {
      "projectId": 10,
      "permission": "view"
    }
  ]
}
```

Validation:

- `email` is required.
- `projectAccess` must contain at least one project.
- every `projectId` must exist in this deployment.
- every `permission` must be one of `view`, `edit`, or `owner`.

Production behavior:

```text
1. Create user in Auth0 for this deployment, or resolve existing Auth0 user.
2. Sync/create local MDV User row.
3. Create UserProject permission rows.
4. Return user and project access.
```

Atomicity:

```text
This operation should be atomic from the API caller's point of view.
If Auth0 user creation succeeds but MDV sync or permission creation fails, the backend should roll back the Auth0-created user where possible.
If rollback fails, the response should make the partial failure explicit.
```

POC/local behavior:

```text
1. Do not call Auth0.
2. Create or resolve local/dev MDV user.
3. Create UserProject permission rows.
4. Return the same response shape.
```

Response:

```json
{
  "user": {
    "id": 5,
    "email": "new.user@example.com",
    "firstName": "New",
    "lastName": "User",
    "isActive": true,
    "isAdmin": false,
    "authId": "auth0|abc123"
  },
  "projectUsers": [
    {
      "userId": 5,
      "projectId": 10,
      "email": "new.user@example.com",
      "permission": "view",
      "canRead": true,
      "canWrite": false,
      "isOwner": false
    }
  ]
}
```

Errors:

| Status | Meaning |
| --- | --- |
| `400` | invalid email, empty project access, invalid project ID, or invalid permission |
| `409` | user already exists and conflict handling is required |
| `502` | Auth0 operation failed |
| `500` | MDV sync or permission write failed; rollback should be attempted if Auth0 user creation already happened |

Open decision:

```text
If the Auth0 user already exists, should this endpoint return the existing user or return 409?
```

Recommended behavior:

```text
Return the existing/synced user when possible and continue assigning project access.
```

## GET `/admin/api/projects`

Purpose:

```text
List projects in this MDV deployment.
```

Response:

```json
{
  "projects": [
    {
      "id": 10,
      "name": "Example Project",
      "path": "/app/mdv/example-project",
      "accessLevel": "editable",
      "isPublic": false,
      "isDeleted": false,
      "updatedAt": "2026-06-18T12:00:00"
    }
  ]
}
```

Open decision:

```text
Should deleted or archived projects be returned by default?
```

## GET `/admin/api/projects/:projectId/users`

Purpose:

```text
Return users assigned to a project and users available to add.
```

Response:

```json
{
  "project": {
    "id": 10,
    "name": "Example Project",
    "path": "/app/mdv/example-project",
    "accessLevel": "editable",
    "isPublic": false,
    "isDeleted": false,
    "updatedAt": "2026-06-18T12:00:00"
  },
  "projectUsers": [
    {
      "userId": 1,
      "projectId": 10,
      "email": "owner@example.com",
      "permission": "owner",
      "canRead": true,
      "canWrite": true,
      "isOwner": true
    }
  ],
  "availableUsers": [
    {
      "id": 2,
      "email": "viewer@example.com",
      "firstName": "",
      "lastName": "",
      "isActive": true,
      "isAdmin": false
    }
  ]
}
```

Errors:

| Status | Meaning |
| --- | --- |
| `404` | project not found |

## POST `/admin/api/projects/:projectId/users`

Purpose:

```text
Grant an existing MDV user access to another project.
```

This endpoint is for users who already exist in the deployment.

Request:

```json
{
  "userId": 2,
  "permission": "edit"
}
```

Response:

```json
{
  "projectUser": {
    "userId": 2,
    "projectId": 10,
    "email": "viewer@example.com",
    "permission": "edit",
    "canRead": true,
    "canWrite": true,
    "isOwner": false
  }
}
```

Errors:

| Status | Meaning |
| --- | --- |
| `400` | invalid user ID or permission |
| `404` | project not found |
| `404` | user not found |

## PATCH `/admin/api/projects/:projectId/users/:userId`

Purpose:

```text
Change a user's permission on a project.
```

Request:

```json
{
  "permission": "owner"
}
```

Response:

```json
{
  "projectUser": {
    "userId": 2,
    "projectId": 10,
    "email": "viewer@example.com",
    "permission": "owner",
    "canRead": true,
    "canWrite": true,
    "isOwner": true
  }
}
```

Errors:

| Status | Meaning |
| --- | --- |
| `400` | invalid permission |
| `404` | project not found |
| `404` | user not found |
| `404` | project membership not found |
| `409` | operation would leave project without an owner |

MVP rule:

```text
Admin cannot demote the last owner of a project.
```

## DELETE `/admin/api/projects/:projectId/users/:userId`

Purpose:

```text
Remove a user's access to a project.
```

MVP meaning:

```text
This removes project access only.
It does not delete the user from Auth0.
It does not remove the user from the deployment.
It does not deactivate the user globally.
```

Response:

```json
{
  "removed": true
}
```

Errors:

| Status | Meaning |
| --- | --- |
| `404` | project not found |
| `404` | user not found |
| `404` | project membership not found |
| `409` | operation would remove the last owner |
