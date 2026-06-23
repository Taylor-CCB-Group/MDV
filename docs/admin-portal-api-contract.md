# Admin Portal API Contract

Status: Draft for team review

Read first: [admin-portal-mvp-plan.md](./admin-portal-mvp-plan.md)

Related architecture doc: [admin-portal-architecture.md](./admin-portal-architecture.md)

## Summary

All Admin Portal APIs live under:

```text
/admin/api/*
```

All responses are JSON.

The admin frontend should use these APIs instead of importing catalog, project-manager, or database internals.

## Open API Questions

- If an email exists in Auth0 but not in this deployment's configured Auth0 database connection, should Admin reuse that identity, create/add the user in the deployment connection, or reject the operation?
- Which Auth0 connection should Admin use for user creation?
- Are Auth0 Management API credentials available inside the MDV app container?
- When Admin creates a new Auth0 user, should the backend trigger an Auth0 invitation/password-reset email, rely on the normal forgot-password flow, or create the user silently?
- Should deleted or archived projects be returned by `GET /admin/api/projects`?

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

### Admin Project Access

```ts
type AdminProjectAccess = {
  projectId: number;
  permission: ProjectPermission;
  canRead: boolean;
  canWrite: boolean;
  isOwner: boolean;
};
```

### Admin Project Member

```ts
type AdminProjectMember = {
  user: AdminUser;
  projectAccess: AdminProjectAccess;
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

Local development note:

```text
Without Auth0, this returns local/dev MDV users.
For local POC and MVP implementation details, see [admin-portal-mvp-plan.md](./admin-portal-mvp-plan.md).
```

## POST `/admin/api/users`

Purpose:

```text
Create/add a user for this deployment with zero or more initial project permissions.
```

Implementation note:

```text
`projectAccess` is an array and may be empty. The backend policy for requiring at least one project assignment is centralized so it can be tightened later without changing the API shape.
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
    },
    {
      "projectId": 12,
      "permission": "edit"
    }
  ]
}
```

Validation:

- `email` is required.
- `projectAccess` is optional; omitted means an empty list.
- every `projectAccess[].projectId` must exist in this deployment.
- every `projectAccess[].permission` must be one of `view`, `edit`, or `owner`.

Production behavior:

```text
1. Create user in Auth0 for this deployment, or resolve existing Auth0 user.
2. Sync/create local MDV User row.
3. Create the initial UserProject permission rows.
4. Return user and project access.
```

Atomicity:

```text
This operation should be atomic from the API caller's point of view.
If Auth0 user creation succeeds but MDV sync or permission creation fails, the backend should roll back the Auth0-created user where possible.
If rollback fails, the response should make the partial failure explicit.
```

Local development behavior:

```text
1. Do not call Auth0.
2. Create or resolve local/dev MDV user.
3. Create the initial UserProject permission rows.
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
    "isAdmin": false
  },
  "projectAccess": [
    {
      "projectId": 10,
      "permission": "view",
      "canRead": true,
      "canWrite": false,
      "isOwner": false
    },
    {
      "projectId": 12,
      "permission": "edit",
      "canRead": true,
      "canWrite": true,
      "isOwner": false
    }
  ],
  "created": true
}
```

Errors:

| Status | Meaning |
| --- | --- |
| `400` | invalid email, empty project access, invalid project ID, or invalid permission |
| `409` | user already exists and conflict handling is required |
| `502` | Auth0 operation failed |
| `500` | MDV sync or permission write failed; rollback should be attempted if Auth0 user creation already happened |

Recommended behavior:

```text
Return the existing/synced user when possible and continue assigning project access.
```

Post-MVP note:

```text
Bulk paste/import of users should be handled as a separate API design, probably backed by a background job because Auth0 Management API calls are rate-limited.
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

## GET `/admin/api/projects/:projectId/users`

Purpose:

```text
Return users assigned to a project.
```

Response:

```json
{
  "members": [
    {
      "user": {
        "id": 1,
        "email": "owner@example.com",
        "firstName": "",
        "lastName": "",
        "isActive": true,
        "isAdmin": false
      },
      "projectAccess": {
        "projectId": 10,
        "permission": "owner",
        "canRead": true,
        "canWrite": true,
        "isOwner": true
      }
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
  "user": {
    "id": 2,
    "email": "viewer@example.com",
    "firstName": "",
    "lastName": "",
    "isActive": true,
    "isAdmin": false
  },
  "projectAccess": {
    "projectId": 10,
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
  "user": {
    "id": 2,
    "email": "viewer@example.com",
    "firstName": "",
    "lastName": "",
    "isActive": true,
    "isAdmin": false
  },
  "projectAccess": {
    "projectId": 10,
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

```text
204 No Content
```

Errors:

| Status | Meaning |
| --- | --- |
| `404` | project not found |
| `404` | user not found |
| `404` | project membership not found |
| `409` | operation would remove the last owner |
