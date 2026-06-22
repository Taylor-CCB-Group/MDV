# Admin Portal Architecture

Status: Draft for team review

## Summary

The Admin Portal is a separate MDV administration surface at `/admin`.

It is responsible for managing users and project access within one MDV deployment/container.

The target responsibility split is:

```text
Auth0 -> identity and login for a deployment
MDV   -> admin status, projects, and project-level permissions
```

Admin should be implemented as a trusted MDV plugin/extension. The end-state target is for the admin plugin to live outside the MDV core repo while still being enabled, disabled, and loaded by MDV through its extension/plugin mechanism.

Related API contract: [admin-portal-api-contract.md](./admin-portal-api-contract.md)

## Goals

- Provide a GUI for adding users to an MDV deployment.
- Require new users to receive access to at least one project during creation/addition.
- Allow admins to manage project-level permissions: `view`, `edit`, and `owner`.
- Allow admins to remove a user's access from a project.
- Keep admin separate from the catalog and project viewer.
- Structure admin as a trusted MDV plugin/extension that can move outside the MDV core repo.
- Support a local POC mode without Auth0.

## Current Understanding

Each MDV deployment/container has its own Auth0 user database and its own MDV database.

When a user is added through Auth0 today, that means the user is added to the Auth0 database for a specific deployment. After that, MDV syncs or creates a matching local user row and grants project access permission.

Project access is meaningful only when the user has permission for at least one project in that deployment.

## Route Ownership

Target route layout:

```text
/                  -> existing catalog app
/project/:id       -> existing project viewer
/admin             -> admin frontend app
/admin/api/*       -> admin backend APIs
```

Rules:

- `/admin` is not a catalog route.
- `/admin` is not mounted from `Dashboard.tsx`.
- `/admin` has its own frontend entrypoint.
- `/admin/api/*` is owned by the admin backend extension.
- The admin frontend should call only `/admin/api/*`.

## Frontend Boundary

The admin frontend should:

- load session status from `/admin/api/session`
- load users from `/admin/api/users`
- load projects from `/admin/api/projects`
- create/add users through `/admin/api/users`
- manage project permissions through `/admin/api/projects/:projectId/users`
- show backend validation errors clearly
- handle empty users/projects/project-members states

The admin frontend should not:

- import catalog internals
- mount inside the catalog dashboard
- assume Auth0 Management API access
- treat frontend permission checks as security

## Plugin/Extension Architecture

Admin should be a modular, trusted MDV plugin/extension.

The goal is modular deployment, external packaging, and a clear host/plugin boundary.

This means:

- admin can be enabled or disabled per deployment
- admin has its own frontend bundle
- admin has its own backend extension routes
- admin declares required permissions
- admin should ultimately live outside the core MDV repo
- admin is treated as a trusted plugin, not arbitrary untrusted third-party code
- the plugin manifest declares required capabilities
- the MDV host enforces route, permission, and security boundaries

Current/incremental shape:

```text
python/mdvtools/dbutils/admin_extension.py
src/admin/admin_index.tsx
src/admin/AdminApp.tsx
src/admin/api.ts
python/mdvtools/templates/admin_build.html
```

End-state shape:

```text
trusted admin plugin package outside the MDV core repo
  -> frontend bundle for /admin
  -> backend extension routes for /admin/api/*
  -> manifest/config declaring routes, permissions, and artifacts
```

Future manifest shape:

```json
{
  "id": "admin",
  "frontend": {
    "route": "/admin",
    "entrypoint": "src/admin/admin_index.tsx",
    "template": "admin_build.html",
    "artifact": "js/admin.js"
  },
  "backend": {
    "routes": ["/admin/api/*"]
  },
  "permissions": [
    "admin:access",
    "admin:users:manage",
    "admin:projects:manage"
  ]
}
```

The target for admin is:

```text
extension/plugin registration through MDV configuration or manifest metadata
developer hot reload during local development
```

## Host API Direction

Because the admin plugin should eventually live outside the MDV core repo, it should communicate through stable host APIs instead of importing arbitrary MDV internals.

Frontend host API should eventually provide:

- session/current-user information
- an authenticated API fetch helper
- route/navigation registration, if needed
- plugin config access

Backend host API should eventually provide:

- auth/session helpers
- user service
- project service
- project access service
- plugin config and manifest access

## MVP User Operations

For the MVP, admin supports:

- create/add user only when at least one initial project permission is provided
- change a user's permission for a project
- remove a user's access from a project

For the MVP, "remove user" means:

```text
remove the user's access from a project
```

In future it could mean:

```text
delete the user from Auth0
remove the user from the deployment
deactivate the user globally
```

Deployment-level removal and Auth0 deletion/deactivation can be designed later as separate operations.

## Permission Model

Project permissions use the existing MDV permission model:


| UI permission | can_read | can_write | is_owner |
| ------------- | -------- | --------- | -------- |
| `view`        | `true`   | `false`   | `false`  |
| `edit`        | `true`   | `true`    | `false`  |
| `owner`       | `true`   | `true`    | `true`   |


## User Creation Rule

The normal admin flow should not create a user with no project access.

Required rule:

```text
Creating/adding a user through the Admin Portal requires at least one initial project permission.
```

Reason:

```text
A user with no project permission has no useful access inside that MDV deployment.
```

Recommended flow:

```text
1. Admin enters user details.
2. Admin selects one or more projects.
3. Admin chooses permission for each selected project.
4. Backend creates the user in Auth0 for this deployment.
5. Backend syncs/creates the MDV User row.
6. Backend creates UserProject permission rows.
7. User can log in and see assigned projects.
```

## Auth0 And MDV Responsibility Split

Auth0 is responsible for identity and login.

MDV is responsible for admin authorization and project authorization.


| Concern                        | Source of truth |
| ------------------------------ | --------------- |
| User identity                  | Auth0           |
| Login/session                  | Auth0           |
| Deployment Auth0 user database | Auth0           |
| Local user record              | MDV             |
| Admin status                   | MDV             |
| Project list                   | MDV             |
| Project membership             | MDV             |
| Project permission             | MDV             |


Important rule:

```text
Auth0 proves who the user is.
MDV decides whether that user is an admin and which projects they can access.
```

This means admin permission should not depend only on Auth0. An authenticated user should also be marked as an admin in MDV for that deployment.

## Admin Authorization

Production admin access should require:

```text
1. User is authenticated through Auth0.
2. User exists or is synced in the MDV deployment database.
3. MDV marks that user as an admin.
```

Current/simple representation:

```python
session["user"]["is_admin"] == True
```

Future representation may become a richer role model:

```text
admin:access
admin:users:manage
admin:projects:manage
```

Backend checks are authoritative.

Frontend checks are only for UX.

Expected auth behavior:


| Case                    | Response        |
| ----------------------- | --------------- |
| No session              | `401`           |
| Logged in but not admin | `403`           |
| Admin user              | request allowed |


## First Admin Bootstrap

There is one bootstrap problem:

```text
Only admins can create users, but the first deployment still needs a first admin.
```

This is not a blocker for the admin MVP.

Acceptable MVP bootstrap approaches:

- setup script creates first admin
- manual DB seed for the first deployment
- CLI command promotes a synced user to admin
- deployment config declares first admin email

Current MVP position:

```text
The first admin may be created outside the Admin Portal through a script, CLI, config, or manual DB seed.
```

## Failure Handling

Creating a user touches more than one system:

```text
Auth0 user database
MDV User row
MDV UserProject permission rows
```

The admin API should be atomic from the caller's point of view.

If Auth0 user creation succeeds but MDV user sync or project permission creation fails, the backend should roll back the Auth0-created user where possible and return an error.

If rollback fails, the backend should return an error that makes the partial failure explicit so it can be repaired manually or retried.

## Last Owner Rule

For the MVP, an admin should not be able to remove or demote the last owner of a project.

This avoids leaving a project with no owner.

If a future global-admin override is needed, it should be designed as a separate flow.

## POC Mode

The POC should be small and local, and should validate the MDV-side project permission workflow without requiring Auth0.

POC summary:

- Auth0 is not called.
- Authentication may be disabled.
- Backend may use a synthetic local admin.
- Users come from local/dev MDV data, with optional backend-seeded demo users.
- Project permission writes use MDV-style `UserProject` rows.
- API shapes should stay close to the final contract.

Detailed POC implementation plan: [admin-portal-poc-plan.md](./admin-portal-poc-plan.md)

## Open Decisions

- Should `User.is_admin` be enough, or do we need a role table?
- Should deleted/archived projects appear in admin project lists?
- What exact repair path should be used if Auth0 rollback fails after MDV sync/permission creation fails?
- Should permission changes be audited?
- Should external plugin registration happen at build time or startup time?

## Review Questions

The team should review and confirm:

- Is the Auth0/MDV responsibility split correct?
- Is admin status correctly treated as an MDV authorization concept?
- Is it correct that adding a user requires at least one project permission?
- Is it correct that MVP "remove user" means remove project access only?
- Is `/admin` as a separate extension surface the right direction?
- Is build-time/startup-time plugin registration enough for admin?

