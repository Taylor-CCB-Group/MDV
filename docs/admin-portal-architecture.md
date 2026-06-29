# Admin Portal Architecture

Status: Draft for team review

Read first: [admin-portal-mvp-plan.md](./admin-portal-mvp-plan.md)

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
- Allow new users to be created with zero or more project assignments.
- Allow admins to manage project-level permissions: `view`, `edit`, and `owner`.
- Allow admins to remove a user's access from a project.
- Keep admin separate from the catalog and project viewer.
- Structure admin as a trusted MDV plugin/extension that can move outside the MDV core repo.
- Support a local POC mode without Auth0.

## Current Understanding

Each MDV deployment/container has its own Auth0 user database and its own MDV database.

When a user is added through Auth0 today, that means the user is added to the Auth0 database for a specific deployment. After that, MDV syncs or creates a matching local user row and grants project access permission.

Project access is meaningful only when the user has permission for at least one project in that deployment, but some deployment users may temporarily exist without project access.

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

This should remain a strict boundary. The Admin frontend should not call catalog or project-manager APIs directly; the Admin backend should expose admin-shaped APIs backed by MDV host services.

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

The target registration model is startup-time plugin registration with build-time/frontend artifact materialisation for production images. The admin plugin package is installed before MDV starts, MDV reads plugin config/manifest at startup, registers `/admin` and `/admin/api/*`, and serves frontend assets that were already built/copied into a Flask-visible location.

## Host API Direction

Because the admin plugin should eventually live outside the MDV core repo, it should communicate through stable host APIs instead of importing arbitrary MDV internals.

Admin should not own MDV domain data for MVP. The plugin owns its frontend code, backend route code, and manifest/config needed to attach to MDV. Users, projects, project permissions, auth/session state, and admin markers remain MDV-owned and should be accessed through MDV-provided host APIs/services.

Frontend host API should eventually provide:

- session/current-user information
- an authenticated API fetch helper
- route/navigation registration, if needed
- plugin config access

Backend host API should eventually provide:

- auth/session helpers
- identity provider/user create-or-resolve service
- user service
- project service
- project access service
- plugin config and manifest access

The current service boundary should remain admin-specific until Auth0-backed user creation is proven. The current implementation already splits Admin database operations from an `AdminIdentityProvider` boundary, so local dev and Auth0-backed identity behavior can differ without changing Admin route handlers. After the Auth0 baseline exists, reusable host-service concepts can be extracted into a generic plugin host API with less guesswork.

Admin should reuse existing MDV services where possible, especially user and project-access services. If current service transaction boundaries prevent the Admin create-user flow from behaving atomically, refactor those services or add transaction-aware variants rather than duplicating business logic in Admin.

## MVP User Operations

For the production MVP, admin supports:

- create/add user with zero or more initial project permissions
- list deployment users
- list existing projects
- change a user's permission for a project
- remove a user's access from a project

The user list should show all deployment users, not only users with project access. The UI should distinguish users with no project access, for example with a project-access count or status indicator, so admins can repair incomplete/manual states.

The first production flow can remain project-centric: select a project, then manage users and permissions for that project. A user-centric access view, where an admin selects a user and sees all projects they can access, should be a follow-up after Auth0-backed user creation.

Project lifecycle management is outside MVP. Admin should not create, delete, archive, restore, or edit projects in the first production increment.

For production MVP, Admin write actions should be logged server-side with the actor, action, target user, target project, old/new permission where relevant, timestamp, and success/failure. A dedicated audit-log table/API/UI can come later.

For the MVP, "remove user" means:

```text
remove the user's access from a project
```

The UI should label this action as "Remove access", not "Remove user" or "Delete user".

In future it could mean:

```text
delete the user from Auth0
remove the user from the deployment
deactivate the user globally
```

Deployment-level removal and Auth0 deletion/deactivation can be designed later as separate operations.

Managing deployment-level admin status is also outside MVP. Admin should not promote or demote other admins in the first production increment; additional admins should come from Auth0 role sync or another approved bootstrap process.

## Post-MVP Bulk User Import

A requested post-MVP feature is a quick paste/import flow where an admin can paste a list of users and create accounts/project access in bulk.

Because existing Auth0 sync code already warns about Management API rate limits, bulk import should probably not be implemented as one long synchronous request. It should be designed as a background job or queued operation with progress, partial-failure reporting, and retry behavior.

## Permission Model

Project permissions use the existing MDV permission model:


| UI permission | can_read | can_write | is_owner |
| ------------- | -------- | --------- | -------- |
| `view`        | `true`   | `false`   | `false`  |
| `edit`        | `true`   | `true`    | `false`  |
| `owner`       | `true`   | `true`    | `true`   |

MVP admins can assign `owner` permission. Last-owner protection is still required so a project cannot be left without an owner.

Admin permission changes should take effect immediately for users. In production/auth mode, Admin writes must refresh the user-project cache after changing project access so the user experience does not depend on a manual script or server restart.

## User Creation Rule

The admin create-user flow currently allows users with no project access.

Current rule:

```text
Creating/adding a user through the Admin Portal accepts zero or more initial project permissions.
```

The backend should keep this policy centralized so it can later require at least one initial project permission without changing the API shape.

Reason:

```text
Some deployment users may need to exist before project access is assigned. A user with no project permission should be visible in Admin but cannot access projects until permissions are added.
```

Recommended flow:

```text
1. Admin enters user details.
2. Admin optionally selects one or more projects.
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

For MVP, Admin should use the existing Auth0 application credentials for Management API calls, matching current `sync_users_to_db()` behavior. A later deployment can add optional dedicated M2M credentials if stronger separation is required.

Admin should create or resolve Auth0 users only in the current deployment's configured `AUTH0_DB_CONNECTION`.

If an email exists elsewhere in Auth0 but not in the configured `AUTH0_DB_CONNECTION`, Admin should create/add the user in `AUTH0_DB_CONNECTION` rather than reusing an identity from another connection.

For MVP, newly-created Auth0 users should receive a generated random password. Admin should not expose that password and should not add new invitation/password-reset email automation. Users can use the existing password reset/forgot-password path if needed; explicit setup-email automation can be added later.

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

This is enough for MVP Admin access.

Existing MDV Auth0 sync code sets `User.is_admin` from an Auth0 role named `admin`. Some deployment notes mention user metadata such as `{"role": "admin"}`, but MVP Admin authorization should follow the current code path and should not treat metadata as authoritative.

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

Preferred MVP position:

```text
The first admin should be bootstrapped through existing Auth0 role sync where possible: assign Auth0 role "admin", sync users into MDV, and rely on User.is_admin=true.
```

Manual DB seed, CLI promotion, or deployment config can remain fallback options where Auth0 role sync is not available.

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

## Local Development Mode

Local development should validate the MDV-side project permission workflow without requiring Auth0.

Local development summary:

- Auth0 is not called.
- Authentication may be disabled.
- Backend may use a synthetic local admin.
- Users come from local/dev MDV data, with optional backend-seeded demo users.
- Project permission writes use MDV-style `UserProject` rows.
- API shapes should stay close to the final contract.

Detailed MVP implementation plan: [admin-portal-mvp-plan.md](./admin-portal-mvp-plan.md)

## Remaining Architecture Questions

The read-first source for implementation status and next steps is [admin-portal-mvp-plan.md](./admin-portal-mvp-plan.md).

Still open:

- Should post-MVP add explicit Auth0 invitation/password-reset email automation?

Settled MVP decisions:

- `User.is_admin` is enough for MVP Admin access.
- Admin status is an MDV authorization concept, while Auth0 proves identity/login.
- Adding or creating a user through Admin accepts zero or more project permissions; the policy can be tightened later without changing the API shape.
- MVP "remove user" means remove project access only.
- `/admin` is a separate extension surface, not a catalog page.
- Admin is a trusted plugin/extension.
- Plugin registration target is startup-time, with build-time frontend asset materialisation for production.
- Auth0 Management API calls use existing Auth0 application credentials for MVP, with optional dedicated M2M credentials later.
- Auth0 user creation uses the current deployment's `AUTH0_DB_CONNECTION`.
- Existing Auth0 emails are resolved only inside `AUTH0_DB_CONNECTION`; if not present there, Admin creates/adds the user in that connection.
- MVP Auth0 user creation uses a generated password, does not expose it in Admin, and does not add new email automation.
- Admin project lists exclude deleted projects by default; deleted/archived visibility can be added later through an explicit filter.
- Admin writes should refresh `user_project_cache` when `ENABLE_AUTH=true`.
- Last-owner conflicts return `409 Conflict`.
- Production MVP should start with server-side Admin write logging, not a full audit-log UI/table.
