# Admin Portal MVP Plan

Status: Draft for implementation review

Read this first if you are implementing or reviewing Admin Portal work. This document explains what is already implemented locally, what the MVP should include, what is intentionally out of scope, and what the next implementation slices are.

Related docs:

- [admin-portal-architecture.md](./admin-portal-architecture.md) - target architecture, plugin model, Auth0/MDV responsibility split
- [admin-portal-api-contract.md](./admin-portal-api-contract.md) - HTTP API request/response shapes and error behavior

## 1. Current Status

We have a local Admin Portal POC that proves the MDV-side project-access workflow without Auth0.

```text
/admin
  -> standalone Admin frontend

/admin/api/*
  -> AdminExtension HTTP endpoints

AdminExtension
  -> admin contracts
  -> MDVAdminServices
  -> MDV DB models: User, Project, UserProject
```

Implemented locally:

- list projects
- list users
- create/reuse local user
- require initial project access
- view project members
- add existing user to project
- change permission: `view`, `edit`, `owner`
- remove project access
- protect last project owner with `409 Conflict`

This confirms the local MDV database workflow. It does not confirm Auth0 user creation, production login, catalog filtering for newly-created users, or production Admin authorization.

## 2. MVP Scope

The production MVP should let an Admin manage user access to existing projects.

In scope:

- list deployment users
- list existing projects
- create or resolve a user for this deployment
- allow zero or more initial project assignments when creating/adding a user
- assign project permissions during user creation when provided
- view project members
- add an existing deployment user to a project
- change a user's project permission: `view`, `edit`, `owner`
- remove a user's access from a project
- protect the last project owner
- refresh project-access cache after writes when `ENABLE_AUTH=true`
- log Admin write actions server-side

Out of scope for MVP:

- project create/delete/archive/restore/lifecycle management
- admin promotion/demotion workflow
- deployment-level remove/deactivate/delete user workflow
- user-centric "projects for this user" view
- audit-log UI/table
- bulk user paste/import
- plugin marketplace UI
- runtime upload/loading of arbitrary plugin code
- untrusted third-party plugin sandboxing

## 3. Decision Status

| Decision | Status |
| --- | --- |
| Admin is a separate app surface at `/admin` | Decided |
| Admin frontend calls only `/admin/api/*` | Decided |
| Admin is a trusted plugin/extension | Decided |
| Untrusted plugin sandboxing | Later |
| Target plugin registration model | Startup-time registration |
| Production frontend artifact handling | Build-time materialisation |
| Admin owns MDV domain data | No |
| Admin uses MDV APIs/services for users/projects/permissions | Decided |
| Service boundary remains admin-specific until Auth0 baseline | Decided |
| MVP Admin authorization marker | `User.is_admin` |
| First admin bootstrap | Auth0 role `admin` sync where possible |
| Admin can assign `owner` | Yes |
| Last-owner violation status | `409 Conflict` |
| Create-user flow requires project access | No, but the backend rule is centralized so this can be tightened later |
| Implemented create-user `projectAccess` shape | Array, may be empty |
| Auth0 create/resolve user | Decided |
| Auth0 invite/password setup flow | Open |
| Auth0 existing email outside deployment connection behavior | Open |
| Auth0 connection name/source | Open |
| Auth0 Management API credentials in app container | Open |
| Deleted/archived project visibility | Open |
| Bulk paste/import | Post-MVP |

## 4. Increment Ledger

```text
Increment 1: Admin read-only shell                 Status: implemented baseline
Increment 2: Local user creation and assignment    Status: implemented baseline
Increment 3A: Project members view                 Status: implemented baseline
Increment 3B: Permission editing/removal           Status: implemented baseline
Increment 3C: Create-user API array shape          Status: implemented baseline
Increment 4: Auth0-backed user creation            Status: planned
```

Implemented Admin endpoints:

```text
GET    /admin
GET    /admin/api/session
GET    /admin/api/users
POST   /admin/api/users
GET    /admin/api/projects
GET    /admin/api/projects/:projectId/users
POST   /admin/api/projects/:projectId/users
PATCH  /admin/api/projects/:projectId/users/:userId
DELETE /admin/api/projects/:projectId/users/:userId
```

Current local verification loop:

```text
open /admin
  -> select project
  -> create local user with initial access
  -> select project in Project Members
  -> see assigned user and permission
  -> change permission
  -> remove project access
```

## 5. Next Implementation Slices

### Increment 3C: Align Create-User API Shape

Status: implemented baseline.

Goal: align the local POC API with the production MVP shape before adding Auth0.

Scope:

- create-user request `projectAccess` is an array
- create-user response `projectAccess` is an array
- zero or more initial project assignments are allowed
- frontend Zod schema and UI use the array shape
- backend parser, service, and tests use the array shape
- local/no-auth behavior remains supported

### Increment 4: Auth0-Backed User Creation

Goal: replace local-only user creation with production Auth0 create/resolve behavior.

Expected flow:

```text
Admin submits email + zero or more project assignments
  -> backend creates or resolves Auth0 user in the deployment
  -> backend syncs/creates local MDV User
  -> backend creates UserProject permission rows
  -> backend refreshes project-access cache when ENABLE_AUTH=true
  -> backend returns success only after the full workflow succeeds
```

User creation should be atomic from the API caller's point of view. If MDV user/project-access creation fails after creating a new Auth0 user, Admin should attempt to roll back the Auth0 user. If rollback fails, return a clear partial-failure error for manual repair.

### Later Slices

- verify production auth/session behavior with `ENABLE_AUTH=true`
- add server-side structured logging for Admin writes
- add user access count/status in the users list
- add user-centric "projects for this user" view
- design bulk paste/import as a background job
- extract reusable host-service concepts after Auth0 baseline is proven
- design generic MDV plugin manifest/host API
- move Admin from in-repo extension to trusted external plugin package

## 6. Plugin Direction

Admin is being designed as a trusted MDV plugin/extension.

Current staged path:

```text
admin-specific seam inside MDV
  -> admin MVP workflows
  -> Auth0-backed baseline
  -> identify reusable plugin concepts
  -> general MDV plugin interface
  -> external trusted admin plugin package
```

Important boundaries:

- Admin frontend should not import catalog/project-manager frontend internals.
- Admin frontend should not call catalog/project-manager APIs directly.
- Admin frontend should use `/admin/api/*` only.
- MDV model imports should stay isolated in the Admin service layer.
- Admin should reuse existing MDV user/project-access services where possible.
- If existing service transaction boundaries block atomic Admin workflows, refactor those services or add transaction-aware variants rather than duplicating business logic.
- Admin owns plugin code, routes, and frontend assets; MDV owns users, projects, permissions, auth/session state, and admin markers.

The target plugin loading model is startup-time registration, not runtime upload/loading. Production still needs build-time materialisation of frontend assets so Flask can serve the Admin bundle. Local development can use Vite hot reload/proxy.

## 7. Auth0 Decisions

Decided:

- Admin should create or resolve users immediately because avoiding manual Auth0 work is the core pain point.
- For MVP production authorization, `User.is_admin` is the Admin access marker.
- First admin bootstrap should use existing Auth0 role sync where possible: assign Auth0 role `admin`, sync users into MDV, and rely on `User.is_admin=true`.
- Admin writes should refresh `user_project_cache` immediately when `ENABLE_AUTH=true`.

Still open:

- If an email exists in Auth0 but not in this deployment's configured Auth0 database connection, should Admin reuse that identity, create/add the user in the deployment connection, or reject the operation?
- Which Auth0 connection should Admin create users in?
- Are Auth0 Management API credentials available inside the MDV app container?
- If Admin creates a new Auth0 user, what exact account-setup flow should be used: Auth0 invitation, Auth0 password-reset ticket/email, normal forgot-password flow, or no email?

## 8. Module Map

- `src/admin/AdminApp.tsx`
  - Admin frontend module.
  - Only talks to `/admin/api/*`.
- `src/admin/api.ts`
  - Admin frontend HTTP module.
  - Validates Admin API responses with Zod schemas before returning typed data.
- `src/admin/schemas.ts`
  - Admin frontend schema module.
  - Defines Zod schemas and inferred TypeScript types for Admin API request/response shapes.
- `python/mdvtools/dbutils/admin_extension.py`
  - Admin route module.
  - Handles Flask routes, request parsing, admin session checks, and HTTP error mapping.
- `python/mdvtools/dbutils/admin_contracts.py`
  - Admin domain interface module.
  - Defines user, project, project-access, and project-member shapes.
- `python/mdvtools/dbutils/admin_services.py`
  - MDV adapter at the Admin seam.
  - Uses MDV `User`, `Project`, and `UserProject`.
- `python/mdvtools/auth/auth0_provider.py`
  - Auth0 login/sync module.
  - Existing sync pulls Auth0 users into MDV and updates `User.is_admin`.
- `python/mdvtools/auth/auth_invitation`
  - Script-like Auth0 user upsert module.
  - Not yet shaped for Admin.
- `python/mdvtools/dbutils/project_manager_extension.py`
  - Existing catalog/project sharing module.
  - Auth-mode oriented and exits early when `ENABLE_AUTH=false`.
