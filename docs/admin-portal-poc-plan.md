# Admin Portal POC Plan

Status: Draft for implementation review

Related docs:

- [admin-portal-architecture.md](./admin-portal-architecture.md)
- [admin-portal-api-contract.md](./admin-portal-api-contract.md)

## 1. POC Goal

Build the Admin Portal in small increments so we can prove the architecture before implementing the full user-management workflow.

The POC does not need real Auth0 integration.

## 2. Increment Overview

```text
Increment 1: Admin read-only shell                 Status: implemented baseline
Increment 2: Local user creation and assignment    Status: next
Increment 3: Permission editing/removal            Status: later
Increment 4: Auth0-backed user creation            Status: later
```

Each increment should update this plan before implementation if the scope changes.

## 3. POC Architecture

```text
/admin
  -> standalone React admin frontend

/admin/api/*
  -> Flask admin APIs registered by AdminExtension

MDV DB
  -> users, projects, user_project permissions
```

Important boundaries:

- Admin is not mounted in catalog.
- Admin is not a `Dashboard.tsx` route.
- Admin frontend calls only `/admin/api/*`.
- Backend performs all real permission checks.
- Frontend permission checks are UX only.

## 4. Auth Model For POC

For local/demo mode:

```text
ENABLE_AUTH=false
```

Backend may return a synthetic local admin:

```json
{
  "id": 0,
  "email": "local-admin@localhost",
  "is_admin": true,
  "synthetic": true
}
```

POC auth limitations:

- no Auth0 login
- no Auth0 Management API calls
- no production admin bootstrap
- no production session hardening beyond route-level checks/tests

## 5. Increment 1: Admin Read-Only Shell

Status: implemented baseline.

Increment 1 proves the admin plugin/extension shape with the smallest useful vertical slice.

Implemented scope:

- `GET /admin`
- `GET /admin/api/session`
- `GET /admin/api/users`
- `GET /admin/api/projects`
- standalone `src/admin` frontend entrypoint
- Flask template for the built admin shell
- Vite admin build entry
- Vite dev rewrite for `/admin`
- Vite dev proxy for `/admin/api`
- focused backend tests for session, users, projects, and empty states

Intentionally not included:

- user creation
- backend demo-user seeding
- project assignment
- permission editing
- project access removal
- Auth0 integration

Current Increment 1 success criteria:

- `/admin` loads as a separate app
- `/admin` is not a catalog route
- `/admin/api/session` works in local mode
- projects are listed when they exist
- users are listed when they exist
- empty user/project states are clear
- backend tests pass
- Vite build emits the admin bundle

## 6. Increment 2: Local User Creation And Assignment

Status: next.

Update this section before implementation.

Goal:

```text
Create a local MDV user without Auth0 and assign that user to at least one project.
```

Likely scope:

- add `POST /admin/api/users`
- require at least one `projectAccess` entry
- create local MDV `User` row
- create initial `UserProject` permission rows
- add frontend form for local user creation with project assignment
- reject user creation with no project access

Questions to resolve before implementation:

- should local user creation allow creating an already-existing email, or should it return the existing user?
- should Increment 2 include backend demo-user seeding, or is the local create-user form enough for demos?
- should local user creation use the same response shape as future Auth0-backed creation?

## 7. Increment 3: Permission Editing And Removal

Status: later.

Update this section before implementation.

Goal:

```text
Manage existing project access after a user exists.
```

Likely scope:

- add `GET /admin/api/projects/:projectId/users`
- add `POST /admin/api/projects/:projectId/users`
- add `PATCH /admin/api/projects/:projectId/users/:userId`
- add `DELETE /admin/api/projects/:projectId/users/:userId`
- show project members
- add existing user to another project
- change permission between `view`, `edit`, and `owner`
- remove project access
- enforce last-owner protection

## 8. Increment 4: Auth0-Backed User Creation

Status: later.

Update this section before implementation.

Goal:

```text
Replace local-only user creation with Auth0-backed user creation for the deployment.
```

Likely scope:

- create or resolve user in Auth0 deployment database
- sync/create local MDV `User` row
- create initial `UserProject` permission rows
- roll back Auth0 user if MDV sync or permission write fails
- expose clear partial-failure errors if rollback fails
- keep frontend API shape stable from Increment 2

## 9. Expected POC Limitations

- no Auth0 user creation until Increment 4
- no first-admin bootstrap
- no audit logging
- no external plugin discovery
- admin still lives inside MDV repo during the POC
