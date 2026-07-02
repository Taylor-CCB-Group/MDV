# Auth0 local development setup

This guide walks through enabling real Auth0 login when running MDV locally with the dev container / `docker-secrets.yml` stack. It captures common pitfalls (missing API authorization, users not in the local database, Management API timeouts) so you do not have to rediscover them.

**Stack assumption:** you are using the dev container or `docker compose -f docker-secrets.yml up`, with the app on `http://localhost:5055`.

**Run Python commands inside the dev container** unless stated otherwise. The container already has dependencies installed and reads env vars from `docker-secrets.yml`. Running scripts on the host often fails because `DB_HOST=mdv_db` and Auth0 env vars are not set there.

---



## Overview


| Layer                              | What it does                                              |
| ---------------------------------- | --------------------------------------------------------- |
| **Auth0**                          | Identity provider — login, users, roles                   |
| `docker-secrets.yml`               | Local compose stack; sets runtime env vars                |
| **MDV Postgres** `users` **table** | Local copy of users MDV knows about                       |
| `user_projects` **table**          | Project permissions (`is_owner`, `can_read`, `can_write`) |


Auth0 proves who someone is. MDV still requires a matching row in the local `users` table before login succeeds. Users are **not** synced automatically on startup — you must run a sync step once (see [Sync Auth0 users into the database](#sync-auth0-users-into-the-database)).

For testing the Admin portal UI without Auth0, use dummy auth instead (`ENABLE_AUTH=0`, `DEFAULT_AUTH_METHOD=dummy`). See [Alternative: dummy auth (no Auth0)](#alternative-dummy-auth-no-auth0).

---



## Part 1 — Auth0 dashboard setup



### 1. Create an Auth0 account and tenant

Sign up at [auth0.com](https://auth0.com) (free tier is enough for local dev). Note your **tenant domain**, e.g. `your-tenant.uk.auth0.com`.

### 2. Create an application

1. **Applications → Create application**
2. Name: e.g. `MDV Local`
3. Type: **Regular Web Application**

Under **Settings → Application URIs**, set:


| Field                  | Value                                                       |
| ---------------------- | ----------------------------------------------------------- |
| Allowed Callback URLs  | `http://localhost:5055/callback`                            |
| Allowed Logout URLs    | `http://localhost:5055/login` (or `http://localhost:5055/`) |
| Allowed Web Origins    | `http://localhost:5055`                                     |
| Allowed Origins (CORS) | `http://localhost:5055`                                     |


Save changes. Copy **Client ID** and **Client Secret** — you will need them for `docker-secrets.yml`.

### 3. Create an API

1. **Applications → APIs → Create API**
2. Name: e.g. `MDV Local API`
3. Identifier: e.g. `https://mdv-local-api` — this must match `AUTH0_AUDIENCE` exactly
4. Signing algorithm: **RS256**



### 4. Authorize the application on your API (easy to miss)

Login will fail or loop if this step is skipped.

1. **Applications → APIs → MDV Local API** (your API)
2. Open the **Application Access** tab
3. Find your **MDV Local** application and **authorize** it
4. Grant both **user access** and **client access** (the UI may show `0` permissions with a green tick — that is normal for the custom API)

Without this, the OAuth flow can complete at Auth0 but MDV token validation or consent behavior may be wrong.

### 5. Authorize the application on Auth0 Management API

User sync uses the **Auth0 Management API** (list users, read roles). Your application must be allowed to call it.

1. **Applications → APIs → Auth0 Management API**
2. **Machine To Machine Applications** (or **Application Access**)
3. Select your **MDV Local** application
4. Enable at least:


| Scope          | Used for                                     |
| -------------- | -------------------------------------------- |
| `read:users`   | Sync users into MDV                          |
| `read:roles`   | Set `User.is_admin` from Auth0 role          |
| `create:users` | Admin portal user creation (if testing that) |
| `delete:users` | Admin portal user deletion (if testing that) |


Authorize and save.

### 6. Create an `admin` role

MDV sets `User.is_admin = true` when the Auth0 role name is exactly `admin` (lowercase).

1. **User Management → Roles → Create role**
2. Name: `admin`
3. Save (no special permissions required on the role for local MDV admin access)



### 7. Create a test user

1. **User Management → Users → Create user**
2. Use **Username-Password-Authentication** (the default database connection)
3. Set email and password
4. **Roles** tab → assign the `admin` role

Note the user’s **User ID** (e.g. `auth0|6a44f83e...`) — useful if you need manual DB insert later.

The connection name `Username-Password-Authentication` must match `AUTH0_DB_CONNECTION` in `docker-secrets.yml` exactly (see [Environment variables](#part-2--mdv-environment-variables)).

---



## Part 2 — MDV environment variables

Edit `docker-secrets.yml` under `mdv_app` → `environment:`.

Use placeholders in docs and git; never commit real client secrets.

```yaml
# Authentication
- ENABLE_AUTH=1
- DEFAULT_AUTH_METHOD=auth0
- FLASK_SECRET_KEY=local-dev-secret-change-me
- LOGIN_REDIRECT_URL=/login

# Auth0
- AUTH0_DOMAIN=your-tenant.uk.auth0.com
- AUTH0_CLIENT_ID=your-client-id
- AUTH0_CLIENT_SECRET=your-client-secret
- AUTH0_CALLBACK_URL=http://localhost:5055/callback
- AUTH0_AUDIENCE=https://mdv-local-api
- AUTH0_PUBLIC_KEY_URI=https://your-tenant.uk.auth0.com/.well-known/jwks.json
- AUTH0_DB_CONNECTION=Username-Password-Authentication
```


| Variable               | Notes                                                                         |
| ---------------------- | ----------------------------------------------------------------------------- |
| `ENABLE_AUTH`          | `1`, `true`, or `yes`                                                         |
| `DEFAULT_AUTH_METHOD`  | Must be `auth0` for real login                                                |
| `FLASK_SECRET_KEY`     | Required for sessions when auth is on                                         |
| `LOGIN_REDIRECT_URL`   | Where unauthenticated users are sent; `/login` is fine                        |
| `AUTH0_DOMAIN`         | Tenant domain only, no `https://`                                             |
| `AUTH0_AUDIENCE`       | Must match API **Identifier** exactly                                         |
| `AUTH0_PUBLIC_KEY_URI` | `https://<domain>/.well-known/jwks.json`                                      |
| `AUTH0_DB_CONNECTION`  | Auth0 database connection name; default is `Username-Password-Authentication` |


`AUTH0_DB_CONNECTION` is required for user sync. Sync only imports users from that connection.

### Apply config changes

Restart the stack so env vars are picked up:

```bash
docker compose -f docker-secrets.yml up -d --force-recreate mdv_app
```

Or restart the dev container from your IDE.

---



## Part 3 — Sync Auth0 users into the database

**Do this before the first login attempt.**

MDV’s login callback stores an Auth0 token, then on the next request looks up the user by `auth_id` in Postgres. If the user is missing, you get a redirect loop: home → `/login` → Auth0 consent → `/callback` → home → …

Sync is **not** run on app startup (by design — Management API rate limits). Run it manually:

### Sync command (inside dev container)

```bash
cd /app/python && uv run python -c "
from mdvtools.dbutils.safe_mdv_app import app
from mdvtools.auth.authutils import get_auth_provider, cache_user_projects

with app.app_context():
    get_auth_provider().sync_users_to_db()
    cache_user_projects()
    print('sync complete')
"
```



### Verify users in Postgres

```bash
docker compose -f docker-secrets.yml exec mdv_db \
  psql -U my_db_user -d my_db_name -c \
  "SELECT id, email, auth_id, is_admin FROM users;"
```

You should see your Auth0 user with `auth_id` like `auth0|...` and `is_admin = t` if the `admin` role was assigned.

### What sync does

- Pulls users from Auth0 (`AUTH0_DB_CONNECTION`)
- Upserts into `users` via `UserService`
- Sets `is_admin` from Auth0 role named `admin`
- Grants **all existing projects** to admin users as owner in `user_projects`

Alternative entry point (also runs sync first):

```bash
cd /app/python && uv run python mdvtools/scripts/manage_project_permissions.py assign \
  --email admin@test.com \
  --project some_project_name \
  --permission owner
```

---



## Part 4 — Log in and test

1. Open `http://localhost:5055/` (or `/login`)
2. Complete Auth0 login
3. Confirm session: `GET http://localhost:5055/admin/api/session` → `"isAdmin": true` for admin users
4. Open Admin portal: `http://localhost:5055/admin`

---



## Troubleshooting



### Redirect loop after Auth0 consent

**Symptom:** Authorize App screen → callback → home → `/login` → Auth0 again.

**Most common cause:** Auth0 user not in local DB. Run [sync](#sync-command-inside-dev-container) and verify `users` table.

**Check app logs** for `User not found` or `Unauthorized access ... Redirecting`.

### Sync fails with `Read timed out (read timeout=5.0)`

**Symptom:** `HTTPSConnectionPool(host='....auth0.com'...): Read timed out`

**Cause:** Docker → Auth0 HTTPS can take **>5 seconds**; `auth0-python` used to default to a 5s timeout.

**Fix:** The codebase sets `AUTH0_API_TIMEOUT = 30.0` in `auth0_provider.py`. Ensure you are on a revision that includes that change, then re-run sync.

**Diagnose network latency** (inside dev container):

```bash
time curl -s -o /dev/null -w "%{http_code}\n" \
  "https://your-tenant.uk.auth0.com/.well-known/openid-configuration"
```

If this takes longer than 5s, sync needs the increased timeout.

**Verify Management API credentials** (inside dev container):

```bash
cd /app/python && uv run python -c "
import os, requests
domain = os.environ['AUTH0_DOMAIN']
resp = requests.post(
    f'https://{domain}/oauth/token',
    json={
        'client_id': os.environ['AUTH0_CLIENT_ID'],
        'client_secret': os.environ['AUTH0_CLIENT_SECRET'],
        'audience': f'https://{domain}/api/v2/',
        'grant_type': 'client_credentials',
    },
    timeout=30,
)
print(resp.status_code, resp.text[:200])
"
```

Expect `200` and an `access_token`. `403` / `unauthorized_client` → fix Management API application access (Part 1 step 5).

### Manual user insert (if sync is blocked)

1. Auth0 Dashboard → Users → copy **User ID** (`auth0|...`)
2. Insert into Postgres:

```sql
INSERT INTO users (
  email, auth_id, is_admin,
  password, is_active, first_name, last_name, administrator
) VALUES (
  'admin@test.com',
  'auth0|PASTE_USER_ID_HERE',
  true,
  '', false, '', '', false
);
```

1. Refresh cache (inside dev container):

```bash
cd /app/python && uv run python -c "
from mdvtools.dbutils.safe_mdv_app import app
from mdvtools.auth.authutils import cache_user_projects
with app.app_context():
    cache_user_projects()
"
```



### Session cookie not persisted (`http://localhost`)

MDV sets `SESSION_COOKIE_SECURE=True`. Some browsers will not store secure cookies over plain HTTP. After `/callback`, check DevTools → Application → Cookies for `localhost:5055`. If empty, login may fail even with a DB user. Workarounds: use HTTPS locally, or add a dev-only env override for secure cookies (not yet in stock config).

### Stale user list in running server

The running Gunicorn process caches permissions in memory. After DB changes, restart:

```bash
docker compose -f docker-secrets.yml restart mdv_app
```



### `projects.owner` is NULL

The `projects.owner` column is **legacy and unused**. Ownership is stored in `user_projects.is_owner`. Empty `projects.owner` is normal and does not break the app.

### Cleaning up local dummy users

If you created test users with `local:...` auth IDs before enabling Auth0, delete them from `user_projects` then `users`, and refresh the cache. Keep your Auth0-synced admin user (e.g. `auth0|...`).

---



## Alternative: dummy auth (no Auth0)

For Admin portal or project-permission UI work without Auth0:

```yaml
- ENABLE_AUTH=0
- DEFAULT_AUTH_METHOD=dummy
```

The Admin extension treats you as a synthetic local admin when auth is off. Switch back to Auth0 when testing real login.

---



## Quick reference — command cheat sheet

All Python commands: **run inside the dev container**.


| Task                         | Command                                                                                                                                   |
| ---------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| Start stack                  | `docker compose -f docker-secrets.yml up -d`                                                                                              |
| Restart app after env change | `docker compose -f docker-secrets.yml restart mdv_app`                                                                                    |
| Sync Auth0 → DB              | See [sync command](#sync-command-inside-dev-container)                                                                                    |
| List users in DB             | `docker compose -f docker-secrets.yml exec mdv_db psql -U my_db_user -d my_db_name -c "SELECT id, email, auth_id, is_admin FROM users;"`  |
| Assign project permission    | `cd /app/python && uv run python mdvtools/scripts/manage_project_permissions.py assign --email USER --project PROJECT --permission owner` |


---



## Related docs

- [DATABASE_BACKEND_MANUAL.md](./DATABASE_BACKEND_MANUAL.md) — auth env vars, permission model, API routes
- [admin-portal-mvp-plan.md](./admin-portal-mvp-plan.md) — Admin portal and Auth0 integration plan
- `docker-secrets.yml` — local dev compose and env template
- `python/mdvtools/auth/auth0_provider.py` — login, validation, sync implementations

