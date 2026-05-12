# Playwright Setup And Run

This is the first document to read when you want to install dependencies, set
up the environment, and run Playwright tests.

## Default Environment

Prefer running Playwright against the app on `http://localhost:5055`.

That is the default path for:

- backend-backed `tests_playwright/project/*`
- the default supported mixed-suite run
- most local verification

Use a local Vite dev server on `http://127.0.0.1:5173` only when you
intentionally want the dev-only catalog path.

## Install Dependencies

From the repository root:

```bash
pnpm install
pnpm exec playwright install --with-deps
```

## Python Environment

Backend-backed project tests require the local Python environment.

Preferred setup:

```bash
pnpm run python-setup
```

If the Python environment is already present, Playwright fixtures will use it.
They do not install or repair Python dependencies during a test run.

## Commands You Should Use

These are the supported commands.

List tests with the low-level runner:

```bash
pnpm run playwright-test --list
```

Run the supported catalog pack:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-catalog --reporter=list
```

Run the dev-only catalog pack:

```bash
pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test-catalog-dev --reporter=list
```

Run backend-backed project preflight:

```bash
pnpm run playwright-preflight-project
pnpm run playwright-preflight-project --diagnostic
```

Run the backend-backed project suite:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project --reporter=list
```

Run one backend-backed project spec:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/chart_creation.spec.ts --reporter=list
```

Run the supported mixed suite:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-all --reporter=list
```

## Commands You Should Not Treat As Default

`pnpm run playwright-test` is the low-level runner.

It is useful for:

- listing tests
- targeted manual runs
- debugging specific files

It is not the recommended default for the mixed suite because it does not
enforce the backend-backed worker policy. If you point it at everything, some
project tests may fail because they run with the wrong execution model.

## What Each Runner Does

- `playwright-test-catalog`
  - runs the supported catalog pack
  - defaults to Chromium
- `playwright-test-catalog-dev`
  - runs the dev-only catalog pack
  - defaults to Chromium
- `playwright-test-project`
  - runs project preflight first
  - defaults to `tests_playwright/project/`
  - defaults to Chromium
  - defaults to one worker
- `playwright-test-all`
  - runs the catalog runner first
  - runs the backend-backed project runner second

## Quick Troubleshooting

Check the backend:

```bash
curl -f http://localhost:5055
```

Run project preflight:

```bash
pnpm run playwright-preflight-project
pnpm run playwright-preflight-project --diagnostic
```

If the Python env is missing or broken:

```bash
pnpm run python-setup
pnpm run playwright-preflight-project --diagnostic
```

If Poetry or `.venv` recreation is blocked on macOS by ACL metadata:

```bash
ls -le python/.venv
chmod -N python/.venv
cd python
poetry install --with dev
```

To reset the local Docker database:

```bash
docker compose -f docker-secrets.yml down -v
docker compose -f docker-secrets.yml up -d
```

## Sandbox Note

Catalog tests are the best candidates for sandbox-first execution.

Backend-backed `tests_playwright/project/*` should be treated as unsandboxed by
default because they depend on:

- reliable access to `localhost:5055`
- real browser launch
- real backend project lifecycle
