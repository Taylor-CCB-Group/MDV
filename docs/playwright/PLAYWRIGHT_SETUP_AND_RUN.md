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

Optional browser-exploration tooling for humans and LLMs:

```bash
pnpm exec playwright-cli --help
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

Run the supported catalog pack:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-catalog --reporter=list
```

Run the dev-only catalog pack:

```bash
pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test-catalog-dev --reporter=list
```

Open the supported catalog pack in UI mode:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-catalog-ui
```

Open the dev-only catalog pack in UI mode:

```bash
pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test-catalog-dev-ui
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

Open the backend-backed project suite in UI mode:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project-ui
```

Run one backend-backed project spec:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/charts/chart_creation.spec.ts --reporter=list
```

Run the supported mixed suite:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-all --reporter=list
```

Do not use the raw low-level Playwright command as the normal entrypoint for
this repo. It does not enforce the backend-backed worker policy, and if you run
everything through it, project test failures are expected.

## Optional `playwright-cli`

`playwright-cli` is available in this repo as an optional browser exploration
and debugging tool. It is not the canonical test runner.

Use it when you want to:

- inspect a local page interactively
- discover stable locators
- reproduce a flow before writing a real test
- debug a paused Playwright test session

Run it through the local dependency:

```bash
pnpm exec playwright-cli open http://localhost:5055 --headed
pnpm exec playwright-cli snapshot
pnpm exec playwright-cli click e12
pnpm exec playwright-cli eval "document.title"
pnpm exec playwright-cli close
```

This repo configures `playwright-cli` to write its session artifacts under:

```text
artifacts/playwright-cli/
```

That keeps YAML snapshots and session logs out of the repo root.

When you want to keep a screenshot or snapshot, give it an explicit filename in
that directory:

```bash
pnpm exec playwright-cli screenshot --filename=artifacts/playwright-cli/screenshots/example.png
pnpm exec playwright-cli snapshot --filename=artifacts/playwright-cli/snapshots/example.yaml
```

When debugging a Playwright test paused with `--debug=cli`, attach to the named
session:

```bash
pnpm exec playwright-cli attach tw-abcdef
```

Keep using the repo wrapper scripts for actual suite execution:

- `pnpm run playwright-test-catalog`
- `pnpm run playwright-test-catalog-dev`
- `pnpm run playwright-test-project`
- `pnpm run playwright-test-all`

## What Each Runner Does

- `playwright-test-catalog`
  - runs the supported catalog pack
  - defaults to Chromium
- `playwright-test-catalog-ui`
  - opens the supported catalog pack in Playwright UI
- `playwright-test-catalog-dev`
  - runs the dev-only catalog pack
  - defaults to Chromium
- `playwright-test-catalog-dev-ui`
  - opens the dev-only catalog pack in Playwright UI
- `playwright-test-project`
  - runs project preflight first
  - defaults to `tests_playwright/project/`
  - defaults to Chromium
  - defaults to one worker
- `playwright-test-project-ui`
  - runs project preflight first
  - opens the backend-backed project suite in Playwright UI
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

The same rule usually applies to `playwright-cli` when you point it at the
local backend.
