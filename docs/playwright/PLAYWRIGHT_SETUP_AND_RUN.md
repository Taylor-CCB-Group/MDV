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

Run all the tests:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-all
```

Run the catalog test suite:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-catalog
```

Run the dev-only catalog pack:

```bash
pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test-catalog-dev
```

Open the catalog test suite in UI mode:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-catalog-ui
```

Run the backend-backed project preflight:

```bash
pnpm run playwright-preflight-project
pnpm run playwright-preflight-project --diagnostic
```

Run the backend-backed project suite:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project
```

Open the backend-backed project suite in UI mode:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project-ui
```

Run one backend-backed project spec:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/charts/chart_creation.spec.ts
```

`playwright-test` now delegates to `playwright-test-all` and is safe to use as
the default mixed-suite command. Do not use `run_playwright_cli.mjs` directly
as the normal entrypoint; it does not enforce the backend-backed worker policy.

## Suppressing the HTML Report

By default Playwright opens an HTML report after a run. To suppress it:

```bash
PLAYWRIGHT_HTML_OPEN=never pnpm run playwright-test
PLAYWRIGHT_HTML_OPEN=never pnpm run playwright-test-project --reporter=list
```

This is useful in CI and LLM-driven runs where opening a browser tab is
unwanted.

## Debugging a Failing Test

Run the failing spec with `--debug=cli`. The test pauses at the start and
prints a session name:

```bash
PLAYWRIGHT_HTML_OPEN=never pnpm run playwright-test-project tests_playwright/project/charts/chart_creation.spec.ts --project=chromium --reporter=list -- --debug=cli
```

Wait until output includes `Debugging Instructions` and a session name like
`tw-abcdef`, then in a second terminal attach `playwright-cli` to that session:

```bash
pnpm exec playwright-cli attach tw-abcdef
```

From there you can take snapshots, inspect elements, and step through the test.
Every action you perform generates the equivalent Playwright TypeScript code in
the output, which you can copy directly into the spec.

After fixing the test, stop the background runner and rerun to verify.

## Optional `playwright-cli`

`playwright-cli` is available in this repo as an optional browser exploration
and debugging tool. It is not the canonical test runner.

Use it when you want to:

- inspect a local page interactively
- discover stable locators
- reproduce a flow before writing a real test
- debug a paused Playwright test session (see above)

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

Keep using the repo wrapper scripts for actual suite execution:

- `pnpm run playwright-test`
- `pnpm run playwright-test-catalog`
- `pnpm run playwright-test-catalog-dev`
- `pnpm run playwright-test-project`

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
- `playwright-test` / `playwright-test-all`
  - equivalent: both delegate to `playwright_all_runner.mjs`
  - runs the catalog runner first
  - runs the backend-backed project runner second (one worker, Chromium)
- `playwright-test-ui`
  - raw pass-through to Playwright UI mode (`--ui`) with no worker or spec
    constraints; intended only for interactive exploration, not as a
    substitute for the structured runners above

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

If `uv` or `.venv` recreation is blocked on macOS by ACL metadata:

```bash
ls -le python/.venv
chmod -N python/.venv
cd python
uv sync --group dev --frozen
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