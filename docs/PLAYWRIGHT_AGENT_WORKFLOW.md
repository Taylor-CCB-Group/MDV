# Playwright Agent Workflow

This note captures how Playwright currently works in this repository from both
human and agent sessions.

## Current State

- The test suite lives in `tests_playwright/` and is configured by the root
  `playwright.config.ts`.
- `@playwright/test` and `playwright` are installed locally in `node_modules`.
- The default `baseURL` is `http://localhost:5055/`, intended for the
  devcontainer or Docker app.
- `TEST_BASE_URL` overrides the target app URL.
- Catalog tests mock the project API in the browser with `page.route`.
- Project tests create real projects, upload `tests_playwright/test-data/scanpy-pbmc3k.h5ad`,
  and require a live backend.
- CI is present in `.github/workflows/playwright.yml`, but it is manual only and
  currently marked as needing review.

## Verified From This Agent Session

These checks were run on April 28, 2026:

```bash
node node_modules/playwright/cli.js --version
```

Result: Playwright `1.57.0`.

```bash
npm run playwright-test -- --list
```

Result: 117 browser test entries were discovered: 39 tests across catalog and
project specs, multiplied across Chromium, Firefox, and WebKit.

The following command started Vite successfully when run outside the filesystem
sandbox:

```bash
npm run dev -- --host 127.0.0.1 --port 5173
```

The following command launched Chromium successfully when run outside the
filesystem sandbox:

```bash
TEST_BASE_URL=http://127.0.0.1:5173 \
  node node_modules/playwright/cli.js test tests_playwright/catalog/catalog_view.spec.ts \
  --project=chromium --reporter=list
```

Observed result: two catalog tests passed and the existing `theme toggle` test
failed because the page remained on `dark` after clicking the toggle. That
failure is an application/test assertion issue, not a Playwright availability
issue.

## Agent Use

Prefer root-level invocations so the config is loaded:

```bash
npm run playwright-test -- --list
TEST_BASE_URL=http://127.0.0.1:5173 npm run playwright-test -- tests_playwright/catalog/ --project=chromium --reporter=list
```

In this Codex workspace, expect two actions to need approval outside the
sandbox:

- starting local servers, because binding to `127.0.0.1` or `0.0.0.0` can fail
  with `listen EPERM`;
- launching Chromium, because macOS browser startup can fail with
  `MachPortRendezvousServer ... Permission denied`.

Avoid plain `npx playwright` in restricted sessions. It can attempt to resolve
the package through `registry.npmjs.org` even when local packages are installed.
Use `npm run playwright-test` or `node node_modules/playwright/cli.js ...`
instead.

## Human Use

For the Docker/devcontainer path:

```bash
npm install
npx playwright install --with-deps
docker compose -f docker-secrets.yml up -d
npm run playwright-test
```

For local Vite-only catalog tests:

```bash
npm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 npm run playwright-test -- tests_playwright/catalog/ --project=chromium --reporter=list
```

Use project tests only against a backend that can create projects and process
uploads:

```bash
TEST_BASE_URL=http://localhost:5055 npm run playwright-test -- tests_playwright/project/ --project=chromium --reporter=list
```

## Mock Data Projects

There are two useful mock-data tracks:

1. Browser-only catalog mocks: add `page.route` mocks in `tests_playwright/utils/routes.ts`
   and small fixtures in `tests_playwright/utils/data.ts`.
2. Real MDV projects: generate AnnData or MDV project zips from Python using
   `python/mdvtools/tests/generate_test_data.py` or
   `python/mdvtools/tests/test_project_factory.py`.

From a worktree, prefer the local Python environment:

```bash
./venv/bin/python -m pytest python/mdvtools/tests -m "not performance"
cd python
../venv/bin/python -m mdvtools.tests.generate_test_data /tmp/mdv-test-project --mock --force
../venv/bin/python -m mdvtools.tests.generate_test_data /tmp/mdv-test-unique --mock --with-unique-column --force
```

Use generated projects for repeatable upload/import/project-view flows, and keep
large generated artifacts out of git unless there is a specific reason to check
in a fixture.

## Iteration Plan

1. Fix the catalog theme toggle test so it asserts the current app behavior or
   waits for the intended theme transition.
2. Add one minimal smoke command to CI that only lists tests and validates the
   Playwright config before running browsers.
3. Rework the manual CI workflow so reports upload from the correct container
   path and the workflow can be run predictably on pull requests when ready.
4. Add a small generated mock-project fixture workflow for import/project tests,
   ideally produced on demand from `test_project_factory.py`.
5. Split fast browser-only catalog tests from backend/project tests in docs and
   CI so agents can run the smallest meaningful check first.
