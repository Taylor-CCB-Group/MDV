# Playwright Agent Workflow

This repository's Playwright tests live in `tests_playwright/` and are
configured by the root `playwright.config.ts`. Run Playwright commands from the
repository root so that config is loaded.

## Commands

List tests:

```bash
npm run playwright-test -- --list
```

Run browser-only catalog tests against a local Vite server:

```bash
npm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 npm run playwright-test -- tests_playwright/catalog/ --project=chromium --reporter=list
```

Run project tests against a backend that can create projects and process
uploads:

```bash
TEST_BASE_URL=http://localhost:5055 npm run playwright-test -- tests_playwright/project/ --project=chromium --reporter=list
```

The npm scripts call the local Playwright install. In restricted agent sessions,
avoid plain `npx playwright` for test execution because it may try the npm
registry even when dependencies are installed.

## Sandboxed Agent Notes

Expect approval to be needed for:

- starting local servers, because binding to localhost can fail with
  `listen EPERM`;
- launching browsers, because macOS browser startup can fail with permission
  errors.

Use the smallest meaningful check first. Catalog tests mock the project API in
the browser with `page.route`; project tests create or upload real projects and
need a live backend.

## Human Setup

For the Docker/devcontainer path:

```bash
npm install
npx playwright install --with-deps
docker compose -f docker-secrets.yml up -d
npm run playwright-test
```

For local Vite-only catalog work, use the Vite command above and set
`TEST_BASE_URL`.

## Mock Data

There are two useful tracks:

- browser-only catalog mocks in `tests_playwright/utils/routes.ts` and
  `tests_playwright/utils/data.ts`;
- real generated MDV projects from Python, using
  `python/mdvtools/tests/generate_test_data.py`,
  `python/mdvtools/tests/test_project_factory.py`, or the synthetic
  SpatialData generator documented in
  `docs/SYNTHETIC_SPATIALDATA_PROJECTS.md`.

From a worktree, prefer the local Python environment:

```bash
./venv/bin/python -m pytest python/mdvtools/tests -m "not performance"
cd python
../venv/bin/python -m mdvtools.tests.generate_test_data /tmp/mdv-test-project --mock --force
```

Keep generated artifacts out of git unless there is a specific reason to check
in a fixture.
