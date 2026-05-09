# Playwright Agent Workflow

This repository's Playwright tests live in `tests_playwright/` and are
configured by the root `playwright.config.ts`. Run Playwright commands from the
repository root so that config is loaded.

For phase-one stabilization scope and suite status tracking, use
`docs/PLAYWRIGHT_STABILIZATION_GUIDE.md` as the primary living reference.

## Commands

List tests:

```bash
pnpm run playwright-test --list
```

Project-test preflight:

```bash
pnpm run playwright-preflight-project
pnpm run playwright-preflight-project --diagnostic
```

Run browser-only catalog tests against a local Vite server:

```bash
pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test tests_playwright/catalog/ --project=chromium --reporter=list
```

Run project tests against a backend that can create projects and process
uploads:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/chart_creation.spec.ts --reporter=list
```

The pnpm scripts call the local Playwright install. In restricted agent sessions,
avoid plain `pnpm exec playwright` for test execution because it may try the npm
registry even when dependencies are installed.

`playwright-test-project` runs light preflight first, then defaults to
`tests_playwright/project/`, `--project=chromium`, and `--workers=1` unless you
override them.

## React Compiler Smoke Check (PR C)

Use this when validating React Compiler rollout safety without enabling it
globally. Keep `VITE_USE_REACT_COMPILER` opt-in and compare off/on runs against
the same project URL.

```bash
# compiler off
VITE_USE_REACT_COMPILER=0 pnpm run build-flask-vite
node scripts/playwright_compiler_probe.mjs http://127.0.0.1:5055/project/191

# compiler on
VITE_USE_REACT_COMPILER=1 pnpm run build-flask-vite
node scripts/playwright_compiler_probe.mjs http://127.0.0.1:5055/project/191
```

Compare probe output for:

- `issues` (console warnings/errors, page errors)
- `chartManagerReadyMs`, `firstChartsMs`, and `chartsSettledMs`
- `summary` payload parity (`chartTypes`, `chartCount`, `rowCount`)

If payload composition differs between runs, treat the timing comparison as
smoke-only and not a strict performance benchmark.

## Sandboxed Agent Notes

Expect approval to be needed for:

- starting local servers, because binding to localhost can fail with
  `listen EPERM`;
- launching browsers, because macOS browser startup can fail with permission
  errors.

Use the smallest meaningful check first. Catalog tests mock the project API in
the browser with `page.route`; project tests create or upload real projects and
need a live backend.

Agents should validate the prepared repo `venv` and backend first. Do not try
to install or repair Python dependencies during a normal Playwright run.

## Agentic Project Testing

For project-view features, prefer the shared helper in
`tests_playwright/utils/projectFixtures` instead of checking in ad hoc project
fixtures.

Recommended pattern:

1. generate a temporary project with the canonical `projectFixtures` adapter;
2. register it with the backend if needed;
3. run the Playwright flow against the real project page;
4. delete both the backend project row and generated files in teardown.

Phase-one fixture hierarchy:

- `syntheticAnndata.ts` is the default backend-backed project fixture
- `syntheticSpatial.ts` is for spatial-specific behavior
- `importZip.ts` is for import-flow coverage only

Preflight for project tests (recommended before any debugging):

```bash
pnpm run playwright-preflight-project
pnpm run playwright-preflight-project --diagnostic
```

If Poetry cannot recreate `python/.venv` with `Permission denied` on macOS,
check for ACL deny-delete entries and clear them:

```bash
ls -le python/.venv
chmod -N python/.venv
cd python
poetry install --with dev
```

Use browser-only route mocks for `catalog/` work. Use `projectFixtures` for
`project/` work that depends on datasource state, chart manager state, save
behaviour, or reload persistence.

## Human Setup

For the Docker/devcontainer path:

```bash
pnpm install
pnpm exec playwright install --with-deps
docker compose -f docker-secrets.yml up -d
pnpm run playwright-test
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
