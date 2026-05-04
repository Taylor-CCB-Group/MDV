# Playwright Tests

Run Playwright from the repository root so it uses `playwright.config.ts`.

## Quick Start

```bash
npm install
npx playwright install --with-deps
docker compose -f docker-secrets.yml up -d
npm run playwright-test
```

The default target is the devcontainer app on `http://localhost:5055`. Override
it with `TEST_BASE_URL`.

For local Vite-only catalog tests:

```bash
npm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 npm run playwright-test -- tests_playwright/catalog/catalog_view.spec.ts --project=chromium --reporter=list
```

The npm scripts call the local Playwright dependency directly. Prefer them over
plain `npx playwright` for test execution, especially in restricted agent
environments.

## Common Commands

```bash
npm run playwright-test -- --list
npm run playwright-test -- tests_playwright/catalog/
npm run playwright-test -- tests_playwright/project/
npm run playwright-test -- --headed
npm run playwright-test-ui
```

## Test Areas

- `catalog/`: catalog view tests with mocked backend APIs.
- `project/`: project view tests that need a real backend.
- `utils/`: shared routes, fixtures, and helpers.
- `test-data/`: small checked-in test inputs.

Project tests cover chart creation and view-management flows. They require a
backend that can create projects and process uploads.

For agent-driven project tests, prefer the shared helper in
`tests_playwright/utils/tempProject.ts`. It creates a temporary project for the
test, imports it through the backend, and cleans it up afterward. The preferred
path uses the Python mock-project generator; if that environment is not
available locally, the helper can fall back to a temporary inline CSV seed so
tests can still exercise the real project page without checked-in one-off test
fixtures.

## Generated Data

Test MDV projects can be generated from Python:

```bash
python -m mdvtools.tests.generate_test_data ~/mdv/test_mock --mock
python -m mdvtools.tests.generate_test_data ~/mdv/test_pbmc3k --scanpy pbmc3k_processed
python -m mdvtools.tests.generate_test_data ~/mdv/test_large --mock --n-cells 1000000 --n-genes 5000
```

For programmatic use, import from
`mdvtools.tests.test_project_factory`. For SpatialData-backed performance
projects, see `docs/SYNTHETIC_SPATIALDATA_PROJECTS.md`.

## Troubleshooting

```bash
docker compose -f docker-secrets.yml ps mdv_app
docker compose -f docker-secrets.yml logs mdv_app
curl http://localhost:5055
npm run playwright-test -- --reporter=list
```

To reset the local Docker database:

```bash
docker compose -f docker-secrets.yml down -v
docker compose -f docker-secrets.yml up -d
```
