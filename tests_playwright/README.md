# Playwright Tests

Run Playwright from the repository root so it uses `playwright.config.ts`.
For the current phase-one stabilization model and suite status table, see
`docs/PLAYWRIGHT_STABILIZATION_GUIDE.md`.

## Quick Start

```bash
pnpm install
pnpm exec playwright install --with-deps
docker compose -f docker-secrets.yml up -d
pnpm run playwright-test-project
```

The default target is the devcontainer app on `http://localhost:5055`. Override
it with `TEST_BASE_URL`.

For local Vite-only catalog tests:

```bash
pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test tests_playwright/catalog/catalog_view.spec.ts --project=chromium --reporter=list
```

The pnpm scripts call the local Playwright dependency directly. Prefer them over
plain `pnpm exec playwright` for test execution, especially in restricted agent
environments.

## Common Commands

```bash
pnpm run playwright-test --list
pnpm run playwright-preflight-project
pnpm run playwright-test-project
pnpm run playwright-test tests_playwright/catalog/
pnpm run playwright-test --headed
pnpm run playwright-test-ui
```

## Test Areas

- `catalog/`: catalog view tests with mocked backend APIs.
- `project/`: project view tests that need a real backend.
- `utils/`: shared routes, fixtures, and helpers.
- `test-data/`: small checked-in test inputs.

Project tests cover chart creation and view-management flows. They require a
backend that can create projects and process uploads.

Three setup flows live under **`tests_playwright/utils/projectFixtures/`** (import
from `projectFixtures` or from `importZip.ts` / `syntheticAnndata.ts` /
`syntheticSpatial.ts` directly). `tests_playwright/utils/tempProject.ts` re-exports
the barrel for older imports.

- **`syntheticAnndata.ts`** — **default backend-backed project fixture**.
  **`createTemporaryProjectViaSyntheticAnndata`**:
  runs `python -m mdvtools.tests.generate_synthetic_anndata_project`, then
  `GET /rescan_projects` and `/projects` id diff. Cleanup removes the generated
  folder under `~/mdv`.
- **`syntheticSpatial.ts`** — **spatial-specific backend-backed fixture**.
  **`createTemporaryProjectViaSyntheticSpatial`**:
  runs `python -m mdvtools.tests.generate_synthetic_spatial_project`, then the
  same rescan/id-diff flow. Matches `docs/SYNTHETIC_SPATIALDATA_PROJECTS.md`.
- **`importZip.ts`** — **import-flow fixture only**. **`createTemporaryProject`**:
  mock `.mdv.zip` via Python `create_test_project_zip`, then `POST /import_project`
  (optional CSV fallback). Issues `PUT /projects/<id>/access` for DB editable mode
  only; disk `state.json` is **not** rewritten by that route.

Cleanup is enabled by default. Playwright-created backend projects delete both
the backend project row and generated files unless the operator explicitly sets
`PLAYWRIGHT_KEEP_PROJECTS=1`.

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
pnpm run playwright-test --reporter=list
```

For backend-backed project tests, run this preflight first:

```bash
pnpm run playwright-preflight-project
pnpm run playwright-preflight-project --diagnostic
```

If Poetry reports `.venv` is broken with `Permission denied` on macOS, check for
an ACL deny-delete entry and clear it:

```bash
ls -le python/.venv
chmod -N python/.venv
cd python
poetry install --with dev
```

`chmod -N` here only updates local filesystem ACL metadata. It is not git-tracked
and does not change application code.

To reset the local Docker database:

```bash
docker compose -f docker-secrets.yml down -v
docker compose -f docker-secrets.yml up -d
```
