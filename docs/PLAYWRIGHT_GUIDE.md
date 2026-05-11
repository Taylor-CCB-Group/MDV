# Playwright Guide

This is the main source of truth for Playwright setup, execution, fixture
choice, and troubleshooting in this repository.

## Goal

Phase one is a small, reliable Chromium-first Playwright foundation for both
manual runs and LLM-driven runs.

Non-goals for phase one:

- broad coverage across all chart types and edge cases
- multi-browser parity
- folding React Compiler probe work into normal app-test runs

## Test Areas

- `tests_playwright/catalog/`
  - browser-only tests
  - route mocks
  - no real backend project lifecycle
- `tests_playwright/project/`
  - real backend-backed tests
  - use `tests_playwright/utils/projectFixtures`
  - create and clean up real temporary projects

## Canonical Commands

Run all commands from the repository root.

List tests:

```bash
pnpm run playwright-test --list
```

Catalog tests against local Vite:

```bash
pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test tests_playwright/catalog/ --project=chromium --reporter=list
```

Project-test preflight:

```bash
pnpm run playwright-preflight-project
pnpm run playwright-preflight-project --diagnostic
```

Canonical backend-backed project run:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project
```

Run one backend-backed project spec:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/chart_creation_single.spec.ts --reporter=list
```

`playwright-test-project` performs light preflight first, then defaults to:

- `tests_playwright/project/`
- `--project=chromium`
- `--workers=1`

## Environment Prerequisites

- `pnpm install`
- `pnpm exec playwright install --with-deps`
- backend running at `http://localhost:5055` for backend-backed project tests
- local Vite dev server for catalog-only Vite runs
- prepared repo Python environment for backend-backed project fixtures

The project fixtures assume the Python environment already exists. They do not
install or repair dependencies during the test run.

## Fixture Hierarchy

Use `tests_playwright/utils/projectFixtures`.

- `syntheticAnndata.ts`
  - default backend-backed project fixture
  - generates a synthetic AnnData-backed project under `~/mdv`
  - rescans backend project registry
- `syntheticSpatial.ts`
  - spatial-specific backend-backed fixture
  - use only when behavior truly depends on spatial/image/region semantics
- `importZip.ts`
  - import-workflow fixture only
  - not the default for new project tests

Legacy:

- `tests_playwright/utils/testUtils.ts` is legacy and should not be used for new
  backend-backed project tests
- `tests_playwright/utils/tempProject.ts` is compatibility re-export only

## Cleanup Policy

Backend-backed Playwright projects clean up by default.

Cleanup should remove:

- the backend project row
- generated files/folders

Keep projects only when explicitly requested:

```bash
PLAYWRIGHT_KEEP_PROJECTS=1 TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/chart_creation_single.spec.ts --reporter=list
```

## Sandbox vs Unsandboxed Execution

Use sandboxed execution for:

- code exploration
- writing or editing tests
- lightweight static checks

Use unsandboxed execution for:

- backend-backed `tests_playwright/project/*`
- browser launches that need host permissions
- anything that needs reliable access to `localhost:5055`

Catalog tests are the best candidates for sandbox-first execution.

Project tests should be treated as unsandboxed by default.

## How To Add A New Project Test

1. Choose the correct fixture family:
   - `syntheticAnndata` for normal backend-backed project behavior
   - `syntheticSpatial` only for spatial-specific behavior
   - `importZip` only for import-flow coverage
2. Use unique names for charts and views.
3. Assert against stable config/state where possible.
4. Always call `cleanup()` in `finally`.
5. Do not use the old catalog-driven setup helpers.

## How To Add A New Catalog Test

1. Use route mocks from `tests_playwright/utils/routes.ts`
2. Keep it independent from backend project lifecycle
3. Prefer focused behavior checks over broad UI sweeps

## Failure Classification

Classify failures before fixing them:

- wrong helper path
- stale selector
- stale setup assumption
- backend/environment issue
- Python fixture dependency issue
- browser permission issue
- actual app regression

Do not jump straight to selector edits.

## Troubleshooting

Backend check:

```bash
curl -f http://localhost:5055
```

Project preflight:

```bash
pnpm run playwright-preflight-project
pnpm run playwright-preflight-project --diagnostic
```

If Python env is missing or broken:

```bash
pnpm run python-setup
pnpm run playwright-preflight-project --diagnostic
```

If Poetry/venv recreation is blocked on macOS by ACL metadata:

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

## Suite Status

| Suite / File | Status | Notes |
| --- | --- | --- |
| `tests_playwright/catalog/catalog_view.spec.ts` | canonical candidate | Route-mock catalog coverage |
| `tests_playwright/catalog/create_project.spec.ts` | canonical candidate | Focused catalog flow |
| `tests_playwright/catalog/import_project.spec.ts` | legacy / needs review | Upload/socket caveats in comments |
| `tests_playwright/catalog/project_operations.spec.ts` | candidate / verify | Keep only if stable and low-maintenance |
| `tests_playwright/catalog/share_project.spec.ts` | legacy / partial | TODO/incomplete behavior |
| `tests_playwright/project/chart_creation.spec.ts` | canonical candidate | Synthetic AnnData pattern |
| `tests_playwright/project/chart_creation_single.spec.ts` | canonical candidate | Minimal backend-backed project test |
| `tests_playwright/project/scatter_view_persistence.spec.ts` | canonical candidate | Persistence via config/state |
| `tests_playwright/project/view_management.spec.ts` | canonical candidate | Project-view behavior |
| `tests_playwright/project/soft_delete.spec.ts` | canonical candidate | Backend-backed mutation/error coverage |
| `tests_playwright/utils/projectFixtures/` | canonical | Main backend-backed fixture layer |
| `tests_playwright/utils/testUtils.ts` | legacy | Old catalog-driven assumptions |

## React Compiler Probe

Keep React Compiler probe work separate from normal Playwright stabilization.

Use the existing compiler probe flow only when you intentionally want compiler
off/on smoke validation.
