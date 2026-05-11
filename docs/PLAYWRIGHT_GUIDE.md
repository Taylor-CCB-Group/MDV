# Playwright Guide

This is the main source of truth for Playwright setup, execution, fixture
choice, and troubleshooting in this repository.

## Goal

The current goal is a small, reliable Chromium-first Playwright suite that is
safe to run, easy to classify, and easy to extend incrementally for both manual
and LLM-driven workflows.

Non-goals for the current phase:

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

List tests with the low-level runner:

```bash
pnpm run playwright-test --list
```

Canonical catalog-only run against local Vite:

```bash
pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test-catalog --reporter=list
```

`playwright-test-catalog` defaults to the currently supported catalog pack, not
every file under `tests_playwright/catalog/`.

Dev-only catalog run for create/import cards:

```bash
pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test-catalog-dev --reporter=list
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

Canonical supported-suite run:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-all --reporter=list
```

Run one backend-backed project spec:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/chart_creation_single.spec.ts --reporter=list
```

`playwright-test-catalog` defaults to:

- the supported catalog pack
- `--project=chromium`

`playwright-test-catalog-dev` defaults to:

- the dev-only catalog pack
- `--project=chromium`

`playwright-test-project` performs light preflight first, then defaults to:

- `tests_playwright/project/`
- `--project=chromium`
- `--workers=1`

`playwright-test-all` runs:

- the catalog runner first
- the backend-backed project runner second

Use `pnpm run playwright-test` only for targeted/manual use. It is not the safe
default for the mixed suite because it does not enforce the backend-backed
worker policy.

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

## Project Execution Policy

- Default lifecycle: one synthetic project per test
- Supported optimization: explicit shared-project setup only for a small number
  of serial, read-mostly specs where setup cost dominates and state leakage is
  controlled
- Keep destructive or mutation-heavy specs per-test
- Keep backend-backed project runs on Chromium and one worker by default

## Project Spec Statefulness

Use these terms to choose the lifecycle model for backend-backed project specs:

- `smoke`
  - minimal confidence that a core project flow works at all
  - usually cheap and a good candidate for shared-project reuse if read-mostly
- `workflow`
  - normal multi-step user behavior
  - may be shared if it does not persist or destructively mutate backend state
- `persistence`
  - verifies state after reload or reopen
  - usually keep isolated per test
- `mutation`
  - changes meaningful project state
  - usually keep isolated per test
- `destructive`
  - removes or invalidates project state
  - keep isolated per test
- `error-path`
  - verifies failures and recovery
  - keep isolated per test

Current project spec classification:

| File | Statefulness | Lifecycle |
| --- | --- | --- |
| `tests_playwright/project/chart_creation_single.spec.ts` | smoke | isolated per test |
| `tests_playwright/project/chart_creation.spec.ts` | workflow, read-mostly | shared serial project |
| `tests_playwright/project/view_management.spec.ts` | workflow, persistence | isolated per test |
| `tests_playwright/project/soft_delete.spec.ts` | mutation, destructive, error-path | isolated per test |

## How To Add A New Project Test

1. Choose the correct fixture family:
   - `syntheticAnndata` for normal backend-backed project behavior
   - `syntheticSpatial` only for spatial-specific behavior
   - `importZip` only for import-flow coverage
2. Use unique names for charts and views.
3. Assert against stable config/state where possible.
4. Always call `cleanup()` in `finally`.
5. Do not use the old catalog-driven setup helpers.
6. Treat new coverage as incremental: add one file or feature area at a time
   from explicit test cases instead of broad speculative generation.

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
| `tests_playwright/catalog/catalog_view.spec.ts` | canonical | Stable route-mock catalog coverage; verified in the catalog runner |
| `tests_playwright/catalog/create_project.spec.ts` | canonical (dev-only) | Requires local Vite where `import.meta.env.PROD` is false; now permission-aware |
| `tests_playwright/catalog/import_project.spec.ts` | canonical (dev-only) | Requires local Vite where `import.meta.env.PROD` is false; currently dialog-level coverage |
| `tests_playwright/catalog/project_operations.spec.ts` | canonical | Info, rename, and delete verified; export removed because the UI entry is not currently rendered |
| `tests_playwright/catalog/share_project.spec.ts` | rewrite-next | Partial coverage with TODOs and incomplete interaction assertions |
| `tests_playwright/project/chart_creation.spec.ts` | canonical | Shared serial project via `projectFixtures`; read-mostly chart workflow coverage |
| `tests_playwright/project/chart_creation_single.spec.ts` | canonical | Minimal backend-backed smoke path for the synthetic AnnData fixture |
| `tests_playwright/project/soft_delete.spec.ts` | canonical | High-value mutation coverage that should stay per-test |
| `tests_playwright/project/view_management.spec.ts` | canonical | Workflow plus persistence checks; keep isolated per test |
| `tests_playwright/utils/projectFixtures/` | canonical | Main backend-backed fixture layer |
| `tests_playwright/utils/testUtils.ts` | retire | Old catalog-driven assumptions; do not use for new backend-backed specs |

Use these classifications literally:

- `canonical`: supported now and expected to stay healthy
- `optimization candidate`: supported now, but may justify explicit lifecycle optimization work
- `rewrite-next`: keep the behavior area, but rewrite the spec shape before expanding it
- `legacy`: kept only as a historical or temporary path; do not expand
- `retire`: not part of the supported future path

Current supported catalog pack:

- `tests_playwright/catalog/catalog_view.spec.ts`
- `tests_playwright/catalog/project_operations.spec.ts`

Current supported dev-only catalog pack:

- `tests_playwright/catalog/create_project.spec.ts`
- `tests_playwright/catalog/import_project.spec.ts`

## React Compiler Probe

Keep React Compiler probe work separate from normal Playwright stabilization.

Use the existing compiler probe flow only when you intentionally want compiler
off/on smoke validation.
