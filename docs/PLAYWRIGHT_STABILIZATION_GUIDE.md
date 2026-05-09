# Playwright Stabilization Guide

This guide is the living reference for phase-one Playwright stabilization.

## Current Goal

Build a small reliable Playwright foundation that works for both humans and
agents with the same command and helper conventions.

## Non-Goals (Phase One)

- Broad coverage across all charts and edge cases.
- Multi-browser parity as a hard requirement.
- Reworking the React Compiler probe into normal app test runs.

## Canonical Commands

Run from the repository root:

```bash
pnpm run playwright-test --list
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test tests_playwright/catalog/ --project=chromium --reporter=list
pnpm run playwright-preflight-project
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project
```

## Environment Prerequisites

- Dependencies installed (`pnpm install`).
- Playwright browser binaries installed locally.
- Backend available at `http://localhost:5055` for project tests.
- Vite dev server available for local catalog-only runs.
- Repo `venv` prepared and usable for backend-backed project tests.

## Preflight (Run Before Project Tests)

Use the dedicated preflight before running backend-backed project tests:

```bash
pnpm run playwright-preflight-project
pnpm run playwright-preflight-project --diagnostic
```

If these checks pass, run project tests from repo root:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project
```

## Decision Tree: Catalog vs Project

- Choose `tests_playwright/catalog/` when route-mock browser-only coverage is
  enough.
- Choose `tests_playwright/project/` when behavior depends on real project
  state, backend save/load, or reload persistence.

## Helper Selection Rules

Import from **`tests_playwright/utils/projectFixtures`** (barrel) or from the
specific module under `tests_playwright/utils/projectFixtures/`:

- Prefer **`syntheticAnndata.ts`** (`createTemporaryProjectViaSyntheticAnndata`) for default backend-backed project tests.
- Prefer **`syntheticSpatial.ts`** (`createTemporaryProjectViaSyntheticSpatial`) only for spatial-specific behavior; see `docs/SYNTHETIC_SPATIALDATA_PROJECTS.md`.
- Prefer **`importZip.ts`** (`createTemporaryProject`) only for import-flow coverage.
- Shared helpers live in **`core.ts`** and **`pythonEnv.ts`**. `tests_playwright/utils/tempProject.ts` is a **deprecated** re-export of the barrel; prefer `projectFixtures` in new code.
- `tests_playwright/utils/testUtils.ts` is legacy for older catalog-driven setup
  assumptions and should not be used for new backend-backed project tests.

## Browser Target

Chromium is the phase-one target (`--project=chromium`).

## Shared Backend Assumptions

- Current stabilization uses a shared backend and database.
- Project tests should isolate via temporary project creation and cleanup.
- Backend-backed project runs should stay on Chromium and one worker in phase one.

## Agent Escalation Notes

- In sandboxed runs, browser launch may require approval even if backend is up.
- Server startup can require approval in restricted sessions.

## Add A New Project Test

1. Create an isolated project with one of the three flows above (`importZip`,
   `syntheticAnndata`, `syntheticSpatial`, all exported from `projectFixtures`).
2. Use unique chart/view names (`Date.now()` suffix).
3. For the import helper, keep `allowCsvFallback: false` for canonical
   backend-backed suites unless the backend contract explicitly supports the
   fallback project-upload flow.
4. Assert against stable app state/config where possible.
5. Always call the handle `cleanup()` in `finally`.
6. Cleanup deletes both generated files and the backend project row by default.
7. Keep projects only by explicitly setting `PLAYWRIGHT_KEEP_PROJECTS=1`.

## Add A New Catalog Test

1. Use route mocks from `tests_playwright/utils/routes.ts`.
2. Keep flows browser-only and independent from backend project state.
3. Prefer focused navigation/assertion coverage over broad UI sweeps.

## Failure Diagnosis Checklist

Classify each failure before fixing:

- Wrong helper path.
- Stale selector.
- Stale setup assumption.
- Backend/environment issue.
- Python fixture-generation dependency issue (for example missing `numpy`,
  `scanpy`, or `h5py`).
- Browser permission issue.
- Actual app regression.

## Common Local Failures (And Fast Fixes)

### Backend is down

Symptoms:

- `TEST_BASE_URL=http://localhost:5055` runs fail immediately.
- `curl http://localhost:5055` does not return `200`.

Fix:

- Start/restart the backend container or local backend process.

### Python virtualenv is broken or missing

Symptoms:

- `projectFixtures` cannot generate mock project archives (Python/subprocess failure).
- errors include `python/.venv/bin/python ENOENT`, `h5py` missing, or `polars` missing.

Fix:

```bash
pnpm run python-setup
pnpm run playwright-preflight-project --diagnostic
```

### macOS ACL blocks Poetry from recreating `.venv`

Symptoms:

- Poetry reports virtualenv is broken and fails with `Permission denied` on `python/.venv`.
- `ls -le python/.venv` shows ACL entries such as `deny delete`.

Fix (local filesystem metadata only):

```bash
chmod -N python/.venv
cd python
poetry install --with dev
```

This ACL fix is local-only (not git-tracked) and does not change app/backend code.

### Playwright-created projects are being retained

Why this happens:

- You explicitly set `PLAYWRIGHT_KEEP_PROJECTS=1`.
- Cleanup failed after bounded retries and the test should surface that failure.

Normal phase-one behavior is to delete both the backend project row and
generated files by default.

## Suite Status

| Suite / File | Status | Notes |
| --- | --- | --- |
| `tests_playwright/catalog/catalog_view.spec.ts` | canonical candidate | Route-mock catalog coverage |
| `tests_playwright/catalog/create_project.spec.ts` | canonical candidate | Focused catalog flow |
| `tests_playwright/catalog/import_project.spec.ts` | legacy / needs review | Upload/socket caveats in comments |
| `tests_playwright/catalog/project_operations.spec.ts` | candidate / verify | Keep only if low-maintenance and stable |
| `tests_playwright/catalog/share_project.spec.ts` | legacy / partial | Contains TODO/incomplete behavior |
| `tests_playwright/project/chart_creation.spec.ts` | canonical candidate | Synthetic AnnData (`syntheticAnndata.ts`), one shared project + serial tests |
| `tests_playwright/project/chart_creation_single.spec.ts` | canonical candidate | Synthetic AnnData via `generate_synthetic_anndata_project.py` + rescan |
| `tests_playwright/project/chart_creation_huge.spec.ts` | optional / heavy | Large synthetic fixture; long timeout |
| `tests_playwright/project/view_management.spec.ts` | canonical candidate | Uses `projectFixtures` (zip import) |
| `tests_playwright/project/scatter_view_persistence.spec.ts` | canonical candidate | Persistence asserted via config/state |
| `tests_playwright/project/soft_delete.spec.ts` | canonical candidate | Strong backend-backed project coverage |
| `tests_playwright/utils/testUtils.ts` | legacy | Old catalog assumptions |
| `tests_playwright/utils/projectFixtures/` | canonical | `importZip`, `syntheticAnndata`, `syntheticSpatial`, `core`, `pythonEnv`; `tempProject.ts` re-exports for compatibility |

## Next Steps After Phase One

- Keep the canonical pack green first, then add targeted tests.
- Promote stable catalog candidates to canonical status after verification.
- Retire stale suites instead of patching outdated structures repeatedly.
