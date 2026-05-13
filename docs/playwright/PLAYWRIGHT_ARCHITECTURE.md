# Playwright Architecture

This document explains how the Playwright test system is structured and how the
main fixtures work.

## Test Areas

- `tests_playwright/catalog/`
  - browser-only tests
  - route mocks
  - no real backend project lifecycle
- `tests_playwright/project/`
  - real backend-backed tests
  - use `tests_playwright/utils/projectFixtures`
  - create and clean up real temporary projects

## Main Modules

Runner layer:

- `scripts/playwright/run_playwright_cli.mjs`
- `scripts/playwright/playwright_catalog_runner.mjs`
- `scripts/playwright/playwright_catalog_ui_runner.mjs`
- `scripts/playwright/playwright_catalog_dev_runner.mjs`
- `scripts/playwright/playwright_catalog_dev_ui_runner.mjs`
- `scripts/playwright/playwright_project_runner.mjs`
- `scripts/playwright/playwright_project_ui_runner.mjs`
- `scripts/playwright/playwright_all_runner.mjs`
- `scripts/playwright/playwright_project_preflight.mjs`

Fixture layer:

- `tests_playwright/utils/projectFixtures/core.ts`
- `tests_playwright/utils/projectFixtures/syntheticAnndata.ts`
- `tests_playwright/utils/projectFixtures/syntheticSpatial.ts`
- `tests_playwright/utils/projectFixtures/importZip.ts`
- `tests_playwright/utils/projectFixtures/index.ts`

UI helper layer:

- `tests_playwright/utils/helpers.ts`

Legacy seam:

- `tests_playwright/utils/testUtils.ts`
- `tests_playwright/utils/tempProject.ts`

## Fixture Hierarchy

Use `tests_playwright/utils/projectFixtures`.

- `syntheticAnndata`
  - default backend-backed project fixture
  - generates a synthetic AnnData-backed project under `~/mdv`
  - triggers `/rescan_projects`
  - opens the project and waits for readiness
- `syntheticSpatial`
  - spatial-specific backend-backed fixture
  - use only when behavior really depends on spatial or image semantics
- `importZip`
  - import-workflow fixture only
  - not the default for general project tests

Legacy:

- `tests_playwright/utils/testUtils.ts` is not the canonical setup path for new
  backend-backed project specs
- `tests_playwright/utils/tempProject.ts` is compatibility re-export only

## Synthetic Project Lifecycle

Default backend-backed project creation flow:

1. generate project files under `~/mdv`
2. trigger `/rescan_projects`
3. locate the project in `/projects`
4. open `/project/:id`
5. wait until the project is ready

Default cleanup flow:

1. delete backend project row
2. remove generated files/folders

Retention override:

```bash
PLAYWRIGHT_KEEP_PROJECTS=1
```

## Execution Policy

- default lifecycle: one synthetic project per test
- supported optimization: explicit shared-project setup only for a small number
  of serial, read-mostly specs
- destructive or mutation-heavy specs should stay per-test
- backend-backed project runs stay Chromium-first and one-worker by default

## Shared-Project Optimization

The shared-project seam now exists in `syntheticAnndata.ts`.

Use it only when all of these are true:

- the spec is serial
- the spec is read-mostly
- each test still gets a fresh Playwright page
- shared backend state does not invalidate the assertions

Current shared-project example:

- `tests_playwright/project/charts/chart_creation.spec.ts`

## Test Type Meanings

Use these terms consistently when classifying project specs:

- `smoke`
  - minimal confidence that a core flow works at all
- `workflow`
  - normal multi-step user behavior
- `persistence`
  - verifies state after reload or reopen
- `mutation`
  - changes meaningful project state
- `destructive`
  - removes or invalidates state
- `error-path`
  - verifies failure handling and recovery

Current supported catalog pack:

- `tests_playwright/catalog/catalog_view.spec.ts`
- `tests_playwright/catalog/project_operations.spec.ts`

Current supported dev-only catalog pack:

- `tests_playwright/catalog/create_project.spec.ts`
- `tests_playwright/catalog/import_project.spec.ts`
