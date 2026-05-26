# Playwright History

This file preserves the important background behind the current Playwright
setup, fixture architecture, and migration path. It is not the day-to-day guide.

For current usage, use:

- `docs/playwright/PLAYWRIGHT_GUIDE.md`
- `docs/playwright/PLAYWRIGHT_ARCHITECTURE.md`
- `docs/playwright/PLAYWRIGHT_WRITING_TESTS.md`

## Why The Current Architecture Exists

### Old project setup was stale

Earlier backend-backed project tests depended on catalog-page assumptions in
`tests_playwright/utils/testUtils.ts`, including a visible `Create new project`
control. As the frontend evolved, those assumptions broke before tests ever
reached the project page.

That led to a deliberate move away from catalog-driven setup toward backend
fixture helpers under `tests_playwright/utils/projectFixtures`.

### Synthetic project generation became preferred

The project originally leaned on import-style helpers and UI upload paths.
Those flows were more brittle than necessary for general backend-backed project
testing.

The fixture architecture expanded to include:

- `importZip.ts`
- `syntheticAnndata.ts`
- `syntheticSpatial.ts`
- `core.ts`
- `pythonEnv.ts`

The resulting hierarchy is:

- `syntheticAnndata` for default backend-backed project tests
- `syntheticSpatial` for spatial-specific behavior
- `importZip` only for import-flow coverage

### Cleanup policy changed

At one stage, fixture cleanup was frequently left off because backend deletion
was flaky and could hang test teardown. That produced leaked projects and extra
manual cleanup burden.

The agreed direction became:

- cleanup on by default
- delete both backend row and generated files
- keep projects only with explicit opt-out (`PLAYWRIGHT_KEEP_PROJECTS=1`)
- bounded retry around backend deletion

### Sandbox is not the canonical project-test execution path

Backend-backed Playwright project tests hit real host backend state and launch
real browsers. In this repo, sandboxed runs can fail because:

- `localhost:5055` is not reachable from sandbox
- browser startup is blocked by sandbox/host permission rules

That is why the current recommendation is:

- sandbox for writing/checking
- unsandboxed execution for real backend-backed project runs

## Problems Encountered During Stabilization

### Python environment issues

Several stabilization attempts were blocked by missing or broken Python deps,
including:

- `h5py`
- `numpy`
- `scanpy`
- `polars`

This led to the conclusion that fixture code should not attempt to repair the
Python environment during test runs. Instead, project preflight should validate
the environment and fail fast.

### CSV fallback instability

The fallback project-upload path proved unreliable against the current backend
contract:

- `/create_project` mismatches
- upload/input attachment mismatches

That is why generated synthetic project paths became the preferred backend test
substrate.

### UI details discovered during test writing

Some important current UI/testing details:

- Add Chart flows may require a title before submission is enabled
- config/state assertions after reload are often more stable than visual-only
  checks

These were not always obvious from the older tests and needed to be rediscovered
during rewrite work.

## Migration Notes

### Rewritten/realigned project specs

The main backend-backed project specs were moved toward the new fixture model:

- `chart_creation_single.spec.ts`
- `chart_creation.spec.ts`
- `view_management.spec.ts`
- `scatter_view_persistence.spec.ts`
- `soft_delete_column.spec.ts`

The direction of travel was:

- away from stale catalog-driven setup
- toward fixture-owned project lifecycle
- toward synthetic AnnData as the default backend-backed substrate

### Historical handoff notes

Two earlier handoff files captured intermediate state during the migration:

- `docs/PLAYWRIGHT_CURSOR_HANDOFF.md`
- `docs/PLAYWRIGHT_CURSOR_HANDOFF_CONTINUATION.md`

Their useful content has been condensed here. They are no longer the primary
reference for the current setup.

### Prompt guidance history

`docs/PLAYWRIGHT_PROMPT_GUIDANCE.md` recorded lessons from earlier LLM runs:

- old runner/path confusion
- sandbox/browser/backend constraints
- stale helper pitfalls

Those practical lessons are now represented in the current Playwright docs.

## Historical Decisions Worth Remembering

- `tests_playwright/utils/testUtils.ts` is legacy, not canonical
- `tests_playwright/utils/tempProject.ts` is compatibility re-export only
- backend-backed project tests should use `projectFixtures`
- synthetic AnnData is the default backend-backed test substrate
- project-test execution should be Chromium-first and one-worker in phase one
- a dedicated project preflight + project runner is better than expecting users
  or LLMs to improvise setup checks

## Remaining Kinds Of Work

The remaining work after the major architecture shift is no longer broad design.
It is mostly:

- verifying rewritten specs
- validating cleanup behavior in practice
- validating teammate/onboarding setup
- keeping suite status current
