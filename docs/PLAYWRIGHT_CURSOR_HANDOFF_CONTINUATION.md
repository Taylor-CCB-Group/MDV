# Playwright Cursor Handoff (Continuation)

This document is a separate continuation handoff and does not replace
`docs/PLAYWRIGHT_CURSOR_HANDOFF.md`.

It captures:

- what has been completed in this chat
- what changed in fixture architecture in parallel work
- what failed, why, and what was fixed
- what remains unresolved
- exact next steps and commands

---

## 1) Current Goal

Stabilize a small, reliable Chromium-first Playwright foundation for:

- backend-backed `tests_playwright/project/*`
- browser-only route-mocked `tests_playwright/catalog/*`

Phase one is not broad coverage; it is reliability and repeatability.

---

## 2) Architecture State (Important)

Project fixtures are now centered on:

- `tests_playwright/utils/projectFixtures/importZip.ts`
- `tests_playwright/utils/projectFixtures/syntheticAnndata.ts`
- `tests_playwright/utils/projectFixtures/syntheticSpatial.ts`
- `tests_playwright/utils/projectFixtures/core.ts`
- `tests_playwright/utils/projectFixtures/pythonEnv.ts`

Compatibility layer:

- `tests_playwright/utils/tempProject.ts` is now a deprecated re-export barrel.

Legacy path:

- `tests_playwright/utils/testUtils.ts` is legacy and should not be used for new
  backend-backed project tests.

---

## 3) What Was Completed In This Chat

### 3.1 Project test rewrites and helper consistency

- `tests_playwright/project/chart_creation.spec.ts` rewritten from old setup.
- `tests_playwright/project/view_management.spec.ts` rewritten from old setup.
- `tests_playwright/project/scatter_view_persistence.spec.ts` added/maintained
  for config/state persistence.
- `tests_playwright/project/helpers.ts` used for shared project-view interactions.

### 3.2 Docs/workflow alignment

- Command-shape clarification made across docs:
  use `pnpm run playwright-test ...` (avoid accidental extra literal separator usage).
- Multiple docs were aligned to `projectFixtures` direction and stabilization workflow.

### 3.3 Type checks and lint checks

- TS follow-up checks run (`pnpm exec tsgo`) and fixed where needed during edits.
- Lints on touched files were clean during this chat’s edit cycles.

---

## 4) Failures Seen, With Classification

### A) Python dependency failures while generating mock project archives

Observed missing modules across runs:

- `h5py`
- `numpy`
- `scanpy`
- `polars`

Classification: backend test fixture environment dependency issue.

### B) CSV fallback path instability

When fallback route engaged, failures included:

- `/create_project` not OK in that backend mode
- upload flow mismatch (`input[type="file"]` not attached)

Classification: helper/backend-contract mismatch.

### C) Teardown delete timeout

Later run failures showed multiple tests timing out in cleanup:

- `apiRequestContext.delete("/delete_project/:id")` timing out under test timeout

Classification: backend/environment teardown reliability issue.

### D) Soft delete dialog behavior drift/flakiness

Failure observed:

- `Add Column` dialog expected to close after `Add`, but remained visible.

Classification: selector/state behavior drift (or async UI timing issue).

---

## 5) Fixes Applied (or Attempted) During This Chat

### Applied

- Installed required Python packages in `venv` for fixture generation.
- Installed missing `polars` after explicit failure surfaced it.
- Smoke-tested `create_test_project_zip(...)` successfully after dependency fix.
- Updated/kept project specs on generated-mock canonical path rather than CSV fallback.
- Added/import-hardening changes in Python package init/test factory to reduce
  eager failures on optional deps.

### Still unresolved after latest user-provided run

- teardown delete endpoint timing out in multiple tests
- soft delete first test dialog close behavior

---

## 6) Current Risk Notes

1. **Shared backend state risk**: if cleanup deletes are flaky and default cleanup
   policy avoids deletion, stale rows/projects accumulate and can affect future runs.
2. **Fixture path split risk**: parallel work introduced synthetic filesystem+rescan
   fixtures; canonical usage needs explicit suite-by-suite choice and documentation.
3. **Long test timeout masking root cause**: some failures are in teardown, not test
   assertions; this can obscure actual pass/fail semantics of test body logic.

---

## 7) Canonical Commands (Current)

From repo root.

List chromium tests:

```bash
pnpm run playwright-test --list --project=chromium
```

Canonical project pack:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test tests_playwright/project/chart_creation.spec.ts tests_playwright/project/view_management.spec.ts tests_playwright/project/scatter_view_persistence.spec.ts tests_playwright/project/soft_delete.spec.ts --project=chromium --reporter=list
```

Minimal catalog pack:

```bash
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test tests_playwright/catalog/catalog_view.spec.ts tests_playwright/catalog/create_project.spec.ts --project=chromium --reporter=list
```

TS follow-up when TS files are touched:

```bash
pnpm exec tsgo
```

---

## 8) What Is Left To Do

1. **Stabilize cleanup policy in fixtures**:
   - decide if canonical tests should default to `deleteOnCleanup: false`
   - if so, document DB cleanup procedure as explicit post-run step
   - if not, debug and harden `/delete_project` behavior in tests
2. **Fix soft_delete first test flake**:
   - verify if dialog requires additional state transition wait
   - validate button action completion signal before asserting dialog hidden
3. **Re-run project canonical pack** until stable.
4. **Run small catalog pack** and classify failures.
5. **Update suite status table** in stabilization guide with observed, not intended, status.

---

## 9) Suggested Next-Agent Prompt

```md
Continue Playwright stabilization using:
- docs/PLAYWRIGHT_CURSOR_HANDOFF_CONTINUATION.md
- docs/PLAYWRIGHT_STABILIZATION_GUIDE.md

Goals:
1) make canonical project pack stable
2) fix soft_delete add-column dialog failure
3) verify canonical catalog subset
4) update suite status table based on actual runs

Constraints:
- Chromium only
- run commands from repo root
- use tests_playwright/utils/projectFixtures as canonical fixture layer
- do not use tests_playwright/utils/testUtils.ts for backend-backed project specs
```

