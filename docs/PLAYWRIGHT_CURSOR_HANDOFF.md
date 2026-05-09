# Playwright Cursor Handoff (Updated)

This is the current handoff artifact for continuing Playwright stabilization with
minimal rediscovery. It supersedes earlier assumptions that project fixtures are
only `tempProject.ts`.

---

## 1) Current Goal

Phase-one goal remains: establish a small, reliable Chromium-first Playwright
foundation for this repo, not broad coverage.

Success criteria:

- canonical project pack is stable enough to run repeatedly
- canonical catalog pack is validated
- docs reflect actual helper architecture and failure modes

---

## 2) Major Direction Change Since Prior Handoff

Backend fixture setup has been expanded/migrated to:

- `tests_playwright/utils/projectFixtures/importZip.ts`
- `tests_playwright/utils/projectFixtures/syntheticAnndata.ts`
- `tests_playwright/utils/projectFixtures/syntheticSpatial.ts`
- `tests_playwright/utils/projectFixtures/core.ts`
- `tests_playwright/utils/projectFixtures/pythonEnv.ts`

`tests_playwright/utils/tempProject.ts` is now a deprecated re-export barrel.

This means future work should treat `projectFixtures` as canonical and
`tempProject.ts` as compatibility only.

---

## 3) What Has Been Completed

### 3.1 Project-test rewrite foundation

- `tests_playwright/project/chart_creation.spec.ts` rewritten away from
  `newProjectSetup` to backend fixture flow.
- `tests_playwright/project/view_management.spec.ts` rewritten away from
  `newProjectSetup`.
- `tests_playwright/project/scatter_view_persistence.spec.ts` added and aligned
  to config/state persistence assertions.
- `tests_playwright/project/helpers.ts` added for shared project-view UI helpers.
- `tests_playwright/utils/testUtils.ts` marked legacy for old catalog assumptions.

### 3.2 Docs aligned to updated command shape and helper direction

- `AGENTS.md`, `docs/PLAYWRIGHT_AGENT_WORKFLOW.md`, `tests_playwright/README.md`
  updated to avoid `pnpm run playwright-test -- ...` argument misuse.
- living guide exists at `docs/PLAYWRIGHT_STABILIZATION_GUIDE.md`.

### 3.3 Environment blockers resolved incrementally

The following missing Python deps were hit and addressed during stabilization:

- `h5py`
- `numpy`
- `scanpy`
- `polars`

Additionally, `mdvtools` import behavior was made more tolerant for missing
optional dependencies:

- `python/mdvtools/__init__.py` now avoids immediate hard-fail for optional
  heavy imports.
- `python/mdvtools/tests/test_project_factory.py` lazy-loads `scanpy`.

---

## 4) Issues Encountered, Classification, and Resolution Status

### Issue A — Wrong Playwright invocation shape in docs/prompts (resolved)

- Symptom: `No tests found` / bad arg parsing when using extra literal `--`.
- Classification: command usage/documentation issue.
- Resolution: docs updated to direct shape, e.g.
  `pnpm run playwright-test tests_playwright/project/...`.

### Issue B — Stale `newProjectSetup()` pattern (resolved by rewrite)

- Symptom: setup fails at catalog assumptions (`Create new project` UI mismatch).
- Classification: stale setup assumption.
- Resolution: project specs migrated to backend fixture helpers.

### Issue C — Missing Python deps blocked zip fixture generation (partially resolved)

- Symptom: `Unable to generate a mock MDV project archive...` with missing modules.
- Classification: environment dependency issue.
- Resolution status: many deps installed/fixed, but this remains sensitive to
  the exact Python env used by test subprocesses. `polars` was the latest blocker.

### Issue D — CSV fallback path unreliable against current backend contract (open/avoid)

- Symptom: `/create_project` and upload flow mismatches; fallback caused failures.
- Classification: helper-path/backend-contract mismatch.
- Current stance: canonical suites should use generated mock/archive path;
  CSV fallback should remain disabled unless backend contract is known-good.

### Issue E — Cleanup delete API hanging and timing out tests (active)

- Symptom: multiple tests timeout while waiting on `request.delete("/delete_project/:id")`.
- Classification: backend/environment issue (teardown path).
- Current mitigation in newer fixture work: prefer `deleteOnCleanup: false` by
  default in zip/synthetic fixture handles; clean temp local files/folders and
  defer DB row cleanup when backend delete is flaky.
- Status: needs re-verification on latest branch state.

### Issue F — `soft_delete` dialog flow regression/flakiness (active)

- Symptom: `Add Column` dialog remains visible after clicking `Add`.
- Classification: possible selector/state drift or async UI behavior change.
- Status: unresolved in latest failing run snapshot.

---

## 5) Latest Known Failing Surface (from provided terminal logs)

When running canonical project pack, major failures were:

- chart creation + scatter persistence + view management timing out, often during
  teardown delete call
- soft delete first test failing to close `Add Column` dialog

One view-management test (`loads default view`) passed in the same run, indicating
core app load path is reachable and not globally broken.

---

## 6) Current Helper Model (Important for Next Agent)

Use `tests_playwright/utils/projectFixtures` directly:

- `importZip.ts`
  - generates mock `.mdv.zip` via Python factory
  - imports through `/import_project`
  - supports optional CSV fallback (not recommended for canonical pack)
- `syntheticAnndata.ts`
  - generates project on filesystem (`~/mdv`)
  - triggers `/rescan_projects`
  - discovers new id via `/projects` diff
- `syntheticSpatial.ts`
  - same filesystem+rescan orchestration for spatial synthetic projects

`tests_playwright/utils/tempProject.ts` is now re-export only.

---

## 7) Canonical Commands

Run from repo root.

List tests:

```bash
pnpm run playwright-test --list --project=chromium
```

Project canonical pack:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test tests_playwright/project/chart_creation.spec.ts tests_playwright/project/view_management.spec.ts tests_playwright/project/scatter_view_persistence.spec.ts tests_playwright/project/soft_delete.spec.ts --project=chromium --reporter=list
```

Catalog pack:

```bash
TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test tests_playwright/catalog/catalog_view.spec.ts tests_playwright/catalog/create_project.spec.ts --project=chromium --reporter=list
```

TypeScript follow-up when TS changed:

```bash
pnpm exec tsgo
```

---

## 8) What Still Needs To Be Done

1. Re-run canonical project pack on the latest `projectFixtures`-based state and
   confirm whether teardown delete timeouts are eliminated with current defaults.
2. Stabilize `soft_delete` add-column interaction (`dialog` close behavior).
3. Verify small catalog canonical pack.
4. Update suite status table in `docs/PLAYWRIGHT_STABILIZATION_GUIDE.md` with
   actual verified statuses, not intended statuses.
5. Add short cleanup playbook for orphaned test projects/rows if
   `deleteOnCleanup` remains default-off in shared backend.

---

## 9) Suggested Failure Classification Rubric (Keep Using This)

For each failing run, classify first:

- wrong helper path
- stale selector
- stale setup assumption
- backend/environment issue
- browser permission issue
- actual app regression

Do not jump straight to selector edits.

---

## 10) Suggested Next Prompt

```md
Continue Playwright stabilization using `docs/PLAYWRIGHT_CURSOR_HANDOFF.md` and `docs/PLAYWRIGHT_STABILIZATION_GUIDE.md`.

Focus:
1) make canonical project pack pass/reliably fail only on real regressions
2) fix soft_delete add-column dialog flake
3) verify catalog canonical pack
4) update suite status table with observed results

Constraints:
- Chromium only
- run from repo root with `pnpm run playwright-test ...`
- use `tests_playwright/utils/projectFixtures` as canonical backend fixture layer
- do not use `tests_playwright/utils/testUtils.ts` for backend-backed project specs
```

