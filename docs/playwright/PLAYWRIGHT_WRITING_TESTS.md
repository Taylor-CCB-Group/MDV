# Playwright Writing Tests

This document explains how to add or update Playwright tests in this repository.

## General Rule

Add tests incrementally.

Do not batch-generate many new files at once. Prefer one target file or one
feature area at a time, then verify it before moving on.

## Optional `playwright-cli` Workflow

`playwright-cli` is useful for writing and debugging tests, but it does not
replace the repo fixture architecture or runner scripts.

Recommended workflow:

1. explore the page with `pnpm exec playwright-cli`
2. use `snapshot`, `click`, `fill`, `eval`, and generated locator/code output
   to understand the UI
3. move the useful parts into a real `@playwright/test` spec
4. adapt them to repo helpers and fixtures
5. validate with the repo wrapper scripts

Example:

```bash
pnpm exec playwright-cli open http://localhost:5055 --headed
pnpm exec playwright-cli snapshot
pnpm exec playwright-cli click e12
pnpm exec playwright-cli --raw generate-locator e12
pnpm exec playwright-cli close
```

Treat the generated code as draft material, not final test code.

## Writing Catalog Tests

Use catalog tests when the behavior is:

- browser-side
- route-mock-friendly
- not about a real backend project lifecycle

Guidelines:

1. Use route mocks from `tests_playwright/utils/routes.ts`
2. Keep the test independent from `projectFixtures`
3. Prefer focused interaction checks over broad end-to-end sweeps
4. If the behavior depends on `import.meta.env.PROD === false`, keep it in the
   dev-only catalog path

## Writing Project Tests

Use project tests when the behavior depends on:

- a real backend-backed project
- view state
- chart state
- persistence after reload
- meaningful project mutation

Guidelines:

1. Choose the correct fixture family:
   - `syntheticAnndata` for most backend-backed project behavior
   - `syntheticSpatial` only for spatial-specific behavior
   - `importZip` only for import-flow coverage
2. Use unique names for charts and views
3. Prefer config/state assertions where possible
4. Always make cleanup explicit through the fixture handle
5. Do not use the old catalog-driven setup helpers

## Choosing The Lifecycle Model

Use these rules:

- `smoke`
  - isolated by default
  - may become shared if there are several related read-mostly tests
- `workflow`
  - can be shared only if it does not persist or destructively mutate backend state
- `persistence`
  - keep isolated
- `mutation`
  - keep isolated
- `destructive`
  - keep isolated
- `error-path`
  - keep isolated

Use shared-project optimization only when:

- the spec is serial
- the spec is read-mostly
- each test still gets a fresh page
- runtime improvement is measurable

## Recommended Reuse Pattern

Prefer this layering:

- `projectFixtures` owns project lifecycle
- `tests_playwright/utils/helpers.ts` owns repeated UI actions
- specs describe behavior only

Good reuse:

- `createSharedSyntheticAnndataSuite(...)`
- `addChartViaUi(...)`
- `createViewViaUi(...)`
- `selectViewViaUi(...)`

Avoid repeating:

- project bootstrap logic
- rescan/delete logic
- repeated locators across many specs

## LLM Guidance

When using an LLM to write tests:

- give one file or one feature area at a time
- provide explicit test cases
- use `playwright-cli` when the model needs to inspect the live UI, find
  locators, or debug a paused browser session
- prefer extending existing helpers over duplicating setup logic
- prefer the documented wrapper scripts for actual test execution, not raw
  Playwright CLI commands
- do not ask the model to repair Python environments during normal test execution

LLMs should use `playwright-cli` for:

- exploration
- locator discovery
- paused-test debugging

LLMs should use repo scripts for:

- final pass/fail validation
- catalog runs
- backend-backed project runs
- mixed-suite runs

## Anti-Patterns

Avoid:

- using `tests_playwright/utils/testUtils.ts` for new backend-backed project specs
- mixing smoke checks and destructive flows in one file
- sharing one backend project across mutation-heavy tests
- treating the raw low-level Playwright runner as the default mixed-suite command
- treating `playwright-cli` as the suite runner
- adding broad speculative coverage before verifying the previous file

## When To Update Other Docs

Update `docs/playwright/PLAYWRIGHT_ARCHITECTURE.md` when:

- a new fixture seam is added
- execution policy changes
- spec classification changes

Update `docs/playwright/PLAYWRIGHT_SETUP_AND_RUN.md` when:

- a supported command changes
- setup prerequisites change
- the default run path changes
