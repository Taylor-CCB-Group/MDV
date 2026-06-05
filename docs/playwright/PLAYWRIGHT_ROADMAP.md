# Playwright Roadmap

Planned improvements for the Playwright test setup. Not implemented yet — use this to align
PRs and design discussions.

Current docs for day-to-day work:

- `PLAYWRIGHT_GUIDE.md` — human setup and run
- `PLAYWRIGHT_AGENT.md` — agent rules and prompt templates
- `PLAYWRIGHT_ARCHITECTURE.md` — fixtures and execution policy
- `PLAYWRIGHT_WRITING_TESTS.md` — how to add tests

---

## 1. Connect CI/CD to Playwright tests

**Goal:** Run the suite automatically on pull requests and main, with clear pass/fail signal
and artifacts on failure.

**Directions to explore:**

- Split jobs by tier: catalog-only (mocked, lighter) vs project tests (Docker backend, Python,
  browsers).
- Reuse existing wrapper scripts (`playwright-test-catalog`, `playwright-test-project`) so CI
  matches local commands.
- Set `TEST_BASE_URL`, `PLAYWRIGHT_HTML_OPEN=never`, and install browsers in the job (default
  Playwright location only).
- Decide when project tests run (every PR vs label/nightly) given Docker + Python + runtime cost.
- Publish traces/screenshots/HTML report as CI artifacts on failure.

**Open questions:**

- Which specs are required to run?
- Shared backend/DB strategy in CI (fresh stack per job vs reusable service).

---

## 2. More test coverage

**Goal:** Cover critical user flows and regressions without one-off manual checks.

**Directions to explore:**

- Prioritize high-value project flows: chart creation, views, column/table actions, persistence
  after reload.
- Keep catalog tests for route-mocked UI that does not need a real project lifecycle.
- One spec per feature area under `tests_playwright/project/<area>/`, starting from
  `_template.spec.ts`.
- Use `createTemporaryProjectViaSyntheticAnndata` by default; shared-project suites only for
  serial read-mostly cases (see `PLAYWRIGHT_ARCHITECTURE.md`).

**Open questions:**

- Minimum bar for new UI features (spec required vs follow-up issue).

---

## 3. Reduce test run time

**Goal:** Faster feedback locally and in CI without sacrificing reliability.

**Directions to explore:**

- Keep project runs Chromium + one worker where backend state is shared; avoid widening
  parallelism until isolation is proven.
- Expand justified shared-project suites (serial, read-mostly) instead of one project per test
  when safe.
- Trim synthetic data size (`minimal` profile, smaller `nCells`/`nGenes`) where assertions allow.
- Run a focused subset during development (`<spec>` or `--grep`) before full suite.
- Cache Playwright browsers and pnpm store in CI; avoid redundant `playwright-setup` steps.
- Measure baseline times per spec/file and target the slowest tests first.

**Open questions:**

- Whether catalog and project suites can run in parallel jobs vs sequential `playwright-test-all`.
- Trade-offs of skipping preflight in CI when env is pre-provisioned.

---

## 4. Better use of `playwright-cli` for exploration and writing tests

**Goal:** Make exploration a reliable step before committing specs, not a dead end in
`artifacts/playwright-cli/`.

**Directions to explore:**

- Standard workflow: explore → snapshot/locators → fresh spec from `_template.spec.ts` →
  validate with `pnpm run playwright-test-project` (documented in `PLAYWRIGHT_AGENT.md`).
- Use `--debug=cli` + `playwright-cli attach` for failing tests (see `PLAYWRIGHT_GUIDE.md`).
- Optional seed specs or tagged flows for common setups (e.g. scatter + legend) once patterns
  stabilize.
- Skill at `.agents/skills/playwright-cli/` for command reference; keep repo rules in
  `PLAYWRIGHT_AGENT.md`.
- Do not commit scratch `run-code` scripts; gitignore `artifacts/playwright-cli/` remains scratch
  only.

**Open questions:**

- Whether to add thin repo helpers promoted from repeated exploration (only when used by ≥2 specs).
- MCP/browser tooling vs `playwright-cli` — single recommended path for agents.

---

## How to update this doc

When a roadmap item ships, move the outcome into `PLAYWRIGHT_GUIDE.md`, `PLAYWRIGHT_AGENT.md`,
or `PLAYWRIGHT_ARCHITECTURE.md` and trim or close the section here.
