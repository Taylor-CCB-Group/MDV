- only run `pnpm exec tsgo` when TypeScript files or TS-facing configs changed
- for Python checks in worktrees, prefer the local `venv`:
- tests: `./venv/bin/python -m pytest python/mdvtools/tests -m "not performance"`
- pyright: from `python/`, run `../venv/bin/pyright`
- avoid `as` casts where possible
- prefer follow-up PRs for wider type-system cleanups
- when composing commit messages, avoid words like "enhance" and favour conciseness and clarity
- for Playwright checks, run from the repo root with the local scripts so `playwright.config.ts` is loaded; use `pnpm run playwright-test ...` as the low-level runner and `pnpm run playwright-test-project ...` for backend-backed `tests_playwright/project/*` runs
- for phase-one Playwright stabilization, target Chromium (`--project=chromium`) and keep backend-backed project runs on one worker
- backend-backed project tests should use `tests_playwright/utils/projectFixtures` and default to cleanup; keep projects only with `PLAYWRIGHT_KEEP_PROJECTS=1`
- in sandboxed agent sessions, expect local server startup and browser launch to require approval, and avoid plain `pnpm exec playwright` for test execution because it may try the npm registry even when dependencies are installed
- for phase-one Playwright stabilization, target Chromium (`--project=chromium`) and keep catalog tests route-mocked while backend-backed project tests use `tests_playwright/utils/projectFixtures` (or the deprecated re-export `tests_playwright/utils/tempProject.ts`)
- do not use `pnpm run playwright-test -- ...`; pass args directly to avoid injecting a literal `--` into Playwright
- `tests_playwright/utils/testUtils.ts` is legacy for old catalog-style setup assumptions and should not be used for new backend-backed `tests_playwright/project/` specs

## Agent skills

### Issue tracker

Issues are tracked in GitHub Issues for this repository. See `.agents/issue-tracker.md`.

### Triage labels

Triage uses the default labels: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, and `wontfix`. See `.agents/triage-labels.md`.

### Domain docs

Domain docs are configured as a single-context layout. See `.agents/domain.md`.
