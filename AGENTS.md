- only run `npx tsgo` when TypeScript files or TS-facing configs changed
- for Python checks in worktrees, prefer the local `venv`:
- tests: `./venv/bin/python -m pytest python/mdvtools/tests -m "not performance"`
- pyright: from `python/`, run `../venv/bin/pyright`
- avoid `as` casts where possible
- prefer follow-up PRs for wider type-system cleanups
- when composing commit messages, avoid words like "enhance" and favour conciseness and clarity
- for Playwright checks, run from the repo root with the local script (`npm run playwright-test -- ...`) so `playwright.config.ts` is loaded; in sandboxed agent sessions, expect local server startup and browser launch to require approval, and avoid plain `npx playwright` for test execution because it may try the npm registry even when dependencies are installed

## Agent skills

### Issue tracker

Issues are tracked in GitHub Issues for this repository. See `.agents/issue-tracker.md`.

### Triage labels

Triage uses the default labels: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, and `wontfix`. See `.agents/triage-labels.md`.

### Domain docs

Domain docs are configured as a single-context layout. See `.agents/domain.md`.
