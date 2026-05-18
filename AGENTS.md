- only run `pnpm exec tsgo` when TypeScript files or TS-facing configs changed
- when TypeScript files change, run `pnpm exec biome lint` on touched paths only (read-only; no `--write`); do not run `jsformat`, `jslint`, or `biome check --write` on `./src` — formatting the tree would touch many unrelated lines; do not remove `biome-ignore` comments without reading the rationale and confirming the underlying constraint is still addressed
- for Python checks in worktrees, prefer the local `venv`:
- tests: `./venv/bin/python -m pytest python/mdvtools/tests -m "not performance"`
- pyright: from `python/`, run `../venv/bin/pyright`
- avoid `as` casts where possible
- prefer follow-up PRs for wider type-system cleanups
- when composing commit messages, avoid words like "enhance" and favour conciseness and clarity
- for any Playwright work, read `docs/playwright/PLAYWRIGHT_AGENT.md` first — it has rules, commands, and prompt templates for agents

## Agent skills

### Issue tracker

Issues are tracked in GitHub Issues for this repository. See `.agents/issue-tracker.md`.

### Triage labels

Triage uses the default labels: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, and `wontfix`. See `.agents/triage-labels.md`.

### Domain docs

Domain docs are configured as a single-context layout. See `.agents/domain.md`.
