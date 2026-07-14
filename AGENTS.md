- only run `pnpm exec tsgo` when TypeScript files or TS-facing configs changed
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

### Spatial image channel controls

For the SpatialData chart image layer panel (`ImageLayerPanel`, `ColorChannelComponents` spatial branch, `image_layer_runtime.ts`):

- Contrast limits, tone, and histogram brush must update **immediately** — do not debounce these writes.
- Persistence is MobX `renderStack` entry props (`channels`, `vivLayerProps`) via `useLayerChannelState` / `patchLayer`; the canvas reads those props through the Render Stack Adapter.
- Per-panel zustand (`VivProvider`) is a one-way UI projection only, not a second persistence path.
- If controls lag or the histogram brush fights itself, fix state ownership/wiring (see `docs/spatialdata-vis-integration.md` § Image channel controls). Debouncing here masks design problems.
- Do not bump `renderStackGeneration` on in-place `props` patches (channels, tone, opacity); that recreates `vivPassthrough` and stalls the spatial renderer.
