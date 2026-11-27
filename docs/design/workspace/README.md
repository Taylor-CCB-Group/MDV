# Workspace redesign

Start here. This folder contains the plan and evidence for bringing the
[workspace prototype](prototype.html) into MDV: one workspace view containing charts from
several data sources, responsive layout, a docked inspector and one colour-token system.

The prototype is the visual starting point, not a decision on the lasting layout model,
persistence format, project/data lifecycle or final mobile experience.

## Documents

| Document | Purpose |
|---|---|
| [PLAN.md](PLAN.md) | Showcase scope, architecture guardrails, implementation sequence, checkpoints and open decisions |
| [RESEARCH.md](RESEARCH.md) | Consolidated inventory of current layout, view, chart-chrome and theme code |
| [LAYOUT-EVALUATION.md](LAYOUT-EVALUATION.md) | Layout-lab evidence and the next grid/tiling/hybrid comparison |
| [prototype.html](prototype.html) | Executable visual reference; open it directly in a browser |
| [ADR 0026](../../adr/0026-preview-one-view-workspace-saved-to-views-v2.md) | Preview boundary, in-memory model and candidate persistence |
| [ADR 0027](../../adr/0027-build-the-workspace-grid-on-the-prototypes-own-engine.md) | In-house engine as the initial baseline behind a replaceable seam |

## Destination

A real-project showcase and the evidence needed to choose the workspace architecture. If the
team agrees with the direction, the next destination is a spec (PRD), ADRs for the durable
decisions and an ordered set of build-slice issues.

## Current status

- **Research completed:** the layout candidates; every current view-layout reader and writer;
  chart chrome, settings and filters; theme and colour ownership; two layout bake-off rounds.
- **Agreed direction:** build a real-project showcase behind a browser preview; choose the
  preview path in `project_bootstrap`; reuse current charts through a compatibility adapter;
  do not make `ChartManager` own the new workspace; leave the pane workspace unchanged.
- **Working choices:** use the in-house grid first because it already hosts real charts;
  deterministically convert pane views into an in-memory workspace model; prioritise desktop
  and responsive laptop/tablet behaviour, with a safe basic phone fallback.
- **Open decisions:** scrolling grid, tiling or hybrid; canonical saved schema and Python
  derivation; longer-term project/data lifecycle; specialised inspector panels; filters,
  theming, undo and rollout.

## Decision history

- ADR 0026 records the preview preference, deterministic pane-view conversion and a candidate
  flat workspace model. Review reopened the durable format, canonical schema, Python
  derivation and detailed responsive rules.
- ADR 0027 records the in-house engine as the showcase baseline behind a narrow seam. Its
  gridstack and react-grid-layout findings stand; review reopened scrolling grid versus tiling
  or a hybrid.
- New durable decisions get new ADRs. An accepted ADR is superseded rather than silently
  rewritten. The working plan can be edited while PR #719 is under review.

## Route

1. Trace `project_bootstrap` and #706; define workspace-session, legacy-chart-adapter and future
   project-data boundaries without implementing `datasources_v2`.
2. Build the bounded first slice in [PLAN.md](PLAN.md): mount real charts through that boundary,
   use the in-house grid and convert pane views in memory.
3. Run the [layout checkpoint](LAYOUT-EVALUATION.md) with shared Playwright scenarios and human
   trials; choose grid, tiling or hybrid.
4. Agree the canonical schema, Python derivation and durable format; then implement save/reload.
5. Settle chart chrome, inspector, tools, filters, theme, undo and rollout with the showcase
   running.
6. Team verdict, then the PRD and independently pickable build slices.

## Working rules

- Work on `docs/workspace-redesign-plan` (PR #719). Nothing merges to `main` until the team has
  tried the showcase and the layout and persistence checkpoints are resolved.
- ADR 0019 applies: extract legacy UI into React in place; restyling is a separate change.
- Keep engine details behind the layout seam and find cards in tests by roles and
  `data-chart-id`, not engine classes.
- A workspace session owns view, layout, selection, dirty state and persistence. Reach existing
  chart machinery through the smallest legacy adapter; do not add workspace orchestration to
  `ChartManager`.
- Preserve the chart-hosting findings in [LAYOUT-EVALUATION.md](LAYOUT-EVALUATION.md).
- Keep project names and data from development servers out of this public repository.
- Every implementation must work on the hosted Flask/Docker server. Static-build support is a
  rollout decision and must be settled before the legacy workspace retires.

## Running the lab

- Front end: `pnpm install`, then `pnpm run dev` (Vite on port 5170, or set `PORT=…`).
- Backend: Flask is expected at `http://127.0.0.1:5055`; run it normally, for example in Docker.
- Open a project with `?layoutLab=inhouse` or `?layoutLab=gridstack` in a development build.
- Checks are in `src/react/prototypes/layoutLab/checks/`. Read
  `docs/playwright/PLAYWRIGHT_AGENT.md` before Playwright work.
- Use three kinds of private test project: spatial/WebGL; a dense view (the lab used 23 charts);
  and a view over three data sources with mixed old pane layouts. Never commit their names.

## Out of scope for the showcase

- Catalog-page redesign.
- A full phone workspace beyond a safe basic fallback.
- The desktop app (`mdv_desktop`) and embedding MDV in other pages.

Before changing the workspace, read this file, [PLAN.md](PLAN.md), ADRs 0026 and 0027, and open
[prototype.html](prototype.html). Then explain what is agreed, what is only a working baseline
and which checkpoint the proposed work belongs to.
