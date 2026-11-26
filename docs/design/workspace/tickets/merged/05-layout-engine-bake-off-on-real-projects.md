Type: prototype
Status: merged
Merged into: [19 Layout engine bake-off and choice](../19-layout-engine-bake-off-and-choice.md)
Blocked by: [01 Layout engine candidates](../01-layout-engine-candidates.md), [02 Everything that reads or writes view layouts](../02-everything-that-reads-or-writes-view-layouts.md)
Map: [Workspace in MDV - Map](../../README.md)

## Question

Using the shortlist from [01 Layout engine candidates](../01-layout-engine-candidates.md) and the hosting contract from [02 Everything that reads or writes view layouts](../02-everything-that-reads-or-writes-view-layouts.md), extend the layout lab so each shortlisted engine hosts the real charts of at least three real projects in one view: one with spatial/WebGL charts, one with many charts, one with several data sources.

Measure and show the team:

- frame rate and chart redraws while dragging and resizing (`perf.ts` for spatial charts)
- layout quality at 1440, 1024, 768 and 390 px
- saving and reloading a layout
- what breaks

## Context

Layout lab: worktree `MDV-ui-refresh`, `src/react/prototypes/layoutLab`, opened with `?layoutLab=grid` or `?layoutLab=tiles`. Measurement script: `artifacts/playwright-cli/layout-lab/measure.mjs`.
