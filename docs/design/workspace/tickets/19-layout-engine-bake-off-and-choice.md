Type: prototype
Status: open
Blocked by: [18 One view, saved layout and screen sizes](18-one-view-saved-layout-and-screen-sizes.md)
Map: [Workspace in MDV - Map](../README.md)
Merges: [05 Layout engine bake-off on real projects](merged/05-layout-engine-bake-off-on-real-projects.md), [06 Choose the layout engine](merged/06-choose-the-layout-engine.md)

## Question

Which layout engine does MDV adopt: gridstack 14 or the in-house engine? Decide it by watching both run, not on paper.

react-grid-layout is out: [01 Layout engine candidates - findings](../research/01-layout-engine-candidates-findings.md) found it copies widths to narrower breakpoints without scaling them, which caused the gaps and staircases seen in the layout lab.

Extend the layout lab so both engines host the real charts of at least three real projects in one view: one with spatial/WebGL charts, one with many charts and one with several data sources. Measure, and show the human side by side:

- frame rate and chart redraws while dragging and resizing (`perf.ts` for spatial charts)
- layout quality at each width settled in [18 One view, saved layout and screen sizes](18-one-view-saved-layout-and-screen-sizes.md) (at least 1440, 1024, 768 and 390 px)
- saving and reloading a layout in the format settled there
- what breaks

The human picks the engine. Record the choice and what it gives up as an ADR.

## Context

- Layout lab: worktree `MDV-ui-refresh`, `src/react/prototypes/layoutLab`, opened with `?layoutLab=grid` or `?layoutLab=tiles`. Measurement script: `artifacts/playwright-cli/layout-lab/measure.mjs`.
- The gaps each engine needs filled: gridstack 14 has no keyboard moves and no documented save of narrow layouts; the in-house engine lacks touch, resizing from the left and top, screen-reader support, chart hosting and tests.
- "Unsaved Grid Stack View" is MDV's integration, not gridstack. A gridstack run must avoid the three causes (early baseline, pixel sizes saved, pushed charts never written back) or it isn't a fair test.
- Settled before this map: a scrolling grid, not tiles.
- From [18 One view, saved layout and screen sizes](18-one-view-saved-layout-and-screen-sizes.md) (ADR 0026): the engine maps a saved 12-column `{x, y, w, h}` layout and derives 6- and 1-column layouts in reading order; nothing narrower is stored. Moving and resizing happen only at 960 px and up, so touch dragging matters less than it did when this ticket was written.
