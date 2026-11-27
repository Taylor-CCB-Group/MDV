Type: prototype
Status: resolved
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

- Layout lab: worktree `MDV-ui-refresh`, `src/react/prototypes/layoutLab`, opened with `?layoutLab=inhouse`, `gridstack` or `rgl` (round two). Measurement script: `artifacts/playwright-cli/layout-lab/measure.mjs`.
- The gaps each engine needs filled: gridstack 14 has no keyboard moves and no documented save of narrow layouts; the in-house engine lacks touch, resizing from the left and top, screen-reader support, chart hosting and tests.
- "Unsaved Grid Stack View" is MDV's integration, not gridstack. A gridstack run must avoid the three causes (early baseline, pixel sizes saved, pushed charts never written back) or it isn't a fair test.
- Settled before this map: a scrolling grid, not tiles.
- From [18 One view, saved layout and screen sizes](18-one-view-saved-layout-and-screen-sizes.md) (ADR 0026): the engine maps a saved 12-column `{x, y, w, h}` layout and derives 6- and 1-column layouts in reading order; nothing narrower is stored. Moving and resizing happen only at 960 px and up, so touch dragging matters less than it did when this ticket was written.

## Bake-off results

Findings: [19 Layout engine bake-off - results](../research/19-layout-engine-bake-off-results.md). Both engines held every layout rule at all four widths and saved and reloaded exactly. They differ at tablet width (gridstack's `list` mode makes every chart full width) and in what we would own versus depend on.

Round two ([19 Layout engine bake-off - round two](../research/19-layout-engine-bake-off-round-two.md)) fixed the lab after the human found it glitchy and unlike the prototype, put all three engines (in-house, gridstack 14, react-grid-layout) inside the prototype's look, and gave each the prototype's reflow. All three pass every check; they differ in smoothness, what the library fights, and how much code is ours. By hand, in-house and gridstack 14 felt almost the same; only in-house could grow a chart past the bottom of the window.

## Answer

**The in-house engine (the prototype's), for now, behind a seam that lets gridstack 14 replace it if the team prefers it.** Recorded as [ADR 0027](../../../adr/0027-build-the-workspace-grid-on-the-prototypes-own-engine.md).

- **Why in-house:** as smooth as gridstack (worst frames 9–24 ms), the only one that grows a chart past the bottom of the window, and the prototype the team liked is its spec. The layout logic is ~65 lines of plain, testable functions; what we own is the gesture code, and gridstack already needed working around in three places (handles, column layouts, bottom-edge scrolling) to match the prototype.
- **MDV's gridstack pain is the wiring, not gridstack:** `GridstackManager` takes row height from the window, makes one grid per data source and saves by watching inline styles. Nothing in MDV patches gridstack itself. Either engine replaces that wiring.
- **The seam (the fallback's price):** the engine takes the charts, the saved 12-column layout and whether editing is allowed, and reports a whole 12-column layout once per finished gesture, plus gesture start and end. The following stay outside the engine:
  - the saved format and the breakpoints
  - `flow`, `compact` and `findSpot`
  - keyboard moves
  - the chart host
  - redraw timing
  - saving

  No engine class names leak out, and tests find cards by `data-chart-id`. The lab already runs all three engines behind this shape.
- **Switch to gridstack 14 if:**
  - the team prefers it after using the skeleton
  - tablets must move or resize charts
  - we want drag-in from a palette or nested grids
  - our gesture bugs cost more than gridstack's workarounds

  The lab's gridstack module (~180 lines) and its known workarounds are listed in the ADR.
- **Gives up:** built-in touch dragging (we add a long-press) and an outside maintainer for the gesture code.
- **Every engine needs:**
  - the chart div moved into the card
  - the pane workspace's own drag and resize removed (`removeDraggable`, `removeResizable`)
  - `isolation: isolate` on the host
  - the grid mounted only after the view's layout loads
  - no chart redraws while its card animates (~350 negative-SVG-size errors otherwise)
  - one redraw per gesture
- **Still to build:** 24 px handles and pointer capture, screen-reader announcements, a touch long-press, and a `findSpot` guard for over-wide charts.

