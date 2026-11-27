Ticket: [19 Layout engine bake-off and choice](../tickets/19-layout-engine-bake-off-and-choice.md)
Map: [Workspace in MDV - Map](../README.md)
Source: layout lab in the `MDV-ui-refresh` worktree (2026-11-26), `src/react/prototypes/layoutLab`, opened with `?layoutLab=inhouse` or `?layoutLab=gridstack`. Script: `artifacts/playwright-cli/layout-lab/bakeoff.mjs`, results and screenshots in `artifacts/playwright-cli/layout-lab/bakeoff/`. Headless Chromium on the Metal GPU, 1440×900, dev server data.

## Set-up

- Both engines host the real chart divs through the same `ChartHost` and start from the same converted layout (ADR 0026: pane by pane, reading order, widths relative to the pane, packed into 12 columns). Both save the v2 shape (`{version: 2, layout}`; localStorage stands in for `views_v2.json`).
- 12 columns from 960 px, 6 from 600 px, 1 below. Moving and resizing only at 12 columns.
- gridstack 14.1.0 is installed as `gridstack14`, so MDV's gridstack 7.3 is untouched. Its narrower layouts use its `list` mode. The in-house engine is the prototype's `compact` and `flow` with pointer handling, about 250 lines in the lab.
- Projects: a spatial view (10 charts including a SpatialData chart), a 23-chart view, and a 16-chart view over three data sources mixing absolute and gridstack panes. Project names are kept out of this public copy.

## Results

| | In-house | gridstack 14 |
|---|---|---|
| Overlaps or charts past the edge at 1440, 1024, 768, 390 px | none, all three projects | none, all three projects |
| Tablet (768 px, 6 columns) | two charts per row, widths halved; page as tall as desktop (23-chart view 9,810 px; 16-chart 6,388 px) | `list` mode makes every chart full width, so tablet looks like phone (17,342 px; 9,860 px). A custom column-layout callback could run our own `flow` instead. |
| Phone (390 px, 1 column) | one column, reading order | same |
| Dragging a chart | 120 fps, worst frame 10–28 ms, one 308 ms frame once on the 23-chart view | 119–120 fps, worst frame 10–17 ms |
| Chart redraws while dragging | none (redrawn once after the drop) | none |
| Resizing a chart | 120 fps; charts redrawn 4–30 times per resize, on each pause | not measured: the script could not grab gridstack's handle (see below) |
| Save, reload, compare | exact on all three; no drift after charts settle | exact on all three; no drift after charts settle |
| Spatial geometry reloads during drags (`perf.ts`) | none | none |

## What broke, and where the fault was

- **Lab, not the engines:** after a view change the engine mounted before the new view's layout was loaded. gridstack, which reads its layout once, then placed every chart at a fallback size; the in-house engine reads it live and hid the same bug. Fixed by not mounting until the view's layout is loaded. A real integration has the same ordering to get right.
- **Both engines:** the pane workspace's own resize handles (`resizer-*` from `Elements.js` `makeResizable`) stay on each chart div, stick out past its edges and catch the pointer meant for the neighbouring chart's handle. Hosting has to remove them.
- **Both engines:** chart content with its own stacking (deck canvases, SVG overlays) sat over the resize handles until the host was isolated (`isolation: isolate`) and the handles raised.
- **gridstack:** resize handles are hidden until the pointer moves over the item, and even with the two fixes above the script's pointer landed on the chart panel, not the handle. Needs trying by hand; not shown to be a gridstack bug.
- **Charts, with either engine:** WebGL contexts lost on the 23-chart view (more WebGL charts than the browser allows), and hundreds of "negative SVG size" errors when charts were drawn while very small during transitions. These appeared under both engines in different steps, so they are chart robustness issues to fix whichever engine is chosen.
- **Script gaps:** one in-house drag on the 16-chart view did not start (the start point did not land on the title bar), and gridstack resizing was never measured. Both need checking by hand.

## Not measured

- Touch, keyboard and screen readers (neither tried).
- Resizing with gridstack.
- A bigger sample: each gesture ran once per project.
