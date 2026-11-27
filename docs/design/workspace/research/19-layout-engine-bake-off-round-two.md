Ticket: [19 Layout engine bake-off and choice](../tickets/19-layout-engine-bake-off-and-choice.md)
Map: [Workspace in MDV - Map](../README.md)
Round one: [19 Layout engine bake-off - results](19-layout-engine-bake-off-results.md)
Source: layout lab in the `MDV-ui-refresh` worktree (2026-11-26), `?layoutLab=inhouse|gridstack|rgl`. Script `artifacts/playwright-cli/layout-lab/bakeoff.mjs`; results `bakeoff/results.json` (in-house, react-grid-layout) and `bakeoff/results-gridstack.json` (gridstack, re-run after fixing the lab's reflow callback). Headless Chromium on the Metal GPU, dev server data.

## Why a second round

The human tried round one by hand and found it glitchy and nothing like the prototype: a chart slid out of its card when dragged, placing a chart was fiddly, the page did not grow or scroll when dragging down, and resizing at the bottom of the window flickered. Round two fixes the lab first so only the engine differs.

## What changed in the lab

- **The pane workspace's own drag and resize are removed from each hosted chart** (`removeDraggable`, `removeResizable`). They had moved the chart inside its card (the black card) and their handles caught the pointer.
- **One look for all three engines:** the prototype's canvas, cards, header, lift shadow, dashed placeholder and handle placement. The lab controls moved to a top bar, off the grid.
- **Charts redraw once, after a gesture,** for every engine (held back while dragging or resizing, as in the prototype).
- **In-house engine:** the prototype's gesture handling ported properly: the card follows the pointer, the others slide, the grid grows while dragging and auto-scrolls near the top and bottom (time-based, ~600 px/s), Esc cancels.
- **gridstack 14:** narrower layouts come from the prototype's `flow`, passed as its column-layout callback (its `list` mode made every chart full width at 6 columns). Its resize handles had to be moved inside the card (by default they sit in the gap below, under the next chart) and kept visible (by default they appear only after a hover).
- **react-grid-layout 2.3:** added back, with the prototype's `flow` handed in as its 6- and 1-column layouts.

The four reported glitches were checked after the changes, on all three engines: no chart leaves its card; holding a dragged chart at the bottom scrolls (in-house also grows the page); resizing at the bottom grows steadily with no flicker.

## Results (three projects, each at 1440, 1024, 768, 390 px)

| | In-house | gridstack 14 | react-grid-layout 2.3 |
|---|---|---|---|
| Overlaps or charts past the edge | none | none | none |
| Tablet (6 columns) | two charts per row (prototype's flow) | same | same, then compacted upwards (a little shorter) |
| Drag | 120–130 fps, worst frame 9–11 ms | 120–128 fps, worst 9–16 ms | 126–157 fps, worst 32–38 ms |
| Resize | 120–132 fps, worst 9–24 ms | 120–126 fps, worst 9–34 ms | 118–143 fps, worst 40–83 ms, one long frame on two projects |
| Chart redraws per gesture | 0 for a move, 1 for a resize | 0–1 | 0–1 |
| Page grows while dragging down | yes | no, scrolls to the current end | no, scrolls to the current end |
| Save, reload, compare | exact, no drift | exact, no drift | exact, no drift |
| Spatial geometry reloaded during gestures | never | never | never |
| Our own code (lab) | ~400 lines: 300 component, 65 engine core, 54 CSS | ~180 lines: 136 component, `flow`, 24 CSS of which 10 override its handles | ~110 lines: 76 component, `flow`, 15 CSS |
| Fighting the library | n/a | handle placement and auto-hide; the column callback's contract (items to place are the 4th argument, pushed into the 3rd) | none found |

All three share the ~55-line `ChartHost` that fits a legacy chart into a card.

## Also found

- **Charts drawn while their card is animating:** on the 16-chart view, loading or reconverting threw ~350–370 "negative SVG size" errors under in-house and gridstack, and none under react-grid-layout, which animates position only. The host should redraw a chart once its card has stopped changing size, not during the transition. Applies to whichever engine is chosen.
- WebGL contexts lost on the 23-chart view under every engine (more WebGL charts than the browser allows), as in round one.
- The script could not start the drag on one react-grid-layout run (16-chart view); everything else was measured.

## Not measured

Touch, screen readers, keyboard moves in gridstack and react-grid-layout (neither has them built in), and how each feels by hand, which is the human's part.
