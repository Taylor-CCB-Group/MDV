# Use the prototype's engine as the showcase baseline behind a replaceable layout seam

**Status:** accepted as the showcase baseline, not as the lasting layout-model decision.
Follows ADR 0026's candidate internal model. The team will compare scrolling grid, tiling and
hybrid arrangements with the walking skeleton in use.

## Context

The one-view workspace (ADR 0026) needs a layout engine. It places each chart on a 12-column
grid, lets people drag a chart by its title bar and resize it from its edges, and lays out the
6- and 1-column versions for narrower screens from the one saved 12-column layout.

The workspace prototype already has one: `compact`, `flow` and `findSpot`, about 65 lines of
layout functions, plus its pointer handling. The candidates were that engine, gridstack 14.1
(MDV runs gridstack 7.3 in the pane workspace today) and react-grid-layout 2.3.

A layout lab hosted the real charts of three projects in all three engines: a spatial view, a
23-chart view, and a 16-chart view over three data sources. All three used the same host
component, the same converted layout, the prototype's look and gesture rules, and the
prototype's `flow` for narrower layouts. Then each was checked automatically at 1440, 1024,
768 and 390 px and tried by hand.

| | In-house (prototype) | gridstack 14 | react-grid-layout 2.3 |
| --- | --- | --- | --- |
| Overlaps or charts past the edge | none | none | none |
| Save, reload, compare | exact | exact | exact |
| Drag, worst frame | 9–11 ms | 9–16 ms | 32–38 ms |
| Resize, worst frame | 9–24 ms | 9–34 ms | 40–83 ms |
| Growing a chart past the bottom of the window | yes: spare rows while dragging, timer scroll | no: scrolls only on pointer moves and only to the current end of the grid | no, same |
| Our own code in the lab | ~400 lines | ~180 lines | ~110 lines |
| Working against the library | n/a | 10 `!important` overrides for its handles, its column-layout callback, React and gridstack both moving the same nodes | none found |

By hand, the in-house engine and gridstack 14 felt almost the same. The one difference was
the bottom edge.

Findings: `docs/design/workspace/LAYOUT-EVALUATION.md`.

The worry with an in-house engine is owning it. The layout logic is small and testable. What we
would own is the gesture code, and gridstack 14 already needed its gestures worked around in
three places to match the prototype. MDV's problems with gridstack today come from how
`GridstackManager` wires it in, not from gridstack itself: row height taken from the window,
one grid per data source, positions saved by watching inline styles. Either engine replaces
that wiring.

## Decision

The first showcase implementation uses the prototype's engine, rebuilt as a typed module with
tests. It sits behind a seam narrow enough that gridstack 14, a tiling engine or a hybrid can
replace it without touching the rest of the workspace. This choice starts implementation; it
does not prejudge the product-level layout model.

### Behaviour

- 12 columns from 960 px, 6 from 600 px, 1 below (ADR 0026).
- At 12 columns, a chart is dragged by its title bar and resized from its edges. Other charts
  move out of the way while it is dragged and are pulled up into free space after it is dropped
  (`compact`).
- At 6 and 1 columns, the layout is the 12-column one in reading order, with widths scaled and
  rows wrapped (`flow`). Nothing narrower is saved, and charts cannot be moved or resized there
  (ADR 0026).
- New and duplicated charts go in the first free place (`findSpot`).
- While dragging, the grid has spare rows below the last chart and scrolls on a timer near the
  top and bottom edges. Esc and `pointercancel` cancel the gesture.
- A chart is redrawn once when a gesture ends, not on every frame.
- Alt+arrows move the focused chart one cell and Shift+Alt+arrows resize it.

### The engine seam

The engine is one module with one job: show the charts in a given 12-column layout and report
the new layout when a gesture ends. Everything else stays outside it, so another engine can be
swapped in.

- **In:** the charts (id and title), the saved 12-column `{x, y, w, h}` layout, and whether
  editing is allowed. The engine works out the column count from its width, using the shared
  breakpoints.
- **Out:** `onLayoutChange(layout)` with a whole 12-column layout, once per finished gesture.
  The engine also calls `gestureStart` and `gestureEnd` so hosts can hold back redraws.
- **Owned outside the engine, and shared by any engine:**
  - the saved format (ADR 0026)
  - the breakpoints and `flow` for narrower layouts
  - `compact` and `findSpot` as plain functions
  - keyboard moves, which take a layout and return a new layout, handed back in as the layout
  - the chart host
  - redraw timing
  - saving
- **Nothing engine-specific leaks out.** No engine class names, node objects or options appear
  outside the module. Tests and other code find cards by `data-chart-id` and roles, so the
  acceptance checks from the walking skeleton judge a replacement engine unchanged.

The layout lab already runs all three engines behind this shape (`charts`, `layout`,
`onLayoutChange`), so the seam is known to fit gridstack 14.

### What hosting real charts needs, whichever engine

The lab found these, and the build has to do each one:

- **Move the chart's own div into the card** and call `setSize` when the card's size settles
  (debounced, and held back during a gesture), as the lab's `ChartHost` does.
- **Remove the pane workspace's drag and resize from hosted charts** (`removeDraggable`,
  `removeResizable`). Left on, a drag moves the chart inside its card, and the `resizer-*`
  handles stick out past the card and catch the pointer meant for a neighbour.
- **Contain each chart's stacking** (`isolation: isolate` on the host) so chart overlays cannot
  cover the card's resize handles.
- **Mount the grid only once the view's layout is loaded.** Mounting with the previous view's
  layout placed every chart at a fallback size.
- **Do not redraw a chart while its card is animating to a new size.** On load, the 16-chart
  view threw ~350 "negative SVG size" errors under any engine that animates width and height.

### Still to build beyond the lab

- Resize handles of at least 24 px (WCAG 2.2 target size), and pointer capture.
- Screen-reader announcements for keyboard moves, and DOM order that follows the layout.
- A long-press before a touch drag, so a drag does not fight page scrolling on touch laptops.
- A guard in `findSpot` for a chart wider than the grid, which currently loops forever.

The layout functions get unit tests. Dragging and resizing get Playwright tests on the
acceptance views (ticket 20).

## Layout checkpoint

Before hardening layout-specific UI, compare the baseline with a small tiling or hybrid
experiment. Shared Playwright scenarios provide geometry assertions, measurements and
screenshots; people try the same workflows and judge feel. Test:

- an adjustable 2/3 arrangement
- moving selected charts between regions
- dense small charts without disruptive movement of unrelated charts
- minimum chart dimensions, with no degeneration into thin slices
- avoiding unnecessarily taller-than-viewport charts
- desktop and narrower laptop/tablet widths
- existing spatial/WebGL chart hosting

The team then records grid, tiling or a hybrid as the durable model. Switching to gridstack 14
remains an implementation option if the grid model wins. It means writing a gridstack 14
version of the engine module behind the same seam. The
lab's version is about 180 lines. It also has to carry these known workarounds:
- resize handles moved inside the card and kept visible
- the column-layout callback fed from `flow`
- spare rows and a timer scroll to grow a chart past the bottom of the window, which still
  stalls if the pointer does not move
- keeping React and gridstack from both owning the card nodes

gridstack 14 and the pane workspace's gridstack 7.3 would sit side by side under a package
alias until the pane workspace is retired.

## Considered and rejected for the initial baseline

- **gridstack 14 now.** It is as smooth and the least risky to own for touch, but it cannot
  grow a chart past the bottom of the window. Matching the prototype needed overrides of its
  handles, its column layouts and its scrolling. Its CSS overrides depend on class names and
  inline styles that have changed across its major versions. It is kept as the fallback above.
- **react-grid-layout 2.3.** The least of our own code, and nothing fought it. But it was the
  least smooth, with worst frames of 32–83 ms while dragging and resizing. It also has no
  growing past the bottom of the window and no keyboard moves.
- **Tiles or a hybrid.** Not selected as the initial baseline because the existing lab already
  proves the in-house grid with real charts. They remain open at the layout checkpoint;
  `react-mosaic` is one implementation to investigate, not a foregone choice.

## Consequences

- For the showcase baseline, the grid is MDV's code to maintain: about 65 lines of layout
  functions and a gesture component.
- The engine seam is a rule for every later slice. Anything that reaches into the engine's DOM
  or state makes the fallback more expensive.
- gridstack 7.3 and `GridstackManager` stay only for the pane workspace, and go when it is
  retired.
- WebGL contexts are lost when a view has more WebGL charts than the browser allows. This
  happens under every engine, and it is a chart fix tracked with the walking skeleton.
