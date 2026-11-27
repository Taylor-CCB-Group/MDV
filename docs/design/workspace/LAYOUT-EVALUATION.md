# Layout evaluation

This document consolidates both layout-lab rounds and defines the next comparison. The lab
lives in `src/react/prototypes/layoutLab/` and opens in development builds with
`?layoutLab=inhouse|gridstack|rgl`.

## Test set

The 2026-11-26 bake-off used three private development projects whose names must not be
committed:

- a spatial view with 10 charts, including a SpatialData chart
- a dense 23-chart view
- a 16-chart view over three data sources mixing old absolute and gridstack panes

The scripts ran headless Chromium on a Metal GPU at 1440, 1024, 768 and 390 px. Every engine
used the same `ChartHost`, converted layout and responsive `flow`.

## Round one

| Check | In-house | gridstack 14 |
|---|---|---|
| Overlap or edge overflow | none | none |
| Tablet | two charts per row | list mode made every chart full width |
| Drag | about 120 fps; worst 10–28 ms | 119–120 fps; worst 10–17 ms |
| Redraw while dragging | none; once after drop | none |
| Resize | about 120 fps; repeated redraw on pauses | not measured; script missed handle |
| Save/reload | exact | exact |
| Spatial geometry reload during drag | none | none |

The first round exposed lab and hosting faults rather than a decisive engine difference:

- the engine mounted before the new view layout was ready
- legacy pane resize handles caught pointers
- chart stacking covered layout handles
- the script did not reliably hit all title bars or gridstack handles
- WebGL limits and negative SVG sizes occurred independently of the layout engine

## Round two

Round two removed legacy drag/resize, shared the prototype styling, delayed redraw until gesture
end and fixed each candidate's responsive wiring.

| Check | In-house | gridstack 14 | react-grid-layout 2.3 |
|---|---|---|---|
| Overlap or edge overflow | none | none | none |
| Drag worst frame | 9–11 ms | 9–16 ms | 32–38 ms |
| Resize worst frame | 9–24 ms | 9–34 ms | 40–83 ms |
| Redraws per gesture | 0 move, 1 resize | 0–1 | 0–1 |
| Grow page while dragging down | yes | no; current end only | no; current end only |
| Save/reload | exact | exact | exact |
| Spatial geometry reload | never | never | never |
| MDV-owned lab code | about 400 lines | about 180 lines | about 110 lines |

In-house and gridstack felt similar by hand. Gridstack needed handle-placement/visibility
overrides and a custom column callback. React-grid-layout required the least adapter code but
had the worst long frames. None of this settled scrolling grid versus tiling.

## Hosting rules that survive any layout choice

- Move the existing chart div into the card.
- Remove legacy pane drag and resize from hosted charts.
- Put `isolation: isolate` on the host so chart overlays cannot cover layout handles.
- Mount only after the current view layout is ready.
- Do not redraw while card dimensions animate; redraw once after a gesture settles.
- Keep layout engine DOM, node objects and class names behind the seam.
- Find cards in tests with roles and `data-chart-id`.
- Spatial geometry must not reload for a move or cosmetic resize. Measure through
  `src/react/spatialdata/perf.ts` and switch the flag off afterwards.

## Next checkpoint: grid, tiling or hybrid

The in-house grid remains the initial showcase baseline because it already hosts real charts.
Before hardening layout-specific UI, time-box a minimal tiling or hybrid candidate behind the
same test-facing seam.

### Shared Playwright scenarios

Read `docs/playwright/PLAYWRIGHT_AGENT.md` before modifying or running the checks. Run the same
scenarios against each candidate using roles and `data-chart-id`, not implementation classes:

1. Resize an adjustable 2/3 arrangement.
2. Move selected charts between its regions.
3. Exercise the dense view and measure displacement of charts unrelated to the gesture.
4. Enforce agreed minimum chart width/height; no chart may degenerate into a thin slice.
5. Avoid individual charts taller than the viewport unless a test explicitly calls for it.
6. Check desktop and narrower laptop/tablet widths.
7. Host spatial/WebGL charts without geometry reloads.
8. Assert no overlap, edge overflow or chart content outside its card.
9. Capture comparable screenshots and geometry measurements.

Playwright supplies repeatability and objective geometry. It cannot decide which interaction
feels better. Team members must try the same workflows by hand.

### Decision

Record one of:

- scrolling grid
- tiling
- hybrid, such as resizable regions containing grids

The decision must explain the tested workflows, human preference, implementation cost,
accessibility implications and responsive behaviour. Record it in a new ADR; ADR 0027 remains
the history of why the showcase started with the in-house grid.
