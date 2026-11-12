# Scatterplot benchmark harness (sketch)

Context: issue #433 compares the classic (`wgl_scatter_plot`) and deck
(`DeckContourScatter`) scatterplots. Filtering now runs in workers, so the
numbers in that issue may be stale. This sketches a repeatable comparison; it
is not implemented.

## Shape

A Playwright spec under `tests_playwright/perf/` (kept out of the default
project run), using `createTemporaryProjectViaSyntheticAnndata` with large
`nCells`. Each scenario runs once per chart type, with the other chart absent,
and reports median and p95 over N repeats to a JSON file in
`artifacts/` (gitignored), so runs can be diffed. Timings are taken inside the
page with `performance.now()` and `requestAnimationFrame`, not from Playwright
round trips.

Parameters: `nCells` in {100k, 1M, 5M}; chart type in `wgl_scatter_plot` and
`DeckContourScatter`; render backend WebGL vs WebGPU for deck (the Browser
setting, via `localStorage` before load). Always headed Chromium with GPU
enabled; headless software rendering makes the frame numbers meaningless.

## Scenarios

1. **Filter change (click latency).** Add a scatter plus a row chart on a
   categorical column. Programmatically set the category filter, then time from
   the call until the scatter has painted the result. Alternate between two
   categories so every repeat does real work.
2. **Pan/zoom with a filter active.** Apply a filter that keeps about 50%, then
   drive a fixed synthetic pan/zoom path (scripted wheel/drag events) for 3s and
   record frame intervals. Same path with no filter as the baseline, since the
   issue says classic collapses only when filtered.
3. **Zoomed-out overdraw.** No filter, view fitted to all data, 3s of idle plus
   a small zoom wiggle; record frame intervals. Targets the deck overdraw
   complaint.
4. **Selection drag.** Draw a rectangle selection over a fixed path. Classic
   commits on mouse-up, deck updates continuously; record per-move update cost
   and total time to final filtered count.

## Measuring "painted"

The two charts need different completion signals, so wrap them behind one
helper the spec calls:

- Common: a `dataStore` listener fires on filter change; after that wait for two
  `requestAnimationFrame`s so the canvas has had a chance to present.
- Classic: `WGLScatterPlot.onDataFiltered` -> `drawChart`. Wrap `drawChart` in
  the page (monkey-patch on the chart instance, same style as the existing
  `_addChart` patch in `chart_creation_single.spec.ts`) to stamp the end time.
- Deck: there is no draw callback on the React component. Options in order of
  preference: (a) a dev-only `performance.mark` in the filtered-indices effect
  in `DeckScatterComponent` (`useFilteredIndices` refresh) and in deck's
  `onAfterRender` for the layer; (b) deck's `_onMetrics` /
  `deck.metrics`. Pick (a) unless it proves too invasive.

Frame intervals for scenarios 2-3 come from a `requestAnimationFrame` loop
recording deltas; report mean, p95 and the share of frames over 33ms. This works
for both charts without chart hooks.

## Open questions

- Should this live as a plain script against a running dev server instead of a
  Playwright spec? Playwright gives fixtures and cleanup for free but is
  heavier to run repeatedly.
- The marks in `DeckScatterComponent` would be shared infrastructure for the
  WebGPU filter engine work, so it is worth deciding on a labelled mark API
  (like `src/react/spatialdata/perf.ts`) before adding any.
- Scripted wheel events do not reproduce Mac trackpad inertia, which is the
  classic chart's reported UX problem; that stays a manual check.
