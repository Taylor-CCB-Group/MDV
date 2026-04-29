# Synthetic SpatialData Projects

This note records the first-pass convention for using generated SpatialData-like
projects in MDV development, agent workflows, and performance investigations.

## Why

We want representative, repeatable projects that can be generated at arbitrary
sizes, opened interactively from the normal project catalog, and used by agents
or tests without relying on large checked-in fixtures.

The immediate package candidate is `dummy-spatialdata`, added to the Python dev
dependency group for Python 3.12 environments. It depends on the SpatialData
ecosystem and `dummy-anndata`. A later pass should review whether
`dummy-anndata` can replace part of `python/mdvtools/tests/mock_anndata.py`.

## First-Pass Scope

Keep this deliberately modest:

- generate spatial projects with configurable cell counts;
- include enough categorical, numeric, text, coordinate, image/region, and table
  columns to exercise common charts;
- create a small set of named views for performance comparison;
- make projects easy to register in the normal catalog;
- make cleanup safe and explicit.

Do not try to cover every chart configuration permutation yet. The first goal is
repeatable comparison and developer confidence, not exhaustive combinatorics.

## Proposed Project Location

Use immediate children of `~/mdv` for interactive projects because the existing
filesystem scanner treats each direct child as a candidate MDV project. Do not
nest generated projects under grouping folders inside `~/mdv`.

Recommended naming:

```text
~/mdv/synth-spatial--<profile>--<run-id>/
```

Examples:

```text
~/mdv/synth-spatial--scatter-table--10k/
~/mdv/synth-spatial--scatter-table--100k/
~/mdv/synth-spatial--scatter-table--1m/
~/mdv/synth-spatial--selection-dialog--100k/
```

The generated MDV project itself should live at the final path. Intermediate
SpatialData stores can live under the project, for example:

```text
~/mdv/synth-spatial--scatter-table--100k/spatial/source.zarr/
```

That keeps the project portable and lets cleanup remove one subtree.

## Registration

For a running database-backed app:

1. Generate the project into a direct child path such as
   `~/mdv/synth-spatial--scatter-table--100k`.
2. Use the existing `/rescan_projects` route, or restart the app, so
   `serve_projects_from_filesystem()` registers filesystem projects missing from
   the database.
3. Generated project names should start with `synth-spatial:` so they are easy
   to find and cleanup by name pattern.

For direct/manual serving, use the project path without DB registration.

Project `state.json` should include provenance:

```json
{
  "provenance": {
    "created_by": "mdvtools synthetic spatialdata generator",
    "generator": "dummy-spatialdata",
    "profile": "scatter-table",
    "n_cells": 100000,
    "seed": 42,
    "cleanup_group": "synth-spatial"
  }
}
```

This matches the existing cleanup script's idea of test-generated projects while
giving us a more specific selector for future cleanup commands.

## Cleanup

Generated projects should be disposable by convention.

Filesystem cleanup:

```bash
rm -rf ~/mdv/synth-spatial--<profile>--<run-id>
```

Database-backed cleanup should prefer existing project-management paths:

- soft-delete via `/delete_project/<id>` when operating through the app;
- use `python/mdvtools/dbutils/cleanup_projects.py` for local development
  database cleanup;
- future improvement: add a selector that filters `state.json.provenance.cleanup_group == "synth-spatial"`.

Avoid deleting all of `~/mdv` from an agent. Only delete direct child project
paths that match the `~/mdv/synth-spatial--...` naming convention.

## Candidate Generator Interface

Use the module command rather than ad hoc notebooks:

```bash
cd python
../venv/bin/python -m mdvtools.tests.generate_synthetic_spatial_project \
  --profile scatter-table \
  --n-cells 100000 \
  --seed 42 \
  --output ~/mdv/synth-spatial--scatter-table--100k \
  --force
```

Current options:

- `--profile`: currently `scatter-table` or `spatial-overview`;
- `--n-cells`: primary scale parameter;
- `--n-genes`: generated expression-like variable count;
- `--image-size`: generated image extent in pixels;
- `--seed`: fixed by default for reproducibility;
- `--output`: output MDV project path, defaulting to a flat
  `~/mdv/synth-spatial--<profile>--<n-cells>` path;
- `--force`: replace an existing generated project at the exact output path.

The first implementation creates a representative SpatialData object with an RGB
image, circle shapes, an AnnData table linked to the shapes, and
`obsm["spatial"]` row coordinates. The explicit `obsm["spatial"]` matrix is
important: MDV's current converter treats tables without it as non-spatial even
if SpatialData table metadata links them to an annotated element.

For larger `scatter-table` samples, the generator switches to an image-annotated
table instead of one shape per row. This keeps the project SpatialData-backed
while avoiding million-feature GeoJSON/shapely overhead that would obscure the
scatter/table performance question. The 1m-row sample generated in this pass is:

```text
~/mdv/synth-spatial--scatter-table--1m/
```

A 10m-row stress sample was also generated with only two genes to keep the
stress centred on row count, scatter plots, and tables:

```text
~/mdv/synth-spatial--scatter-table--10m/
```

Planned but not implemented yet:

- `--n-regions`;
- `--cleanup`;
- automatic DB registration;
- profile-specific performance comparison views.

## Initial Profiles

### `scatter-table`

Purpose: compare core scatter and table implementations at different row counts.
The generator currently creates the backing data; the comparison views below are
the next layer to add.

Views:

- `Deck scatter + React table`
- `WGL scatter + legacy table`
- `Deck scatter + legacy table`
- `WGL scatter + React table`

Charts/components of interest:

- `src/react/components/DeckScatterReactWrapper.tsx`
- `src/charts/WGLScatterPlot.js`
- `src/charts/TableChart.js`
- `src/react/components/TableChartReactWrapper.tsx`

### `selection-dialog`

Purpose: measure the cost of adding selection/filter UI and linked selection
state.

Views:

- `Scatter/table baseline`
- `Scatter/table + selection dialog`

Component of interest:

- `src/react/components/SelectionDialogComponent.tsx`

### `spatial-overview`

Purpose: keep one representative spatial view with regions/images/coordinates so
conversion regressions are visible.

Views:

- `Spatial overview`
- `Spatial overview + table`

## Performance Measurements

Start with browser-level metrics that agents can collect consistently:

- time to first visible project content;
- time until chart manager and charts are ready;
- row count and chart count;
- interaction latency for pan/zoom/filter/select where practical;
- browser console errors;
- optional trace and screenshot artifacts under `output/playwright/`.

Use Playwright from the repo root:

```bash
TEST_BASE_URL=http://127.0.0.1:5174 npm run playwright-test -- tests_playwright/project/ --project=chromium --reporter=list
```

For focused performance work, add dedicated specs rather than overloading the
functional project tests.

A lightweight Playwright profiler is available for ad hoc row-chart filter
timing against an already-served project:

```bash
node scripts/profile_mdv_interactions.mjs \
  --url http://localhost:5055/project/191 \
  --views "react view,classic view" \
  --categories type_0,type_1,type_2 \
  --out output/playwright/project-191-row-filter-profile.json
```

The profiler records:

- time from category click to the datastore `filtered` event;
- time to materialise `DataStore.getFilteredIndices()`;
- time after two animation frames, as a coarse visual-settle marker;
- chart/view metadata and console warnings/errors.

First 1m-row observation: category filtering itself was tens of milliseconds,
while React chart updates were dominated by filtered-index materialisation. The
first cold React click took much longer than subsequent clicks, so repeat runs
and warm/cold separation should be part of the next profiling pass.

## UI-Authored Comparison Views

The Python chart helper API is useful for older charts, but it should not be the
source of truth for newer React-backed chart configs. For comparison views, use
the normal UI path first and inspect the saved `views.json`.

Current probe result from `~/mdv/synth-spatial--ui-chart-probe--25`:

- UI label `Row Chart` saved as `row_chart`;
- UI label `2D Scatter Plot` saved as `wgl_scatter_plot_dev`;
- UI label `2D Scatter Plot (Classic)` saved as `wgl_scatter_plot`;
- UI label `Table` saved as `table_chart_react`;
- UI label `Table (Classic)` saved as `table_chart`.

The new scatter is implemented by `DeckScatterReactWrapper.tsx`, but the
persistent chart type is currently `wgl_scatter_plot_dev`, not `DeckScatter`.
That makes UI-authored config a safer starting point than hand-written Python
chart config for these comparison views.

Recommended agent/human loop:

1. Generate or copy a flat project under `~/mdv`, for example
   `~/mdv/synth-spatial--ui-chart-probe--25`.
2. Serve it editably:

   ```bash
   ./venv/bin/python -c 'from mdvtools.mdvproject import MDVProject; p=MDVProject("/Users/ptodd/mdv/synth-spatial--ui-chart-probe--25"); p.set_editable(True); p.serve(open_browser=False, port=5062, websocket=False)'
   ```

3. Use Playwright against `http://127.0.0.1:5062/` to add charts through the Add
   Chart dialog.
4. Save the current view through the app, then inspect
   `~/mdv/synth-spatial--ui-chart-probe--25/views.json`.
5. If the saved config is stable and useful, promote that view shape into a
   reusable fixture or generator step.

One caveat from the first probe: the Add Chart defaults chose
`quality_score x quality_score` for both scatter implementations and only
`cell_type` for the table implementations. That is good enough to establish the
saved config shape, but performance comparison views should explicitly set
`x/y` for scatter plots and include a broader table column list.

The current generated comparison view uses `x/y` for both scatter charts and all
non-internal cell columns for both table charts.
Generated `scatter-table` projects also include split `react view` and
`classic view` views so Playwright profiling can compare one chart family at a
time.

## Interview Questions

These are the decisions to settle before writing the first generator:

1. What row-count ladder should represent "small", "normal", "large", and
   "stress" for day-to-day work? A starting guess is `10k`, `100k`, `1m`.
2. Do we care more about initial render time, steady-state interaction latency,
   memory growth, or all three for the first benchmark?
3. Should generated projects be registered automatically when a DB-backed app is
   running, or should generation only print the rescan/open instructions?
4. Should test-generated catalog entries be visible to everyone in local auth
   mode, or owned by the current admin user after rescan?
5. Which charts are mandatory in the first comparison view, and which should be
   follow-up once the harness is stable?
6. What pass/fail threshold is useful initially: absolute budgets, regression
   against a stored baseline, or just structured measurement artifacts?
7. Should generated SpatialData source stores be retained inside the MDV project
   for debugging, or optionally omitted after conversion to save disk?
8. What naming convention do you want in the catalog: terse machine names
   (`synth-spatial:scatter-table:100k`) or more readable names
   (`Synthetic spatial scatter/table 100k`)?

## Subsequent Passes

1. Inspect `dummy-spatialdata`'s API directly after installing it and map its
   generated elements to MDV's spatial converter.
2. Build the smallest `generate_synthetic_spatial_project` command.
3. Add a cleanup selector for `cleanup_group == "synth-spatial"`.
4. Add a Playwright performance smoke spec for one small generated project.
5. Review `dummy-anndata` against `mock_anndata.py` and remove only the local
   code it can replace cleanly.
