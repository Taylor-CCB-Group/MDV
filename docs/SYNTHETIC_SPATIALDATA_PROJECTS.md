# Synthetic SpatialData Projects

This note records the reusable conventions for generated SpatialData-backed MDV
projects. The aim is repeatable local projects for chart development,
performance profiling, and agent-driven checks without committing large
fixtures.

For scatter/table performance background, see
[issue #433](https://github.com/Taylor-CCB-Group/MDV/issues/433).

## At A Glance

- Use generated projects for repeatable local chart, spatial, and performance
  work without checking large fixtures into git.
- Keep the first scope modest: representative columns and common chart families,
  not every obscure chart configuration.
- Put interactive generated projects directly under `~/mdv`; nested folders in
  `~/mdv` are not supported by the project scanner.
- Treat generated projects as disposable. Use exact path cleanup only, and keep
  provenance metadata so later cleanup can be automated safely.
- For new React-backed chart configs, prefer creating the view through the UI
  and then inspecting `views.json` before encoding the shape in Python.
- Two serving modes are useful for performance work: the running container
  (main-branch frontend) and the local Python venv (current-branch frontend).
  See [Serving Generated Projects](#serving-generated-projects).

## Project Convention

Use flat, disposable paths:

```text
~/mdv/synth-spatial--<profile>--<size>--<layout>/
```

Examples:

```text
~/mdv/synth-spatial--scatter-table--100k--single/
~/mdv/synth-spatial--scatter-table--1m--single/
~/mdv/synth-spatial--scatter-table--staggered-regions/
```

The generated MDV project should live at that final path. SpatialData source
stores may live inside the project, for example:

```text
~/mdv/synth-spatial--scatter-table--100k--single/spatial/source.zarr/
```

That keeps the project portable and lets cleanup remove one subtree.

Generated projects include `state.json.provenance.cleanup_group =
"synth-spatial"`. Agents should only remove exact direct-child paths matching
the `~/mdv/synth-spatial--...` convention, and should avoid broad cleanup of
`~/mdv`.

Example provenance:

```json
{
  "provenance": {
    "created_by": "mdvtools synthetic spatialdata generator",
    "generator": "dummy-spatialdata",
    "profile": "scatter-table",
    "n_cells": 100000,
    "n_genes": 50,
    "image_size": 512,
    "seed": 42,
    "n_coordinate_systems": 1,
    "coordinate_system_cell_counts": [100000],
    "cleanup_group": "synth-spatial"
  }
}
```

See [Serving Generated Projects](#serving-generated-projects) for how to make
a generated project visible to the running backend or to the local Python venv.

## Identifying Existing Synthetic Projects

The directory naming convention is the quickest filter:

```bash
ls -d ~/mdv/synth-spatial--*/ 2>/dev/null
```

The provenance metadata is the authoritative filter and also catches projects
generated with a custom `--output` path:

```bash
for f in ~/mdv/*/state.json; do
    if jq -e '.provenance.cleanup_group == "synth-spatial"' "$f" >/dev/null 2>&1; then
        dirname "$f"
    fi
done
```

Relevant `state.json` fields are documented above and include `profile`,
`n_cells`, `n_genes`, `n_coordinate_systems`, and
`coordinate_system_cell_counts`, so a one-liner can list size/shape per
project.

## Serving Generated Projects

There are two useful modes; pick deliberately when measuring performance.

### 1. Through the running container (main-branch frontend)

The local docker compose stack defined in `docker-secrets.yml` bind-mounts
host `~/mdv` to the container's `/app/mdv`, so any direct child of `~/mdv`
becomes visible to the running backend on `localhost:5055`. The frontend the
container serves is whatever is built into its image, i.e. the **main**
branch, not the current worktree. The Vite dev server on `localhost:5170`
proxies its API calls to `localhost:5055`, so it pairs naturally with the
container backend when iterating on frontend code.

After generating a project, register it without restarting the container:

```bash
curl -fsS http://127.0.0.1:5055/rescan_projects -o /dev/null
```

That route is `GET` and responds with a `302` redirect to `/`; the backend log
line `RESCAN PROJECTS` confirms it ran. Once registered, the project is
available at `http://localhost:5055/project/<id>` (and the same path through
`localhost:5170` via the proxy).

### 2. Through the local Python venv (current-branch frontend)

Use this mode to verify the **current branch's** frontend against a generated
project (apples-to-apples performance comparisons, build-vs-dev checks, etc.).

```bash
pnpm run build-flask-vite
venv/bin/mdvtools serve ~/mdv/synth-spatial--scatter-table--100k--single --port 5057
```

`build-flask-vite` writes the bundle to `python/mdvtools/static/`, which the
local `mdvtools serve` CLI (using `serverlite.py` and `templates/page.html`)
serves directly. Use a port that does not clash with the container (5055) or
Vite (5170). Multiple concurrent instances on different ports are fine.

`python/mdvtools/mdv_desktop.py`, which would otherwise scan `~/mdv` and serve
every project, currently hardcodes `project_dir = "/app/mdv"` for the Docker
deployment and is not usable as-is for local multi-project serving. Until that
is parameterised, prefer per-project `mdvtools serve` invocations locally.

## Cleanup

Filesystem-only cleanup of synthetic projects (use this when the project was
never registered with the DB-backed backend, or you also intend to run the DB
purge step below):

```bash
rm -rf ~/mdv/synth-spatial--*
```

DB-backed cleanup goes through `cleanup_projects.py`, which removes both the
filesystem directory and the related rows in `projects`, `files`, and
`user_projects`. The `test-generated` selector matches any project whose
`state.json` has `provenance` (which includes synthetic spatial projects).
Combine it with `--filter-name` to restrict the purge:

```bash
# Preview
venv/bin/python -m mdvtools.dbutils.cleanup_projects \
    --list test-generated --filter-name "synth-spatial--*"

# Purge
venv/bin/python -m mdvtools.dbutils.cleanup_projects \
    --purge test-generated --filter-name "synth-spatial--*" --confirm
```

See `python/mdvtools/dbutils/CLEANUP_SCRIPT_README.md` for selector and
filter semantics. Note that `cleanup_projects` connects to the database the
same way `mdv_server_app` does, so it must be run with credentials that match
the running backend (typically by executing it inside the running container,
not from the host venv).

Future improvement: add a dedicated cleanup selector keyed on
`state.json.provenance.cleanup_group == "synth-spatial"` so cleanup is not
conflated with other test-generated projects.

## Generator

Run from `python/` using the local venv:

```bash
../venv/bin/python -m mdvtools.tests.generate_synthetic_spatial_project \
  --profile scatter-table \
  --n-cells 100000 \
  --seed 42 \
  --output ~/mdv/synth-spatial--scatter-table--100k--single \
  --force
```

Important options:

- `--profile`: currently `scatter-table` or `spatial-overview`;
- `--n-cells`: total row count when per-coordinate-system counts are not given;
- `--n-genes`: generated expression-like variable count;
- `--n-coordinate-systems`: number of image-backed coordinate systems;
- `--coordinate-system-cell-counts`: comma-separated per-region counts such as
  `1k,5k,10k,25k,50k,100k,250k,500k,1m,2m`; when supplied, `--n-cells` and
  `--n-coordinate-systems` are inferred;
- `--image-size`: generated image extent in pixels;
- `--output`: flat `~/mdv` project path. When omitted, the default path
  includes both total size and coordinate-system layout to avoid collisions
  between single-region and staggered projects;
- `--force`: replace the exact output path.

For small single-region projects, the generator uses `dummy-spatialdata` shapes.
For larger `scatter-table` projects, it uses an image-annotated table rather
than one shape per row, keeping the project SpatialData-backed without
million-feature GeoJSON overhead.

The project includes enough categorical, numeric, coordinate, image/region, and
table columns to exercise common scatter, row-chart, and table workflows. A
later pass should review whether `dummy-anndata` can replace part of
`python/mdvtools/tests/mock_anndata.py`.

The staggered multi-region stress shape is:

```bash
../venv/bin/python -m mdvtools.tests.generate_synthetic_spatial_project \
  --profile scatter-table \
  --n-genes 2 \
  --image-size 1024 \
  --coordinate-system-cell-counts 1k,5k,10k,25k,50k,100k,250k,500k,1m,2m \
  --output ~/mdv/synth-spatial--scatter-table--staggered-regions \
  --force
```

That produces 3,941,000 total cells with small, medium, and large region-local
working sets in one project.

## Views And Chart Types

Generated `scatter-table` projects include:

- `Scatter/table comparison`: row chart, classic scatter, React/deck scatter,
  classic table, and React table;
- `react view`: React/deck scatter plus React table;
- `classic view`: classic scatter plus classic table.

The UI-authored chart type mapping currently used for these views is:

- `Row Chart` -> `row_chart`;
- `2D Scatter Plot` -> `wgl_scatter_plot_dev`;
- `2D Scatter Plot (Classic)` -> `wgl_scatter_plot`;
- `Table` -> `table_chart_react`;
- `Table (Classic)` -> `table_chart`.

Although the newer scatter is implemented by
`DeckScatterReactWrapper.tsx`, the persisted chart type is currently
`wgl_scatter_plot_dev`.

When expanding generated views, use this loop:

1. Generate or copy a small flat project under `~/mdv`.
2. Serve it editably with `MDVProject(...).set_editable(True)`.
3. Use Playwright or the browser UI to add charts through the normal Add Chart
   dialog.
4. Save the view, inspect `views.json`, and promote only stable useful config
   into the generator.

## Profiling

Use Playwright from the repo root. For focused ad hoc row-chart timing against
an already-served project:

```bash
node scripts/profile_mdv_interactions.mjs \
  --url http://localhost:5055/project/191 \
  --views "react view,classic view" \
  --categories type_0,type_1,type_2 \
  --out output/playwright/project-191-row-filter-profile.json
```

The profiler records:

- time from row-chart category click to the datastore `filtered` event;
- time to materialise `DataStore.getFilteredIndices()`;
- time after two animation frames as a coarse visual-settle marker;
- chart/view metadata and browser console warnings/errors.

Early observations from 1m-row profiling:

- category filter updates themselves were usually tens of milliseconds;
- React/deck chart updates were dominated by filtered-index materialisation;
- the first cold React click was much slower than later clicks, so warm/cold
  runs should be separated.

## Follow-Up Leads

Issue #433 captures the broader comparison:

- classic scatter re-renders external filter changes quickly, but panning and
  zooming can degrade badly once filters are active;
- React/deck scatter has more overhead on filter updates, and continuous
  selection-region updates are more expensive than the classic mouse-up flow;
- deck rendering is generally reasonable, but zoomed-out overdraw with many
  points per pixel needs attention;
- classic scatter has visual legibility and trackpad/camera-interaction issues;
- deck.gl's WebGPU path is worth exploring as it matures.

Memory/performance areas to audit:

- `DataStore.filterBuffer`, `Dimension.filterArray`, and
  `DataStore.getFilteredIndices()` allocations;
- duplicated row-index/filter-state buffers in React/deck paths;
- region-local filtering so spatial views do not materialise datasource-wide
  filtered indices unnecessarily;
- a lazy gating implementation: the current startup path eagerly creates a
  `GateManager`, which creates an empty `__gates__` multitext membership buffer
  of `rows * 24 * 2` bytes. Avoid metadata-only generated-project workarounds;
  fix this in application code with gating regression coverage.

Useful code entry points for that audit:

- `src/react/components/DeckScatterReactWrapper.tsx`;
- `src/charts/WGLScatterPlot.js`;
- `src/charts/TableChart.js`;
- `src/react/components/TableChartReactWrapper.tsx`;
- `src/react/components/SelectionDialogComponent.tsx`;
- `src/datastore/DataStore.js`;
- `src/datastore/Dimension.js`;
- `src/react/gates/GateManager.ts`.

Future passes:

- add an explicit cleanup selector for `cleanup_group == "synth-spatial"`;
- add a small Playwright performance smoke spec;
- review whether `dummy-anndata` can replace part of
  `python/mdvtools/tests/mock_anndata.py`.
