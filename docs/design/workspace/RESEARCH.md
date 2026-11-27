# Workspace research

This is the consolidated evidence behind the workspace plan. It replaces the former research
and ticket files so the active plan is navigable. Git history retains the exhaustive
line-by-line inventories. Recheck code references before implementing because the branch
continues to move.

## Executive summary

- A current view is organised by data-source pane. Layout is per pane; charts carry pixel
  `size`/`position` and often grid `gssize`/`gsposition` as well.
- The dev-server sample contained 197 distinct views: 74% had more than one pane, with 379
  absolute and 162 gridstack panes; 14 distinct views mixed both. Median chart count was 7 and
  the maximum was 23.
- The false “unsaved grid view” behaviour is explained by MDV's integration, not a necessary
  gridstack bug: the baseline is captured before layout settles, pixel sizes depend on the
  window, and collision-pushed charts are not reliably written back.
- All charts share `BaseChart` chrome. Most settings can be described by `getSettings()`, but
  spatial, Viv, genome browser, table and gate charts have custom UI.
- Theme ownership is fragmented across legacy project CSS, a drifted catalog copy, Tailwind /
  shadcn, two MUI themes and roughly 470 hard-coded colours. Light/dark starts from the OS and
  is not persistently owned by one system.

## Current view and layout model

A view is broadly:

```text
{
  dataSources: {
    <name>: { layout, panelWidth, highlight }
  },
  initialCharts: {
    <name>: [chart config]
  },
  links,
  viewImage,
  name
}
```

`ChartManager._init` creates one Split.js pane per `dataSources` key. A missing pane layout
means absolute. Gridstack uses its default 12 columns. Placement fields are:

| Field | Meaning |
|---|---|
| `size` | chart width/height in pixels; required by the base schema |
| `position` | absolute left/top in pixels; used but absent from the schema |
| `gssize` | width/height in grid units |
| `gsposition` | x/y in grid units |
| `location` | legacy overloaded placement, converted when a chart is added |

Important front-end ownership:

- `src/charts/ChartManager.js`: builds panes and menus, loads views, creates charts, switches
  absolute/gridstack layouts and serialises view state.
- `src/charts/GridstackManager.ts`: creates one grid per data source; wraps chart DOM; observes
  size and style; writes `gssize`/`gsposition`; owns compact and lock UI.
- `src/charts/ViewManager.ts`: changes views and coordinates saved-state comparisons.
- `src/charts/BaseChart.ts`: writes pixel size and handles full screen.
- `src/utilities/Popout.ts`: removes and restores layout wrappers.
- `src/charts/schemas/ChartConfigSchema.ts`: only partially describes the persisted shape.

Python has many independent readers and writers: default/template views, spatial prototypes,
ChatMDV normalisation, conversion utilities, `MDVProject.set_view`, examples and fixtures.
Most writers omit `layout`, producing absolute panes. Default views, ChatMDV,
`create_view_with_all_datasources`, the image prototype and parts of the spatial template write
gridstack. The spatial template mixes gridstack and absolute panes.

Contracts a new format must account for include rename/delete, column deletion's fail-closed
scan, JSON/Icechunk history, thumbnails, static output, links, highlights and Python-generated
views. A formal schema must cover fields current Zod schemas omit and provide a derivation or
validation path for Python rather than another hand-maintained model.

## Why current gridstack views drift

Three mechanisms are visible in MDV:

1. `ViewManager` stores the saved baseline before `GridstackManager` observers and deferred
   chart sizing have settled.
2. Pixel `size` is retained for gridstack charts, while row height depends on container/window
   state and on one shared manager value.
3. MDV watches inline `style`, but gridstack 7 writes collision-pushed positions to attributes;
   pushed charts can retain stale `gsposition`.

Gridstack provides a `change` event and grid-unit `save()` that MDV does not use. Gridstack 14
also changed its rendering and responsive APIs, but upgrading alone would not fix baseline
timing or persistence ownership.

## Layout candidates

### Gridstack 14

- Native collision, locking, touch and breakpoint support; existing DOM can be adopted with
  `makeWidget`.
- No built-in keyboard move/resize and no simple persisted set of all breakpoint layouts.
- Upgrading from MDV's 7.3 touches `oneColumnSize`, `float`, touch behaviour and integration
  assumptions.
- Approximately 97 KB minified / 27 KB gzip for `gridstack-all.js` in the researched version.

### react-grid-layout 2.3

- Responsive layout objects and comparatively little adapter code.
- Needs a React wrapper around legacy chart DOM and did not provide keyboard movement.
- The layout lab found substantially worse drag/resize worst frames than the other candidates.

### In-house prototype engine

- Small `compact`, `flow` and `findSpot` core already matches the prototype and can grow the
  canvas while dragging.
- Requires MDV to own pointer/touch gestures, accessibility, handles, tests and maintenance.
- The layout lab proved that it can host real chart DOM when kept behind a narrow seam.

The engine comparison did not answer the product-level question of scrolling grid versus
tiling or a hybrid. That is the next checkpoint in [LAYOUT-EVALUATION.md](LAYOUT-EVALUATION.md).

## Chart chrome, settings, tools and filters

`BaseChart` supplies the common title bar and actions such as settings, duplicate and download.
Only six registered types add title-bar icons. Gridstack adds lock; `ChartManager` adds pop-out
and close.

The context menu contains chart actions independently of the workspace layout. Most of the 51
registered chart types are representable through `getSettings()`. Exceptions and specialised
paths include spatial layers, Viv, genome browsers, table configuration and gates. The
showcase therefore cannot assume that docking the generic settings renderer replaces every
dialog.

`ChartManager.addMenuIcon` adds per-data-source or nominal main-menu DOM. Its only live source
caller found in the inventory was `GridstackManager`; `hyperion.js` was not imported.
`chart.addMenuIcon` has real per-chart callers, so extension actions need an explicit migration
path rather than a React-state replacement that ignores callers expecting returned elements.

Selection tools are not one uniform abstraction. `SelectionOverlay` serves Deck scatter and
Viv/spatial region charts, while other charts have their own mechanisms.

There are roughly a dozen filter mechanisms. Most can be described in words, but few expose a
saved, uniform contract. A workspace filter presentation needs at least:

- a human summary
- the owning data source
- the passing rows or count
- a clear action

Linked-data-source and colour-legend filters, plus charts that cannot describe a filter, remain
explicit design questions in [PLAN.md](PLAN.md).

## Theme and colour ownership

Four systems do not share one source of truth:

1. Legacy variables in `src/charts/css/charts.css`, used by project and static pages.
2. A drifted variable set in `src/catalog/catalog_index.css`.
3. Tailwind/shadcn tokens in `src/tailwind.css` and the Tailwind configuration.
4. MUI themes: `ThemeProvider.tsx` for catalog/login and `MaterialWrapper` for project UI.

Canvas, WebGL, SVG and image-export code also reads colours directly in JavaScript. Replacing
CSS variables alone cannot provide complete light/dark behaviour.

On the project page `ChartManager.setTheme` applies the `<html>` class; catalog/login use their
own provider. Both begin from the OS setting. A manual choice is not consistently persisted and
an OS change can reset it. Static pages depend on the project CSS/theme path rather than the
catalog provider.

Data colours are a separate domain: column palettes live in datasource configuration and
`DataStore`, with linked-source colour synchronisation. They should not be conflated with UI
tokens.

The token decision must establish one published source, mappings for legacy variables and MUI,
ownership of the stored light/dark preference, and a migration approach for hard-coded drawing
colours.

## Evidence gaps

- Browser reproduction of every current gridstack persistence bug.
- Touch, keyboard and screen-reader behaviour across layout candidates.
- Complete handling of specialised settings and every filter mechanism.
- Final static-build support for the workspace format.
- The long-term `datasources_v2` / SpatialData / AnnData hydration lifecycle.

These are not reasons to expand the initial showcase. They are constraints and checkpoints in
[PLAN.md](PLAN.md).
