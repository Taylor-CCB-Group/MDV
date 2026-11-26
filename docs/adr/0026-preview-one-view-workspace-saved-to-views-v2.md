# Preview a one-view workspace behind a preference, saving to `views_v2.json`

**Status:** proposed. Settles the view-shape questions for the workspace redesign (the
workspace prototype); the layout engine is chosen separately.

## Context

A view is split into one pane per data source. `ChartManager._init` builds a pane, a menu bar
and a layout (absolute or gridstack) for each key of `view.dataSources`, and a chart's data
source is recorded only by which `initialCharts[dataSource]` array it sits in. On the dev
server, 74% of the 197 distinct views have more than one pane, 379 panes are absolute and 162
gridstack, and 14 views mix the two. Absolute panes store pixel `position` / `size`, so a
layout made on a wide screen runs off a narrow one.

The redesign puts every chart on one grid, with the data source chosen in Add chart and
settings in a docked inspector (`docs/design/spatial-tables/04-views-datasources-chartmanager.md`,
option B2). ADR 0023 makes each table its own data source, so views will mix more of them.

`views.json` has many readers and writers: `save_state`, `MDVProject.set_view`, `rename_view`,
`delete_view`, the Python template writers, ChatMDV, `/get_view`, the catalog thumbnail,
`soft_delete_column` (ADR 0003), `column_name_report.py`, `project_json_history.py` and the
static build. Changing its shape in place would touch all of them before the new workspace has
been tried by anyone.

## Decision

### An experimental preference

The new workspace is a browser setting in `src/lib/browserSettings.ts`, marked `preview` so it
is listed while Feature previews is on, as Convert legacy charts is. Off, views open in the
pane workspace exactly as today. On, every view opens in the new workspace. The pane workspace
stays until the team retires it.

### One flat grid

A view's charts share one grid with no section per data source. Each chart card shows its data
source. The pane menu bars are replaced as in the prototype: Clear all filters in the filter
bar, one row count per data source in the status bar, settings in the inspector, and Add chart
asks for the data source first. The absolute/gridstack layout menu goes.

**A view's data sources** are the keys of the view's `dataSources` when it opens, including a
pane that had no charts, plus any data source a chart is added for. None is dropped while the
view is open. Highlighted rows are saved per data source, filter events reach `ChartManager`
listeners only for data sources that have been shown (ADR 0019), and an interaction link adds
an image chart only for a data source in the view; all three keep working on this set.

### Converting a pane view

A view with no v2 entry is converted in memory as it opens. Charts are taken pane by pane in
the old left-to-right order, and within a pane top to bottom, then left to right. Each keeps
its width relative to its pane and is packed row by row into 12 columns with no gaps. The same
view always converts to the same layout. As with legacy charts, an alert says the view was
converted and that saving keeps it.

### Saving to `views_v2.json`

Saving in the new workspace writes the view to `views_v2.json`, keyed by the same view name.
It never writes `views.json`. With the preference on, a view opens from `views_v2.json` when it
has an entry there, and is converted from `views.json` otherwise. With the preference off, the
`views.json` view opens unchanged. From the first save the two copies are independent.

```json
{
  "version": 2,
  "dataSources": { "cells": { "highlight": [12, 40] }, "genes": {} },
  "charts": [
    { "id": "a1", "dataSource": "cells", "config": { "type": "wgl_scatter_plot" } },
    { "id": "b2", "dataSource": "genes", "config": { "type": "row_chart" } }
  ],
  "layout": {
    "a1": { "x": 0, "y": 0, "w": 8, "h": 6 },
    "b2": { "x": 8, "y": 0, "w": 4, "h": 6 }
  },
  "links": [],
  "viewImage": "..."
}
```

- Each chart names its data source; a chart config holds no position.
- `layout` is in grid units only: 12 columns, and rows whose height the format version fixes.
  No pixels are stored.
- Only the 12-column layout is saved. Narrower layouts are always derived from it (below).

A "Copy layout JSON" action lets people trying the preview send back arrangements that worked
or broke.

### Everything else that reads views

| Code | During the preview |
| --- | --- |
| Python template writers, ChatMDV, `set_view` | keep writing `views.json`; their views are converted when opened |
| `rename_view`, `delete_view` | act on both files, so a view name always means one view |
| `soft_delete_column` | scans both files; ADR 0003's scan fails closed and must see v2 charts |
| `project_json_history.py`, and Icechunk history (ADR 0023) | cover both files |
| View list | a view that exists only in `views_v2.json` is listed only while the preference is on |
| Catalog thumbnail, static build | keep reading `views.json` |

### Screen sizes

The grid has 12 columns from 960 px, 6 from 600 px and 1 below. The 6- and 1-column layouts
are derived from the saved 12-column layout in reading order and are never stored.

| Width | Columns | What works |
| --- | --- | --- |
| 960 px and up | 12 | everything |
| 600–959 px | 6 | everything except moving and resizing charts; the inspector slides over the grid; added charts go at the end |
| under 600 px | 1 | viewing, filtering, selecting and full screen |

## Why a preference, when ADR 0019 rejected flags

ADR 0019 rejected keeping two implementations behind a flag for a single extraction, which is
easy to revert as one commit. This replaces the whole workspace and the view shape, and its
design is still settling, so it needs to be tried on real projects before the pane workspace
goes. The preference follows a pattern the app already has, and saving to a separate file means
nobody with the preference off is affected by anyone who has it on.

## Considered and rejected

- **Keep old views on the pane workspace and open only new views in the new one**, with a
  `kind` field (the design note's proposal). Existing views could not be tried in the new
  workspace, and they are what has to work.
- **Convert in memory and never save.** Settled first, then reversed: the preview has to be
  usable for real work, and a format nobody saves is never exercised.
- **Add a layout block inside each `views.json` view**, keeping charts under
  `initialCharts[dataSource]`. Both workspaces would share one chart list, but the new one's
  flat model would be squeezed into the per-data-source shape, and every reader of
  `views.json` would meet a block it does not understand.
- **Save a tablet layout once someone rearranges at 6 columns** (as the prototype does). It
  gives tablets full editing, but the two layouts drift: later desktop edits no longer reach
  the tablet one, and adding or removing a chart needs rules for both. A tablet layout can be
  added later as a new key without breaking saved views.
- **Keep each chart's position in its config**, as `gsposition` is today. It mixes layout into
  chart settings and leaves no single place for a second layout if one is ever added.
- **Rewrite every `views.json` in one batch at upgrade time.** Nothing can go back, and every
  reader has to change at once.
- **Have the Python writers write v2 now.** More work before the format has been tried; their
  views convert on open meanwhile.

## Consequences

- A view can have two copies that differ. People with the preference off see the
  `views.json` copy. This ends when `views.json` is retired, which needs a conversion of the
  remaining views and the readers in the table above moved to v2.
- Converted views that are never saved keep converting on every open, so the converter must
  stay deterministic and is the upgrade path when the pane workspace is retired.
- The format is independent of the layout engine; the engine maps `{x, y, w, h}` onto its own
  model.
- Tablets cannot move or resize charts. Charts added on a tablet go at the end of the layout.
- The static build cannot save, and it shows the `views.json` view.
- Both files are last-writer-wins, as `views.json` already is.
