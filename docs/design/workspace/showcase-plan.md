Map: [Workspace in MDV - Map](README.md)
Builds on: [18 One view, saved layout and screen sizes](tickets/18-one-view-saved-layout-and-screen-sizes.md) (ADR 0026), [19 Layout engine bake-off and choice](tickets/19-layout-engine-bake-off-and-choice.md) (ADR 0027)
Status: draft, 2026-11-27

# Workspace showcase plan

## Goal

Real MDV, with real projects and real charts, wearing the [workspace prototype](prototype.html)'s look and feel. It needs to be good enough for the team to answer one question: **is this the right direction for MDV?**

It is not a conversion of MDV. Anything that doesn't change how the workspace looks or feels is left out or stubbed, and the long tail waits for the team's verdict.

**Time:** about a week (five working days), on a branch, never merged until the team has tried it and we have iterated.

## What the showcase achieves

Someone opens a real project with the Workspace preview on and sees:

- **One view, one grid.** Every chart of the view, from every data source, on one 12-column canvas in the prototype's cards. No panes and no per-data-source menu bars.
- **The prototype's handling.**
  - Drag a chart by its header and the others slide aside, with a dashed placeholder showing where it lands.
  - Resize from the edges, including past the bottom of the window.
  - Maximise a chart, lock it, duplicate it, delete it.
  - Alt+arrows move the focused chart.
- **The prototype's frame:**
  - app bar with project title, saved state, menus and a light/dark toggle
  - filter bar
  - tool rail
  - docked inspector
  - status bar with view tabs and row counts
- **Settings in the inspector, not a floating dialog.** Selecting a chart shows its settings in the side panel.
- **Add chart from one place.** Pick the type and the data source, and the chart lands in the first free spot.
- **Saved and reloaded.** Save writes the one-grid layout to `views_v2.json`. Reloading gives the same layout. The old `views.json` is never touched.
- **Light and dark.** The prototype's colours, inside the workspace.
- **Narrower screens.** The grid reflows to 6 columns on tablets and 1 on phones (view only), as in the prototype.

What the team should come away able to judge:
- the overall feel
- one grid against panes
- the inspector against the floating dialog
- the filter bar
- the saved format
- whether the effort is worth it

## Prototype feature by feature

| Prototype piece | In the showcase | How |
|---|---|---|
| Grid, cards, drag, resize, placeholder, auto-scroll | Yes | The layout lab's in-house engine, already working on real charts |
| Card header: type icon, title, data source badge, filtering pill, actions | Yes | New header around the chart's own body. The chart's existing icons move into the header's menu |
| Maximise, lock, duplicate, delete, hide | Yes | Maximise and lock are workspace state. Duplicate and delete call the existing `ChartManager` actions |
| App bar: title, saved state, menus | Yes | Menus call existing actions: save view, save as, upload data, add chart, new view, reset filters |
| Light/dark toggle | Yes | Prototype colour tokens, scoped to the workspace |
| Inspector: chart settings | Yes | The existing settings component, docked instead of floating |
| Inspector: chart list and data sources | Yes | List of the view's charts with hide/show. Data source cards with filtered counts and columns |
| Filter bar | Yes, simplified | Rows passing filters per data source, a chip for each chart we can tell is filtering, and Reset all. The full chip-per-filter picture waits (about a dozen filter mechanisms today) |
| Tool rail: select, lasso, pan, zoom | Partly | Select and pan map to the existing "what dragging does" setting. Lasso only where the chart already supports it. Zoom left out |
| Status bar: view tabs, counts | Yes | View tabs switch views through the existing view manager |
| Add chart (type, then data source, then fields) | Yes, reused | The existing Add chart dialog in the prototype's styling, placing the new chart with `findSpot` |
| Toolbar acting on the selection | Basic | Chart title, data source, settings, duplicate, delete. Per-type quick fields wait |
| Undo and redo | No | Layout-only undo is cheap and can be added if wanted. Undo of chart settings is a project of its own |
| Command palette (⌘K), keyboard shortcuts sheet | No | Nice to have |
| Share, account avatar | No | Not part of the look-and-feel question |
| Responsive chrome (inspector as a drawer, rail at the bottom on phones) | Basic | The grid reflow is done. The chrome gets the prototype's main breakpoints only |
| Design notes and the prototype's control bar | No | Prototype-only |

## How MDV changes

**With the Workspace preview off (the default): nothing changes.** The pane workspace, `views.json`, the static build and every existing view behave exactly as today.

**With it on** (Settings, Feature previews, Workspace preview, held in this browser only):
- Opening a project shows the new workspace instead of the panes.
- A view saved in the new format opens from `views_v2.json`. Any other view is converted from its pane view as it opens, and a note says saving keeps it.
- Saving writes only `views_v2.json`. Turning the preview off shows the pane view as it was.
- Charts are the same chart objects. Only their frame, settings panel and placement change. Filtering, highlighting and selection between charts work as today.

**For the code:** the new workspace lives in its own folder, and the existing code gets a few small hooks. Nothing in the pane workspace is removed or rewritten.

## Code changes

All on this branch, `docs/workspace-redesign-plan` (PR #719). The layout lab already sits here in `src/react/prototypes/layoutLab/`: the showcase builds from it, and the lab is deleted once the showcase has replaced it.

### New: `src/react/workspace/`

| File | What it does |
|---|---|
| `layout/engine.ts` | `compact`, `flow`, `findSpot`, breakpoints, from the lab. Unit-tested |
| `layout/Grid.tsx` | The in-house engine behind the ADR 0027 seam: takes charts and the 12-column layout, reports the new layout after each gesture |
| `ChartCard.tsx` | Card frame and header, and the lab's `ChartHost`, which moves the chart's own div into the card and resizes it after a gesture |
| `Workspace.tsx` | The frame: app bar, toolbar, filter bar, rail, canvas, inspector, status bar |
| `Inspector.tsx` | Tabs for chart settings (the existing `SettingsDialogComponent` docked), the chart list and data sources |
| `FilterBar.tsx`, `ToolRail.tsx`, `ViewTabs.tsx` | As in the prototype, wired to the existing data stores, drag setting and view manager |
| `convert.ts` | Pane view to v2 (ADR 0026), from the lab's `convertPaneView` |
| `viewStore.ts` | Loads and saves v2 views, and tracks unsaved changes |
| `workspace.css` | The prototype's tokens and styles, scoped to the workspace root, with the dark variant |

### Changed, small hooks

| File | Change |
|---|---|
| `src/lib/browserSettings.ts` | Add `workspacePreviewSetting` (off/on), listed under Feature previews, like Convert legacy charts |
| `src/modules/project_bootstrap.tsx` | With the preview on, mount the workspace once `ChartManager` is ready, and hide the pane layout |
| `src/charts/ChartManager.js` | Let the workspace take a view's charts when it loads, place added and duplicated charts, and produce v2 state for saving |
| `src/charts/BaseChart.ts` | Settings opens in the inspector when the workspace is active, otherwise the floating dialog as now. Remove the pane drag and resize from hosted charts |
| `python/mdvtools/mdvproject.py` | `views_v2` read and write next to `views`, with `set_view_v2` |
| `python/mdvtools/server.py`, `serverlite.py` | Save and load routes for v2 views. `rename_view` and `delete_view` act on both files (ADR 0026), plus a test |

## Day by day

| Day | Work | Done when |
|---|---|---|
| 1 | Branch. Preview setting. Mount the workspace with the lab's grid and cards. Pane view conversion. Panes hidden | A real project opens with the preview on, showing every chart in one grid you can drag and resize |
| 2 | The frame: app bar and menus, status bar with view tabs, tool rail, light/dark tokens. Card header with maximise, lock, duplicate, delete | It looks like the prototype around real charts, in both themes |
| 3 | Inspector: docked settings, chart list, data sources. Toolbar on the selection | Selecting a chart edits it in the side panel. No floating dialog |
| 4 | Saving and loading v2 (Python routes and front end). Add chart. Filter bar. Narrow widths | Save, reload, same layout. Add a chart from any data source. Filter chips and Reset all work |
| 5 | Run on the test projects (spatial, 23 charts, 3 data sources). Fix what looks wrong. Quick automatic check (no overlaps, save and reload unchanged). Demo build | The team can run it |

## What gets checked

- On the three lab projects, at 1440, 1024, 768 and 390 px: no overlaps, no chart outside its card, save and reload unchanged (the lab's script, pointed at the workspace).
- `pnpm exec tsgo`, and the Python tests for the new routes.
- Spatial charts don't reload their geometry on a drag (`perf.ts`, turned off afterwards).
- A walk through every menu item and card action, in light and dark.

## Risks

- **Real charts misbehaving in the new frame.** This derailed the last attempt. The lab found and fixed the main causes (legacy drag, stray resize handles, overlays, redraws mid-gesture). Anything new is fixed in the card host, not per chart.
- **Settings components that assume a floating dialog.** Spatial, Viv, genome browser, table and gate charts have their own dialogs. They keep opening those in the showcase.
- **Hiding the pane layout cleanly.** The lab drew over the panes. The showcase has to stop them taking space and events, without breaking chart creation, which still goes through the panes.
- **The week slipping.** Days 1–3 are the look and feel. If time runs short, the filter bar and narrow widths shrink first, not the frame.

## After the showcase

- The team runs the branch and reacts. We iterate on the branch, as many rounds as needed, before anything merges.
- Decisions the team weighs in on with the real thing in front of them:
  - the saved format
  - the inspector for charts with their own dialogs
  - how filters are shown
  - MDV-wide theming
  - undo
  - when (and whether) the preview merges, then becomes the default
- Only then: the PRD and build issues for the rest, if the verdict is yes.

## Changes to the map if this plan is accepted

- [20 Walking skeleton](tickets/20-walking-skeleton.md) becomes this showcase.
- [21 Chart chrome, inspector, tools and filters](tickets/21-chart-chrome-inspector-tools-and-filters.md) and [22 Theme, undo and rollout](tickets/22-theme-undo-and-rollout.md) shrink to choices made while building (following the prototype) and questions for the team afterwards.
- The destination becomes "a showcase the team can judge", with the PRD and build issues only after a yes.
- Notes: nothing merges to main until the team has tried the branch and we have iterated.
