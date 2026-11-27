Map: [Workspace redesign](README.md)
Builds on: [ADR 0026](../../adr/0026-preview-one-view-workspace-saved-to-views-v2.md)
and [ADR 0027](../../adr/0027-build-the-workspace-grid-on-the-prototypes-own-engine.md)
Status: draft, 2026-11-27

# Workspace showcase plan

## Goal

Real MDV, with real projects and real charts, wearing the
[workspace prototype](prototype.html)'s look and feel. It must help the team answer two
questions:

1. **Product:** is this the right workspace direction for MDV?
2. **Architecture:** can it be introduced through a bounded compatibility adapter without
   making `ChartManager` own more of the application?

It is not a conversion of MDV. Anything that doesn't change how the workspace looks or feels is left out or stubbed, and the long tail waits for the team's verdict.

**Time:** about a week (five working days), on a branch, never merged until the team has tried it and we have iterated.

## What the showcase achieves

Someone opens a real project with the Workspace preview on and sees:

- **One view, one initial grid.** Every chart of the view, from every data source, on the
  in-house grid in the prototype's cards. No panes and no per-data-source menu bars. Tiling or
  a hybrid remains a checkpoint, not a rejected option.
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
- **Converted without commitment.** Pane views convert deterministically into an in-memory
  workspace model. Copy layout JSON is available for feedback. Durable save/reload follows the
  persistence checkpoint.
- **Light and dark.** The prototype's colours, inside the workspace.
- **Narrower screens.** Desktop and responsive laptop/tablet behaviour are part of the
  showcase. Phones need a safe basic fallback; a full mobile workspace can follow later.

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
| Grid, cards, drag, resize, placeholder, auto-scroll | Yes, baseline | The layout lab's in-house engine, already working on real charts; compared with a small tiling or hybrid experiment |
| Card header: type icon, title, data source badge, filtering pill, actions | Yes | New header around the chart's own body. The chart's existing icons move into the header's menu |
| Maximise, lock, duplicate, delete, hide | Yes | Maximise and lock are workspace state. Duplicate and delete go through the legacy chart adapter |
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
| Responsive chrome | Basic | Desktop, laptop and tablet breakpoints; phones get a non-breaking single-column/read-only fallback |
| Design notes and the prototype's control bar | No | Prototype-only |

## How MDV changes

**With the Workspace preview off (the default): nothing changes.** The pane workspace, `views.json`, the static build and every existing view behave exactly as today.

**With it on** (Settings, Feature previews, Workspace preview, held in this browser only):
- Opening a project shows the new workspace instead of the panes.
- A pane view converts deterministically into the in-memory workspace model as it opens.
- Until the persistence checkpoint is resolved, the old `views.json` is never changed and Copy
  layout JSON is the only export.
- Charts are the same chart objects. Only their frame, settings panel and placement change. Filtering, highlighting and selection between charts work as today.

**For the code:** `project_bootstrap` selects the preview path early. A workspace session owns
the current view, layout, selected chart, dirty state and persistence boundary. Existing chart
machinery is reached through a small legacy adapter. Nothing in the pane workspace is removed
or rewritten, and new workspace orchestration does not accumulate in `ChartManager`.

## Architecture guardrail before implementation

Spend one short session tracing `project_bootstrap`, the path implicated by #706 and the
minimum existing chart operations the showcase needs. Record:

- what must be ready before the workspace mounts
- what the workspace session owns
- the minimal operations of a legacy chart adapter: list, host, create, duplicate, delete,
  settings and resize notification
- what state should survive if a future project-data adapter replaces the active stores, for
  example when changing SpatialData samples

This is not an implementation of `datasources_v2`, sample switching or a `ChartManager`
replacement. It is a boundary that lets the showcase reuse today's charts without making
today's god-object the new architecture.

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
| `convert.ts` | Pane view to the internal workspace model, from the lab's `convertPaneView` |
| `WorkspaceSession.ts` | Owns the current view, layout, selection, dirty state and a persistence boundary |
| `LegacyChartAdapter.ts` | The minimum compatibility operations over existing chart machinery |
| `workspace.css` | The prototype's tokens and styles, scoped to the workspace root, with the dark variant |

### Changed, small hooks

| File | Change |
|---|---|
| `src/lib/browserSettings.ts` | Add `workspacePreviewSetting` (off/on), listed under Feature previews, like Convert legacy charts |
| `src/modules/project_bootstrap.tsx` | Select the workspace preview path early and construct its compatibility adapters |
| `src/charts/ChartManager.js` | Only the smallest hooks the legacy adapter proves necessary; no workspace state or persistence ownership |
| `src/charts/BaseChart.ts` | Settings opens in the inspector when the workspace is active, otherwise the floating dialog as now. Remove the pane drag and resize from hosted charts |
| Python view writers and routes | No change before the persistence checkpoint |

## Day by day

| Day | Work | Done when |
|---|---|---|
| 0 | Trace the integration boundary and record the workspace-session and legacy-adapter responsibilities | Implementation can start without assigning new workspace ownership to `ChartManager` |
| 1 | Preview setting. Mount the workspace with the lab's grid and cards through the adapter. In-memory pane-view conversion. Panes hidden | A real project opens with the preview on, showing every chart in one grid you can drag and resize |
| 2 | The frame: app bar and menus, status bar with view tabs, tool rail, light/dark tokens. Card header with maximise, lock, duplicate, delete | It looks like the prototype around real charts, in both themes |
| 3 | Inspector: docked settings, chart list, data sources. Toolbar on the selection | Selecting a chart edits it in the side panel. No floating dialog |
| 4 | Add chart, filter bar and responsive laptop/tablet widths. Run the layout checkpoint | Shared evidence supports a team choice of grid, tiling or hybrid |
| 5 | Run on the test projects (spatial, 23 charts, 3 data sources), fix what looks wrong and produce the demo build | The team can run it and decide whether to continue |

Durable save/reload is the next slice after the persistence checkpoint. It is required before
the preview can merge, but it does not block judging the initial workspace.

## Layout checkpoint

Keep the in-house grid as the starting implementation, then time-box a minimal tiling or hybrid
comparison before hardening layout-specific UI. Run the same Playwright scenarios against both
through roles and `data-chart-id`, never engine-specific selectors:

- resize an adjustable 2/3 arrangement
- move selected charts between regions
- exercise a dense view of small charts and measure unrelated-chart displacement
- enforce minimum usable chart dimensions and prevent thin slices
- avoid making individual charts taller than the viewport without a clear reason
- check desktop and narrower laptop/tablet widths
- host existing spatial/WebGL charts without geometry reloads

Playwright supplies repeatable geometry assertions, measurements and screenshots. People try
the same workflows and decide which feels better. The outcome can be a grid, tiling or a
hybrid; the seam lets the workspace survive that choice.

## Persistence checkpoint

Before adding durable save/load, agree:

- the canonical TypeScript/JavaScript schema for workspace state
- how Python models or writers are generated from, or validated against, that schema
- the versioned on-disk shape and whether `views_v2.json` is still the right isolation boundary
- conversion and migration rules

The workspace uses an internal model behind a repository boundary, so this decision does not
block the initial visual and interaction showcase. It does block a mergeable preview.

## What gets checked

- On the three lab projects, at 1440, 1024 and 768 px: no overlaps, no chart outside its card,
  no chart collapsed below the agreed minimum, and stable placement after gestures.
- At 390 px, the basic fallback does not overflow or expose unusable editing controls.
- The layout checkpoint's shared Playwright scenarios and screenshots run against each candidate.
- `pnpm exec tsgo`. Python route tests are added when persistence is implemented.
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
  - scrolling grid, tiling or a hybrid
  - the saved format
  - schema ownership and Python derivation
  - the longer-term project/data lifecycle
  - the inspector for charts with their own dialogs
  - how filters are shown
  - MDV-wide theming
  - undo
  - when (and whether) the preview merges, then becomes the default
- Only then: the PRD and build issues for the rest, if the verdict is yes.

## Decision backlog

These questions are deliberately kept here, next to the running showcase, rather than split
into ticket files. Settle them with the real UI in front of the team. Record durable decisions
as ADRs.

### Inspector and chart chrome

- Which settings can the docked inspector render from `getSettings()`?
- What happens for spatial, Viv, genome browser, table and gate charts with custom dialogs?
- What do popped-out charts use when no workspace inspector is present?
- What do the Chart, Charts and Data tabs own, including multi-selection?
- Where do plugin actions, annotations and the old per-data-source menu actions move?
- Which current actions belong in the card menu, toolbar, app menus or a later command palette?

### Tools and filters

- Which of Select, Lasso, Pan, Zoom, Polygon and Modify belong in the shared rail?
- How does the workspace show which charts support and respond to the active tool?
- What filter description contract must a chart provide: summary, data source, passing rows and
  clear action?
- How are linked-data-source and colour-legend filters represented?
- What does the filter bar show for charts that cannot yet describe their filters?

### Theme, undo and rollout

- What becomes the one token source across legacy CSS, Tailwind/shadcn and MUI?
- Where is light/dark stored, and which existing variables map to the new tokens?
- Does first-release undo cover only layout, or also chart creation, settings and filters?
- What evidence is required before the preview becomes the default?
- What evidence is required before the pane workspace and `views.json` can be retired?
- Can the static build consume the new workspace format before retirement?

Nothing merges to main until the team has tried the branch, the layout and persistence
checkpoints are resolved, and the plan has been iterated from that feedback.
