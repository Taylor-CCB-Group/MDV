Ticket: [02 Everything that reads or writes view layouts](../tickets/02-everything-that-reads-or-writes-view-layouts.md)
Map: [Workspace in MDV - Map](../README.md)
Source: MDV repo at `main` `b3c5a868` (read 2026-11-26). Paths are relative to the repo root. Line numbers are for that commit. Dev server counts are from `http://localhost:5055` on 2026-11-26.

## Summary

- A view is `{ dataSources: {<ds>: {layout, panelWidth, highlight}}, initialCharts: {<ds>: [chart config]}, links, viewImage, name }`. One split pane is built per key of `dataSources`. Each pane has its own menu bar, content div, layout menu and (if gridstack) its own GridStack instance.
- Layout is per pane, not per view. `layout` missing means `absolute` in the front end. Python writers mostly leave it out.
- Chart placement lives on each chart config: `size` and `position` in pixels (absolute), `gssize` and `gsposition` in grid units (gridstack). Both sets are often present on the same chart.
- On the dev server: 280 views in 37 projects. 74% of views have more than one pane. 70% of panes are effectively absolute. 29 views mix absolute and gridstack panes. Charts per view: min 0, median 7, max 23.

## 1. Data model as the code reads it

| Key | Where it lives | Values | Default when missing |
|---|---|---|---|
| `dataSources[ds].layout` | view | `"absolute"` or `"gridstack"` | `"absolute"` (`ChartManager.js:460`, `:479`; docs `docs/extradocs/views.md`) |
| `dataSources[ds].panelWidth` | view | percent, 0 to 100 | equal widths, if any pane lacks it (`ChartManager.js:676-684`) |
| `dataSources[ds].highlight` | view | row indices | none (`ChartManager.js:1168-1186`, `:1578-1585`) |
| `dataSources` itself | view | object | built from `initialCharts` keys, all absolute (`ChartManager.js:657-666`) |
| no view at all | view | none | every datasource gets a gridstack pane (`ChartManager.js:650-656`) |
| `size` | chart | `[w, h]` px | 300 x 300 (`GridstackManager.ts:494-501`) |
| `position` | chart | `[left, top]` px | first free cell from `findFreeSpace` (`GridstackManager.ts:503-514`) |
| `gssize` | chart | `[w, h]` grid units | gridstack measures the div (`GridstackManager.ts:178-185`) |
| `gsposition` | chart | `[x, y]` grid units | gridstack auto-places during view load (`GridstackManager.ts:187-191`, `:302-307`) |
| `location` `{x,y,width,height}` | chart | legacy grid cells | converted to `size`/`position` on add (`ChartManager.js:1613-1624`) |

Grid facts. GridStack runs with its default 12 columns. The code computes a column count at `GridstackManager.ts:52` but never passes it (line 68 is a comment). `cellHeight` is the content div height divided by `round(height / 150)` (`:51-55`). It is one value shared by all panes ("hack", line 55). So a grid unit's pixel height depends on the pane height when the first grid was made. `x` is clamped to 11 when converting from pixels (`:183`).

## 2. Front end readers and writers

| File:line | What it does | R/W |
|---|---|---|
| `src/charts/ChartManager.js:166-187` | Builds `this.dataSources` and `dsIndex`. The `DataSource` object later carries `contentDiv`, `menuBar`, `menuBarRoot`, `charts`. | W |
| `src/charts/ChartManager.js:212-229` | Main `menuBar` and main `contentDiv` (the area split into panes). | W |
| `src/charts/ChartManager.js:258` | Creates the single `GridStackManager`. | W |
| `src/charts/ChartManager.js:275-276` | `layoutMenus`: one "Change layout" menu per datasource. | W |
| `src/charts/ChartManager.js:443` | Interaction link adds an image chart only if `dsPanes[ods.name]` exists. | R |
| `src/charts/ChartManager.js:458-473` | `_setUpChangeLayoutMenu`: Absolute / Grid Stack choices; reads current `layout`. | R |
| `src/charts/ChartManager.js:476-503` | `changeLayout`: destroys gridstack or strips drag/resize, sets `view.layout`, re-applies layout to each chart that is not popped out. | R/W |
| `src/charts/ChartManager.js:589-634` | `_loadView`: loads first view, calls `_init`, then saves the view if it has no `viewImage`. | R/W |
| `src/charts/ChartManager.js:647-666` | `_init`: sets `viewData`; default gridstack panes when no view; legacy `initialCharts`-only views get empty `dataSources` entries. | R/W |
| `src/charts/ChartManager.js:668-673` | Clears every datasource's `contentDiv`, `menuBar`, `menuBarRoot`. | W |
| `src/charts/ChartManager.js:674-691` | Reads `panelWidth`s, calls `splitPane` with one pane per `dataSources` key, in object key order. | R |
| `src/charts/ChartManager.js:694-741` | Per pane: `dsPanes[name]`, `ds.charts = []`, `ds.menuBar` (702), `_setUpMenu` (721), `ds.contentDiv` (723), layout menu (738). | W |
| `src/charts/ChartManager.js:777-823` | Loads every chart in `initialCharts`, with `_inInit` true (794-827). | R |
| `src/charts/ChartManager.js:1122-1135` | `getState`: writes `panelWidth` back by parsing Split.js's `calc(N% ...)` width string. | R/W |
| `src/charts/ChartManager.js:1146-1165` | `getState`: for absolute panes sets `config.position` from `div.offsetLeft/Top`; gridstack panes rely on the observers. | W |
| `src/charts/ChartManager.js:1167-1186` | `getState`: deep-copies `viewData`, sets `initialCharts`, writes `highlight`. | W |
| `src/charts/ChartManager.js:1223-1238` | `addMenuIcon`: adds icons to `dsIndex[ds].menuBar` (used by gridstack "compact"). | W |
| `src/charts/ChartManager.js:1485-1501` | `_setUpMenu`: mounts React `DataSourceMenuBar` in the pane's `menuBar` with `flushSync`, so gridstack measures the right height. | W |
| `src/charts/ChartManager.js:1595-1650` | `addChart`: converts legacy `location`, calls `positionChart`, creates the absolute div in `ds.contentDiv`. | R/W |
| `src/charts/ChartManager.js:2039` | `_addChart` hands the new chart to `_makeChartRD`. | W |
| `src/charts/ChartManager.js:2245-2263` | `setChartsAsGrid`: moves all chart divs into a pixel grid. No caller found in `src`. | W |
| `src/charts/ChartManager.js:2283-2285` | `_popOutChart` calls `popoutChart`. | W |
| `src/charts/ChartManager.js:2291-2314` | `_makeChartRD`: gridstack panes go to `gridStack.manageChart(chart, ds, _inInit)`; others get `makeDraggable`/`makeResizable`, `onresizeend` calls `setSize`. | R/W |
| `src/charts/GridstackManager.ts:47-86` | `getGrid`: one GridStack per datasource on `ds.contentDiv`; adds "compact layout" icon to the pane menu. | W |
| `src/charts/GridstackManager.ts:90-107` | `compact`: writes `gssize`/`gsposition` for every chart in the pane. | W |
| `src/charts/GridstackManager.ts:109-152` | `destroy`: converts gridstack boxes to absolute px, writes `size`/`position`, clears `gssize`/`gsposition`. | W |
| `src/charts/GridstackManager.ts:164-191` | `manageChart`: derives w/h/x/y from the div, or uses saved `gsposition`/`gssize`. | R |
| `src/charts/GridstackManager.ts:233-287` | ResizeObserver and MutationObserver write `gssize`/`gsposition` after any move or resize (debounced 20 ms). | W |
| `src/charts/GridstackManager.ts:309-329` | Per-chart lock button. | W |
| `src/charts/GridstackManager.ts:342-344`, `383-418` | `gridstackPopoutCallback` and `removeLayout`: unwrap the chart from gridstack. | W |
| `src/charts/GridstackManager.ts:428` | `Config` type: `size`, `position`, `gssize`, `gsposition`. | R |
| `src/charts/GridstackManager.ts:429-490` | `findFreeSpace`: scans visible chart rects in a pane for a free 2 x 1 cell. | R |
| `src/charts/GridstackManager.ts:492-525` | `positionChart`: initial px box from `size`/`position`, or from `gssize`/`gsposition` if the pane is gridstack. | R |
| `src/charts/ViewManager.ts:13-27` | `ViewDataDataSource` and `View` types (`layout?: string`). | R |
| `src/charts/ViewManager.ts:141-178` | `convertLegacyChart`: re-adds the chart with its own config, so placement fields carry over. | R/W |
| `src/charts/ViewManager.ts:238-271` | `changeView`: destroys gridstack panes, removes charts, empties `contentDiv`, calls `_init`. | R/W |
| `src/charts/ViewManager.ts:274-287` | `createImageofView`: screenshots `cm.contentDiv` (all panes) for `viewImage`. | R |
| `src/charts/ViewManager.ts:305-333` | `saveView`: `getState` plus `viewImage`, fires `state_saved`. | W |
| `src/charts/ViewManager.ts:337-394` | `addView`: new view gets `dataSources[name] = {}` (no layout, so absolute). Clone copies current state. | W |
| `src/charts/ViewManager.ts:399-426` | `deleteView`: sends `view: null`. | W |
| `src/charts/ViewManager.ts:493-538` | `hasUnsavedChanges`: deep-compares the whole view, so any placement change counts. | R |
| `src/charts/BaseChart.ts:280-283` | Duplicate chart: copies the full config, including `gsposition`/`gssize`/`position`. | W |
| `src/charts/BaseChart.ts:289-357` | Full screen: restores `oldSize` on exit; re-runs `manageChart(..., remanage=true)` for gridstack panes (325-326). | R/W |
| `src/charts/BaseChart.ts:875` | `removeLayout?()` hook declared. | R |
| `src/charts/BaseChart.ts:1231-1244` | `getConfig`: no layout logic. Comment notes the `size`/`position` issue outside absolute. | R |
| `src/charts/BaseChart.ts:1273-1292` | `setSize` / `_setDimensions`: writes `config.size` from the div rect. | W |
| `src/utilities/Popout.ts:4-123` | Pop-out: calls `gridstackPopoutCallback` (15), strips drag/resize (20-21); on close re-appends to the pane and calls `_makeChartRD` (109-112). | R/W |
| `src/utilities/Elements.js:34-50` | `splitPane`: Split.js, horizontal only, gutter 5 px, `sizes` from `panelWidth`. | W |
| `src/utilities/Elements.js:237`, `354`, `368`, `428` | `makeResizable`, `removeResizable`, `removeDraggable`, `makeDraggable` for absolute panes. | W |
| `src/charts/schemas/ChartConfigSchema.ts:50-61` | `GridStackConfigSchema`: `gssize` and `gsposition` optional int tuples. | R |
| `src/charts/schemas/ChartConfigSchema.ts:74-87` | `BaseConfigSchema`: `size` required. `position` is not in the schema. | R |
| `src/charts/schemas/ChartConfigSchema.ts:631-656` | `ViewConfigSchema`: `layout` enum, `panelWidth` 0-100. No `highlight`, `links`, `viewImage`. | R |
| `src/charts/schemas/ChartConfigSchema.ts:668` | `ChartManagerConfigSchema.gridstack` boolean. No reader found. | R |
| `src/charts/chartConfigUtils.ts:165-173` | Validates each chart against its schema; failures are only logged. | R |
| `src/charts/legacyCharts.ts:43-51` | `chartSize`: reads `size` (or `location`) to scale point radius when converting a classic scatter. | R |
| `src/charts/charts.d.ts:168-180` | `DataSource` type: `contentDiv`, `menuBar`, `menuBarRoot`, `charts`. | R |
| `src/react/components/DataSourceMenuBar.tsx:94`, `:111` | Pane full screen uses `menuBar.parentElement` as the pane. | R |
| `src/react/components/DataSourceMenuBar.tsx:123-127` | "Change layout" icon opens `cm.layoutMenus[name]`. | R |
| `src/charts/dialogs/AddViewDialog.tsx:21-33` | Pre-ticks datasources that are keys of `viewData.dataSources`. | R |
| `src/charts/dialogs/AddChartDialogReact.tsx:402` | Adds a chart with no placement; `positionChart` finds space. | W |
| `src/charts/columnRemovalUtils.ts:150-151` | Loads every view and reads `initialCharts[ds]`. | R |
| `src/react/chartLinkHooks.ts:20-28` | Reads `viewData.links`. | R |
| `src/dataloaders/DataLoaderUtil.ts:103-155` | `viewLoader`: static build reads from `views.json`; server posts to `get_view`. | R |
| `src/modules/projectRuntime.ts:46-60` | Fetches all of `views.json` at start-up and passes it to the loader. | R |
| `src/modules/project_bootstrap.tsx:263-284` | On `state_saved`, mounts `ProjectStateHandler`. | W |
| `src/react/ProjectStateHandler.tsx:22-26` | Posts the state to `save_state` (skipped for static). | W |
| `src/charts/dialogs/ChatAPI.ts:682`, `:816` | After chat creates a view, adds its name and switches to it. | R |
| `src/legacy/hyperion.js:342-447`, `1290`, `1352-1357` | Legacy: extra `viewData.spatial_stats` key; own `views.json` load and `save_state`. No other file references it. | R/W |
| `src/modules/desktop_index.js:45-54`, `93-94`, `121` | Desktop entry: own view load and save. Out of scope per map. | R/W |

## 3. Python readers and writers

| File:line | What it does | Layout written |
|---|---|---|
| `python/mdvtools/mdvproject.py:318-321` | New project: `views.json` = `{}`. | none |
| `python/mdvtools/mdvproject.py:360-366` | `views` property: whole-file read and write. | n/a |
| `python/mdvtools/mdvproject.py:62-74`, `88-108` | Catalog thumbnail: first `viewImage` in `all_views` order. | R |
| `python/mdvtools/mdvproject.py:525-560` | Finds charts that reference a column, across all views' `initialCharts`. | R |
| `python/mdvtools/mdvproject.py:877-880` | Interactions: view with `initialCharts` only. | none (absolute) |
| `python/mdvtools/mdvproject.py:1179-1186` | Image datasource update: view from `image_view_prototype`. | gridstack |
| `python/mdvtools/mdvproject.py:1556-1572` | Delete datasource: deletes every view whose `initialCharts[name]` is non-empty. | R/W |
| `python/mdvtools/mdvproject.py:1884-1911`, `2067-2092` | `add_datasource` / `_polars`: table at `size [792,472]`, `position [10,10]`; new view is `initialCharts` only. | none (absolute) |
| `python/mdvtools/mdvproject.py:2110-2139` | `create_view_with_all_datasources`: one pane per datasource, `layout: gridstack`, one table chart with px `size`/`position`, no `gs*`. | gridstack |
| `python/mdvtools/mdvproject.py:2141-2154` | `create_table_plot` / `convert_plot_to_json`. | n/a |
| `python/mdvtools/mdvproject.py:2435-2449` | `save_state`: `set_view(currentView, view)`. | passes through |
| `python/mdvtools/mdvproject.py:2510-2565` | `get_view`, `set_view`: requires `initialCharts` (KeyError if missing), adds chart ids, `null` deletes. | passes through |
| `python/mdvtools/mdvproject.py:2567-2595`, `2597-2610` | `rename_view` (also sets `view["name"]`), `reorder_views`. | n/a |
| `python/mdvtools/mdvproject.py:2964-3003` | `get_centroid_plot`: sets `size` from ROI x scale. | px |
| `python/mdvtools/server.py:108-150` | `classify_view_change` reads `currentView` and `view` of a save. | R |
| `python/mdvtools/server.py:322-353` | Serves `views.json` raw (GET). | R |
| `python/mdvtools/server.py:397-409` | `POST get_view`; records a `view_open` usage event. | R |
| `python/mdvtools/server.py:459-478` | `POST save_state`. | W |
| `python/mdvtools/serverlite.py:120-125`, `183-188` | Same two routes, lite server. | R/W |
| `python/mdvtools/project_protocols.py:32` | `get_view` protocol. | R |
| `python/mdvtools/charts/view.py:13-25` | `Panel` and `View` TypedDicts: `layout` Literal, `panelWidth`, required `initialCharts`. | type |
| `python/mdvtools/charts/base_plot.py:6-15` | Every Python chart class writes px `size` and `position`. No `gs*`. | px |
| `python/mdvtools/charts/base_plot.py:32-33` | `set_position` sets `self.position`, not `plot_data["position"]`, so it has no effect on output. | n/a |
| `python/mdvtools/template_views/create_default_views.py:92-111` | `_chart`: every chart gets `size [400,300]`, `position [0,0]` and real `gsposition`/`gssize`. | both |
| `python/mdvtools/template_views/create_default_views.py:270-295` | `Grid`: left-to-right 12-column packer. | gridstack |
| `python/mdvtools/template_views/create_default_views.py:1203`, `1263`, `1368-1389`, `1483`, `1536-1554` | Five views. Panes are gridstack; widths 100, or 70/30, or 70/15/15. | gridstack |
| `python/mdvtools/template_views/create_default_views.py:1605-1665` | Writes `views.json` and `state.json` directly, not through `set_view`. Backs up `views.json` first. | W |
| `python/mdvtools/chart_prototypes.py:13-235` | Prototype charts with `size`, `position`, `gsposition`, `gssize`. No importer anywhere in the repo outside `static`. | both |
| `python/mdvtools/image_view_prototype.py:7-94` | One gridstack pane, `panelWidth 100`, one chart with `size`, `gsposition [0,0]`, `gssize [12,5]`. | gridstack |
| `python/mdvtools/spatial/spatial_view_template.json:48-131` | `cells` gridstack 69.08%, `genes` absolute 30.92%. Genes charts carry both px and `gs*`. | mixed |
| `python/mdvtools/spatial/project_helpers.py:29-65` | Loads that template, swaps `"cells"`/`"genes"` names by string replace, `set_view`. | mixed |
| `python/mdvtools/spatial/conversion.py:1725-1743` | Calls it once per spatial group. | mixed |
| `python/mdvtools/spatial.py:1-46`, `48-104`, `107-165` | Centroid, interaction-matrix and ROI views: px `position`/`size`; `panelWidth` 20/80; no `layout`. | absolute |
| `python/mdvtools/network_helpers.py:268-343` | Ligand network view: px placement, `panelWidth`, no `layout`. | absolute |
| `python/mdvtools/conversions.py:359` | AnnData/MuData default view: `initialCharts` only. | absolute |
| `python/mdvtools/conversions.py:362-383` | Builds `new_views` with two 50% panes added to every view, but never saves them. | none |
| `python/mdvtools/conversions.py:1040-1063`, `1210-1234` | Regulamentary and CapSeqUM projects from `templates/views/*.json`: px only. | absolute |
| `python/mdvtools/conversions.py:1822-1850` | Project merge: re-keys `initialCharts` and `dataSources` by datasource rename; keeps layout. | passes through |
| `python/mdvtools/spatial/xenium.py:250-280` | Default view `initialCharts` only; or adds two 50% panes to every view and saves. | absolute |
| `python/mdvtools/csv_to_static.py:208-209`, `261-262` | Static export view: `initialCharts` only. | absolute |
| `python/mdvtools/templates/views/capsequm.json`, `regulamentary.json` | Legacy chart lists with px `position`/`size`. | absolute |
| `python/mdvtools/llm/column_field_resolve.py:376-401` | `ensure_view_gridstack_layout`: forces `layout: gridstack` on every pane with charts; overrides `absolute`. | gridstack |
| `python/mdvtools/llm/langchain_mdv.py:1246-1276` | After the generated script runs: normalise params, prune, force gridstack, `set_view`. | gridstack |
| `python/mdvtools/llm/guidance.py:557-591` | Prepends a summary TextBox at `size [900,320]`, `position [10,10]`, no `gs*`. | px only |
| `python/mdvtools/llm/templates.py:426` | Prompt tells the model to set `initialCharts` only; ChatMDV applies gridstack. | n/a |
| `python/mdvtools/llm/local_files_utils.py:8`, `langchain_mdv.py:415` | Feeds `test_projects/RAG_examples/ANNDATA_examples` (29 files) to the model as examples. | n/a |
| `python/mdvtools/test_projects/RAG_examples/**` | 194 files call `set_view` with `initialCharts` only and px `position` (88 use `[10, 10]`). None use `gs*` or `layout`. | absolute |
| `python/mdvtools/llm/chat_cli.py:106-110`, `262-279` | Counts charts; renames the generated view by copy then delete. | R/W |
| `python/mdvtools/llm/verification.py:115`, `189`, `285-291` | Reads `initialCharts` to check charts. | R |
| `python/mdvtools/llm/code_manipulation.py:318` | Checks a parsed view name exists. | R |
| `python/mdvtools/examples/flexible_network_example.py:81-153`, `ligand_network_example.py:177-190`, `rename_delete_columns_example.py:120` | Example scripts writing views. | absolute / gridstack |

Note: `python/mdvtools/static/` was ignored as built output. It also contains `load_network_examples.py`.

## 4. Tests and fixtures that pin layout

| File:line | What it pins |
|---|---|
| `python/mdvtools/tests/test_column_field_resolve.py:300-319` | `ensure_view_gridstack_layout` adds missing panes as gridstack and overrides absolute; keeps `panelWidth`. |
| `python/mdvtools/tests/test_spatial_view_template.py:24-46` | Spatial default view has the SpatialData chart and panes `obs`, `var`. Does not check layout. |
| `python/mdvtools/tests/test_template_views.py:166-184` | Default views get the expected names. |
| `python/mdvtools/tests/test_soft_delete_column.py:212-219` | `create_view_with_all_datasources` table columns. |
| `python/mdvtools/tests/generate_synthetic_spatial_project.py:349-503` | Fixture project: gridstack panes 100 / 0, charts with `gsposition`/`gssize`. |
| `python/mdvtools/tests/generate_dotplot_label_project.py:147` | Fixture view: `initialCharts` only. |
| `src/tests/react/DataSourceMenuBar.test.tsx:9`, `:111` | Layout icon calls `layoutMenus[name].show`. |
| `tests_playwright/project/charts/density_grid_selection.spec.ts:26`, `59-60` | Adds a chart with `size [700,700]`, `position [10,10]`; notes `positionChart` throws before the view loads. |
| `tests_playwright/project/charts/scatter_selection_tools.spec.ts:34`, `61-62` | Same pattern. |
| `uploads/views.json` | One view, `initialCharts` only, two empty panes. |
| `examples/basic_example.js:85-101`, `examples/10Mand1_points.js:45` | `only_view` with px `size`/`position`. |

No unit test covers `GridstackManager`, `changeLayout`, `splitPane` or `getState` placement.

## 5. How saved views use layouts

### 5a. Repo JSON

| Source | Views | Panes | Layout per pane | Charts |
|---|---|---|---|---|
| `spatial/spatial_view_template.json` | 1 | 2 | `cells` gridstack, `genes` absolute (mixed) | 3 (1 + 2) |
| `uploads/views.json` | 1 | 2 (from `initialCharts`) | missing, so absolute | 0 |
| `templates/views/regulamentary.json` (+ browser) | 1 | 1 | missing, so absolute | 6 |
| `templates/views/capsequm.json` (+ browser) | 1 | 2 | missing, so absolute | 6 (4 + 2) |
| `image_view_prototype.py` | 1 | 1 | gridstack | 1 |

Code-built views in the repo: `create_default_views` and ChatMDV write gridstack. `create_view_with_all_datasources` writes gridstack. All other Python writers leave `layout` out, so the front end shows them absolute.

### 5b. Dev server

Method: `GET /projects`, then `GET /project/<id>/views.json` and `/state.json` for each of the 46 projects. This is the raw file route (`server.py:322`). It avoids `POST get_view`, which records a usage event. Pane layout is "effective": missing counts as absolute.

| Measure | All views | Distinct views (copies merged, `viewImage` ignored) |
|---|---|---|
| Projects | 46 (37 with at least one view; 35 distinct names) | |
| Views | 280 | 197 |
| Views with 0 / 1 / 2 / 3 / 4 / 5 panes | 2 / 72 / 162 / 37 / 1 / 6 | 1 / 57 / 112 / 23 / 1 / 3 |
| Views with more than one pane | 206 (74%) | 139 (71%) |
| Views with charts in more than one pane | 181 | 120 |
| Views all absolute / all gridstack / mixed | 175 / 74 / 29 | 136 / 46 / 14 |
| Panes absolute / gridstack | 379 / 162 | 277 / 92 |
| Panes with `layout` written: gridstack / absolute / missing | 162 / 88 / 291 | |
| Panes holding charts: absolute / gridstack | 338 / 149 | 253 / 80 |
| Charts per view, min / median / mean / max | 0 / 7 / 6.7 / 23 | 0 / 7 / 6.2 / 23 |
| Charts per non-empty pane, min / median / max | 1 / 2 / 21 | 1 / 2 / 21 |
| Empty views (no charts) | 18 | 11 |
| Views with no `dataSources` key | 9 | |
| Views with `panelWidth` on every pane | 258 | |
| Views whose `panelWidth`s do not sum to 100 | 5 (all in one project, sum 71.43) | |
| Views in `views.json` but not in `all_views` | 6 | |

Chart count histogram (all 280 views): 0:18, 1:12, 2:27, 3:24, 4:23, 5:12, 6:10, 7:28, 8:29, 9:20, 10:14, 11:43, 12:6, 14:2, 15:5, 16:5, 17:1, 23:1. Largest: one view with 23 charts.

Placement fields on the 1,866 charts:

| Pane layout | `gs*` | `position` | `size` | Charts |
|---|---|---|---|---|
| absolute | no | yes | yes | 1,122 |
| absolute | yes | yes | yes | 23 |
| gridstack | yes | no | yes | 460 |
| gridstack | yes | yes | yes | 239 |
| gridstack | no | yes | yes | 22 |

Every chart has `size`. No chart uses legacy `location`.

Top-level view keys seen: `initialCharts` 280, `dataSources` 271, `viewImage` 209, `name` 74, `links` 70.

Mixed views (29, 14 distinct). The pattern is a gridstack cells or spatial pane next to an absolute genes, protein or RNA pane. Most come from the spatial template (`cells` gridstack, `genes` absolute) and from spatial atlas projects that pair a gridstack cells pane with absolute gene panes. Project names are kept out of this public copy.

## 6. Contracts a new layout format must keep

Who reads each field today, so who breaks if it changes.

| Field | Readers that depend on it |
|---|---|
| `dataSources` keys (pane list and order) | `ChartManager._init` builds one pane per key, in key order (`:674-741`). `AddViewDialog` (`:23`). `getState` highlight write (`:1168-1186`). Interaction links check `dsPanes` (`:443`). Python merge re-keys it (`conversions.py:1835-1840`). `xenium.py` and `conversions.py` add panes by key. |
| `dataSources` may be missing | `_init` rebuilds it from `initialCharts` (`:659-664`). 9 server views and most Python writers rely on this. |
| `dataSources[ds].layout`, missing = absolute | `_init`, `changeLayout`, `_makeChartRD`, `getState`, `positionChart`, `BaseChart` full screen, `ViewManager.changeView`/`addView`. Python: `ensure_view_gridstack_layout`, `create_default_views`, `create_view_with_all_datasources`, `image_view_prototype`, spatial template, `charts/view.py` Literal, `ViewConfigSchema` enum. 291 server panes have no `layout`. |
| `dataSources[ds].panelWidth` (percent) | `splitPane` sizes (`:676-691`); `getState` parses it back from Split.js `calc(...)` (`:1131-1133`). Written by spatial, network, xenium, conversions, default views, template. |
| `dataSources[ds].highlight` | `_restoreHighlights` (`:1578-1585`), `getState` (`:1180-1186`). Not in `ViewConfigSchema`. |
| `initialCharts[ds]` as an array of chart configs keyed by datasource | Every Python writer; `set_view` requires the key and walks it for ids (`mdvproject.py:2525-2530`); `delete_datasource` (`:1569-1571`); column reference scan (`:556`); `columnRemovalUtils.ts:151`; ChatMDV normalise, prune, verification, guidance, chat_cli; `_init`. A chart's datasource is the key it sits under. |
| `size` `[w,h]` px | Required by `BaseConfigSchema`. `positionChart`, `BaseChart.setSize`/`_setDimensions`, full screen `oldSize`, `legacyCharts.chartSize` (radius scaling), `GridstackManager.destroy`. Every Python chart class writes it. All 1,866 server charts have it. |
| `position` `[l,t]` px | `positionChart`, `getState` (absolute panes), `GridstackManager.destroy`. Written by `BasePlot`, spatial, network, conversions, RAG examples, Playwright tests. Not in the zod schema. |
| `gssize` `[w,h]`, `gsposition` `[x,y]` grid units, 12 columns | `GridstackManager.manageChart`, `compact`, observers, `positionChart`. Written by `create_default_views`, spatial template, `image_view_prototype`, synthetic test fixture. 722 server charts carry them. Pixel meaning depends on pane height via shared `cellHeight`. |
| Gridstack charts without `gsposition` get auto-placed | `manageChart(..., autoPosition = _inInit)` (`ChartManager.js:2295`, `GridstackManager.ts:302-307`). ChatMDV views and `create_view_with_all_datasources` depend on this (22 server charts). |
| Legacy `location {x,y,width,height}` | Converted on add (`ChartManager.js:1613-1624`), `legacyCharts.ts:46`. The UCSC browser uses `config.location` for a genome range (`UCSCBrowserComponent.tsx:49-91`), so the name is shared. |
| `links`, `viewImage`, `name`, `spatial_stats` (legacy) | `getState` deep-copies `viewData`, so unknown keys survive a save (`:1167`). `viewImage` feeds the gallery (`ViewManager.ts:289-297`) and catalog thumbnail (`mdvproject.py:62-74`). `rename_view` writes `name`. |
| Per-pane DOM: `ds.contentDiv`, `ds.menuBar`, `dsPanes` | `GridstackManager` (`:49`, `:119`, `:438`, `:447`), `addMenuIcon` (`:1225-1226`), `addChart` (`:1650`), `DataSourceMenuBar` pane full screen (`:94`), `getState` pane width, `ViewManager.createImageofView` (whole `contentDiv`). A datasource with no `contentDiv` is skipped by `getState` (`:1129`). |
| `layoutMenus[ds]` | `DataSourceMenuBar.tsx:126`, unit test `DataSourceMenuBar.test.tsx:111`. |
| `chart.removeLayout`, `div.gridstackPopoutCallback` | `Popout.ts:15`, `BaseChart.ts:875`, `GridstackManager.destroy` (`:145`). |
| Save path | `getState` -> `state_saved` -> `ProjectStateHandler` -> `POST save_state` -> `MDVProject.save_state` -> `set_view` writes `views.json` as sent. Static builds read `views.json` and never save. |
| Unsaved-change detection | `hasUnsavedChanges` compares the whole view, so any field that changes on load (for example `size` from `setSize`) reads as a change (`GridstackManager.ts:288-295` comment). |

## Could not verify

- Behaviour in a browser. Nothing was run in the app. All front end facts are from reading code.
- Whether GridStack's 12-column default is what users see. It depends on the installed `gridstack` version's defaults; the code passes no `column`.
- Dev server projects with the same name are counted separately in "All views". The "Distinct views" column merges identical view content only.
- Views on other servers or in user folders outside this dev server.
