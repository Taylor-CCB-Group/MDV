Type: research
Status: resolved
Blocked by: none
Map: [Workspace in MDV - Map](../README.md)

## Question

For every chart type, what is in its title bar and menus today, and where does it come from?

- `BaseChart` icons: reset filter, more, settings, duplicate, info, full screen, unpin, pop out, close.
- Per-type `addMenuIcon` calls, and callers of `ChartManager.addMenuIcon` (`GridstackManager`, `src/legacy/hyperion.js`, others).
- Which charts' settings are fully described by `getSettings()` GuiSpecs, and which have custom panels (for example the spatial `ImageLayerPanel`, MobX control islands).
- Which charts own selection tools (`SelectionOverlay` modes: Pan, Rectangle, Polygon, Lasso, Modify; spatial tools).
- For each chart type: how its filter is stored and cleared (`removeFilter`, `resetButton`), and whether there's enough state to describe the filter in words.

## Context

Feeds [09 Extension points for plugins](merged/09-extension-points-for-plugins.md), [10 What the inspector shows](merged/10-what-the-inspector-shows.md), [11 Global tool rail vs chart-owned tools](merged/11-global-tool-rail-vs-chart-owned-tools.md) and [12 How charts describe their filters](merged/12-how-charts-describe-their-filters.md).

## Answer

Every chart gets the same BaseChart title bar: reset filter (hidden until a filter is set), more, settings, duplicate, info and full screen. ChartManager then adds pop out and close, gridstack adds a lock icon, and pinning adds an unpin icon. Only six types add their own icons: both tables, genome_browser, viv_scatter_plot, VivMdvRegionReact and SpatialDataMdvRegionReact.
Most of the 51 registered types are fully described by getSettings() GuiSpecs. The custom UIs are the spatial layer dialog (ImageLayerPanel and the other layer panels), the Viv channel dialogs, genome browser tracks, the table column dialogs, gate dialogs, and controls inside the chart body (selection_dialog, geneNetwork).
SelectionOverlay (Pan, Rectangle, Polygon, Freehand, Modify) is used only by the Deck scatter family, VivMdvRegionReact and SpatialDataMdvRegionReact. Classic charts draw their own brushes or filter on click.
There are about a dozen filter mechanisms. Category, range, sankey, category-pair, multi-line brush, search, selection-dialog and colour-legend filters hold enough state to put in words. Polygon, index, dot-plot and genome-browser filters do not, or only partly. Only some filters are saved in the config.
ChartManager.addMenuIcon is used live only by GridstackManager. hyperion.js is not imported anywhere, and its "_main" target has no element.

Findings: [03 Chart chrome inventory - findings](../research/03-chart-chrome-inventory-findings.md)
