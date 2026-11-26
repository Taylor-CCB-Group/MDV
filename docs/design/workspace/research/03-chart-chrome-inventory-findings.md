Ticket: [03 Chart chrome inventory](../tickets/03-chart-chrome-inventory.md)
Source: MDV repo at `main` (b3c5a868), read 2026-11-26. All paths are relative to the repo root. Line numbers are for that commit.

## 1. Title bar and BaseChart icons

`BaseChart` builds the title bar in its constructor. The bar is a `div.ciview-chart-titlebar` holding a title `div.ciview-chart-title` and a `div.ciview-chart-menuspace` (src/charts/BaseChart.ts:177-196). `config.title_color` sets the bar background (BaseChart.ts:202-204). Every icon goes through `BaseChart.addMenuIcon`, which appends a `<span>` with an `<i class="ciview-chart-icon ...">` to `menuSpace`, attaches a tippy tooltip, and returns the span (BaseChart.ts:493-530). Icons appear left to right in the order they are added.

| Order | Icon | Class / tooltip | Added at | When it shows | What it does |
|---|---|---|---|---|---|
| 1 | Reset filter | `fas fa-sync`, "remove filter" | BaseChart.ts:207-214 | Hidden at start. Shown by each chart when it sets a filter, or by `updateResetButtonVisibility` when `colorLegendFilter` or `getFilter()` is set (BaseChart.ts:663-669) | Calls `removeFilter()`, then `clearColorLegendFilter`, then hides itself |
| 2 | More | `fas fa-bars`, "more" | BaseChart.ts:274-275 | Always | Opens the context menu (see section 2) |
| 3 | Settings | `fas fa-cog`, "settings" | BaseChart.ts:277 | Always | Opens one `SettingsDialogReactWrapper` per chart (BaseChart.ts:1105-1111) |
| 4 | Duplicate | `fas fa-copy`, "duplicate chart" | BaseChart.ts:279-284 | Always | Deep-copies `getConfig()`, new id, title + " copy", `chartManager.addChart` |
| 5 | Info | `fas fa-info`, tooltip is `config.legend` or "No description" | BaseChart.ts:287 | Always | Tooltip only, no click handler. The "Chart Legend" setting updates its `aria-label` (BaseChart.ts:953-956) |
| 6 | Full screen | `fas fa-expand`, "Full Screen" | BaseChart.ts:347-360 | Always | `requestFullscreen` on the chart div. A document `fullscreenchange` listener swaps the icon to `fa-compress` and back (BaseChart.ts:290-346) |
| 7 | Chart-specific icons | various | subclass constructors | Per type (see section 3) | |
| 8 | Lock position | `fas fa-unlock` / `fa-lock`, "lock position" | src/charts/GridstackManager.ts:317-329 | Only in a gridstack layout. Added by `manageChart`, which `_makeChartRD` calls before the popout and close icons (src/charts/ChartManager.js:2039, 2291-2296) | Toggles gridstack `locked` and the `gridLock` class. Removed by `removeLayout` (GridstackManager.ts:396) |
| 9 | Pop out | `fas fa-external-link-alt`, "popout" | ChartManager.js:2041-2045 | Always. Nothing hides it after pop-out (src/utilities/Popout.ts has no reference to it) | `popoutChart(chart)` opens a new window (ChartManager.js:2283-2285) |
| 10 | Close | `fas fa-times`, "remove chart" | ChartManager.js:2046-2052 | Always | `chart.remove()`, removes the div, deletes from `charts`, removes links, fires `chart_removed` |
| 11 | Unpin | `fas fa-thumbtack`, "unpin chart" | BaseChart.ts:1097-1103 | Added only after "pin chart" is chosen from the More menu | Removes itself and calls `unpinChart()` |

Other title-bar facts:

| Fact | Source |
|---|---|
| In absolute layout the title div is the drag handle | ChartManager.js:2304-2308 |
| In gridstack layout the drag handle is `.ciview-chart-title` too | GridstackManager.ts:40, 60 |
| The More icon can carry a coloured pip with extra tooltip text. Only the SpatialData chart sets one | BaseChart.ts:534-543; src/react/components/SpatialDataMDVReact.tsx:382-386 |
| A chart that fails to build gets an error card with its own close icon (`createMenuIcon`, not `addMenuIcon`) | ChartManager.js:1751-1770 |
| Pinned state lives in `chart.isPinned` only. It is not saved in config, so the unpin icon never appears on load | BaseChart.ts:1112, 1137-1146 |

## 2. More (context) menu items

The menu is built on each open: `getContextMenu()` then three items that BaseChart always appends (BaseChart.ts:247-272, 1120-1156).

| Item | Condition | Source |
|---|---|---|
| create svg image / create png image | chart has `getImage` | BaseChart.ts:1122-1135 |
| pin chart | chart has `pinChart` and is not pinned | BaseChart.ts:1137-1146 |
| download data | chart has `getChartData` | BaseChart.ts:1148-1154 |
| items from `addToContextMenu()` | chart (or outside code) sets it | BaseChart.ts:1155-1157 |
| convert to ... | `getLegacyConversionTarget` finds a replacement | BaseChart.ts:250-256 |
| debug / report chart | always. Carries the menu pip | BaseChart.ts:257-265 |
| copy config JSON to clipboard | always | BaseChart.ts:266-270 |

Who provides the optional hooks:

| Hook | Charts | Source |
|---|---|---|
| `getImage` | all `SVGChart` subclasses (this includes the WGL family), genome_browser, Deck scatter family, compound/split dot plot, multi_line_chart_react | src/charts/SVGChart.js:508; src/charts/GenomeBrowser.js:425; src/react/components/DeckScatterReactWrapper.tsx:144; src/react/components/dotplot/DotPlotReactWrapper.tsx:228; src/react/components/multiline/MultiLineReactWrapper.tsx:214 |
| `pinChart` | row_chart, wordcloud, ring_chart (CategoryChart); bar_chart; sankey_chart; table_chart; WGLChart subclasses (wgl_scatter_plot, density_scatter_plot, image_scatter_plot, viv_scatter_plot, wgl_3d_scatter_plot, violin_plot, box_plot) | src/charts/CategoryChart.js:84; HistogramChart.js:187; SankeyChart.js:225; TableChart.js:365; WGLChart.js:216 |
| `getChartData` | heat_map only | src/charts/HeatMap.js:69 |
| `getContextMenu` override | density_scatter_plot / image_scatter_plot add "Set as Default Image" when the store has regions and no viv. wgl_scatter_plot and subclasses add "Adjust <group> position" (opens `OffsetDialog`) when the store has matching `offsets` | src/charts/DensityScatterPlot.js:146-162; src/charts/WGLScatterPlot.js:615-637 |
| `addToContextMenu` | set from outside by `src/legacy/hyperion.js` ("Set as Default Image") | src/legacy/hyperion.js:1386-1399 |

## 3. Registered chart types and their extra title-bar items

Registrations are side-effect imports in src/charts/registerChartModules.ts. `table_chart_react` is also lazily loaded (src/charts/lazyChartRegistrations.ts:1-3). "Flags" are from the registration object. "Extra chrome" means anything beyond section 1 and 2.

| Type key | Name | Class | Registered at | Flags | Extra chrome |
|---|---|---|---|---|---|
| bar_chart | Histogram | HistogramChart | src/charts/HistogramChart.js:237 | | none |
| row_chart | Row Chart | RowChart | src/charts/RowChart.js:242 | | none |
| wordcloud | Word Cloud | RowChart | RowChart.js:253 | | none |
| ring_chart | Pie Chart | RingChart | src/charts/RingChart.js:87 | | none |
| stacked_row_chart | Stacked Row Chart | StackedRowChart | src/charts/StackedRowChart.js:171 | | none |
| heat_map | Heat Map | HeatMap | src/charts/HeatMap.js:413 | | none |
| dot_plot | Dot Plot | DotPlot | src/charts/DotPlot.js:714 | | none (fraction legend in content) |
| violin_plot | Violin Plot | ViolinPlot | src/charts/ViolinPlot.js:427 | | none |
| box_plot | Box Plot | BoxPlot | src/charts/BoxPlot.js:99 | | none |
| sankey_chart | Sankey Diagram | SankeyChart | src/charts/SankeyChart.js:241 | | none |
| multi_line_chart | Multi Line Chart (Classic) | MultiLineChart | src/charts/MultiLineChart.js:339 | legacy | none |
| single_heat_map | Interaction Matrix | SingleHeatMap | src/charts/SingleHeatMap.js:371 | | none |
| single_series_chart | Single Interaction Chart | SingleSeriesChart | src/charts/SingleSeriesChart.js:209 | | none |
| cell_radial_chart | Radial Connectivty Map | CellRadialChart | src/charts/CellRadialChart.js:582 | | none |
| cell_network_chart | Spatial Connectivity Map | CellNetworkChart | src/charts/CellNetworkChart.js:954 | | none |
| flexible_network_chart | Network Graph | FlexibleNetworkChart | src/charts/FlexibleNetworkChart.js:879 | | none |
| custom_box_plot | Abundance Box Plot | CustomBoxPlot | src/charts/CustomBoxPlot.js:158 | | none |
| density_box_plot | Density Box Plot | CustomBoxPlot | CustomBoxPlot.js:171 | | none |
| deeptools_heatmap | DeepTools HeatMap | DeepToolsHeatMap | src/charts/DeepToolsHeatMap.js:256 | | none |
| category_heatmap | Category Heatmap | CategoryHeatMap | src/charts/CategoryHeatMap/CategoryHeatMap.js:114 | | none |
| row_summary_box | Row Summary Box | RowSummaryBox | src/charts/RowSummaryBox.js:113 | | none |
| image_table_chart | Image Table | ImageTableChart | src/charts/ImageTableChart.js:175 | | none |
| ImageScatterChart | Image Scatter Plot | ImageScatterChart | src/charts/ImageScatterChart.ts:271 | | none |
| table_chart | Table (Classic) | TableChart | src/charts/TableChart.js:516 | legacy | "Change To Edit Mode" toggle (`fa-edit` / `fa-table`), "Download data"; "Add Column" (`fa-plus`) added when edit mode is entered (TableChart.js:116-134, 259-263) |
| table_chart_react | Table | TableChartReact | src/react/components/TableChartReactWrapper.tsx:317 | | "Download data"; "Add Column" only when project permission is `edit` (TableChartReactWrapper.tsx:116-130) |
| genome_browser | Genome Browser | GenomeBrowser | src/charts/GenomeBrowser.js:613 | | "Show/Hide Tracks" (`fa-eye-slash`) and a location text `<input>` appended straight into `menuSpace`. Both are added in a `.then` after a dynamic import, so after the constructor returns (GenomeBrowser.js:17-42) |
| igv_browser | IGV Browser | IGVBrowser | src/charts/GenomeBrowsers/IGVBrowser/IGVBrowser.tsx:489 | | none |
| ucsc_browser | UCSC Browser | UCSCBrowser | src/charts/GenomeBrowsers/UCSCBrowser/UCSCBrowser.tsx:85 | | none |
| sv_circos_plot | SV Circos Plot | SVCircosPlot | src/charts/GenomeBrowsers/SVCircosPlot/SVCircosPlot.js:35 | | none |
| wgl_scatter_plot | 2D Scatter Plot (Classic) | WGLScatterPlot | src/charts/WGLScatterPlot.js:833 | legacy | none |
| wgl_3d_scatter_plot | 3D Scatter Plot (classic) | WGL3DScatterPlot | src/charts/WGL3DScatterPlot.js:171 | legacy, always converted on load (src/charts/legacyCharts.ts:188) | none |
| density_scatter_plot | Density Scatter Plot | DensityScatterPlot | src/charts/DensityScatterPlot.js:317 | legacy | none |
| image_scatter_plot | Centroid Plot | DensityScatterPlot | DensityScatterPlot.js:337 | | none |
| viv_scatter_plot | Viv Scatter Plot | VivScatterPlot | src/charts/VivScatterPlot.js:315 | allow_user_add false | "Alter Channels" (`fa-palette`) when `config.viv` is set, opens legacy `ColorChannelDialog` (VivScatterPlot.js:204-209) |
| VivMdvRegionReact | Spatial Image Viewer | VivMdvReact | src/react/components/VivMDVReact.tsx:296 | spreads viv_scatter_plot, allow_user_add true | "Alter Channels" (`fa-palette`), opens `ColorChannelDialogReactWrapper` (VivMDVReact.tsx:144-149) |
| SpatialDataMdvRegionReact | SpatialData.js viewer | SpatialDataMdvReact | src/react/components/SpatialDataMDVReact.tsx:521 | `required` needs spatial regions | "Manage Layers" (`fa-layer-group`), opens `SpatialLayerDialogReactWrapper` (SpatialDataMDVReact.tsx:390-398). Menu pip from render status (SpatialDataMDVReact.tsx:382-386) |
| DeckScatter | 2D Scatter Plot (new) | DeckScatterReact | src/react/components/DeckScatterReactWrapper.tsx:400 | allow_user_add false | none in title bar. On-canvas `SelectionOverlay` |
| DeckDensity | Scatter Plot + Density | DeckScatterReact | DeckScatterReactWrapper.tsx:419 | allow_user_add false | same |
| DeckContourScatter | 2D Scatter Plot | DeckScatterReact | DeckScatterReactWrapper.tsx:435 | spreads DeckDensity | same |
| DeckScatter3D | 3D Scatter Plot | DeckScatterReact | DeckScatterReactWrapper.tsx:446 | | same |
| wgl_scatter_plot_dev | 2D Scatter Plot (no density) | DeckScatterReact | DeckScatterReactWrapper.tsx:470 | legacy | same |
| wgl_3d_scatter_plot_dev | 3D Scatter Plot (Deck) | DeckScatterReact | DeckScatterReactWrapper.tsx:486 | legacy, allow_user_add false | same |
| DeckSplatter | Splatter Plot | DeckSplatterReact | src/react/components/SplatterPlotReactWrapper.tsx:526 | preview | none |
| compound_dot_plot | Dot Plot (Compound) | DotPlotReact | src/react/components/dotplot/DotPlotReactWrapper.tsx:444 | | none |
| split_dot_plot | Dot Plot (Compound) | DotPlotReact | DotPlotReactWrapper.tsx:489 | allow_user_add false | none |
| multi_line_chart_react | Multi Line Chart | MultiLineReact | src/react/components/multiline/MultiLineReactWrapper.tsx:417 | | none |
| geneNetwork | Gene / Protein Info | GeneNetworkChartWrapper | src/react/components/GeneNetworkChart.tsx:77 | | none in title bar. Search bar in content |
| highlightedFeature | Highlighted Feature | HighlightedFeatureChartWrapper | src/react/components/HighlightedFeatureChart.tsx:43 | allow_user_add false | none |
| selection_dialog | Selection Dialog | SelectionDialogReact | src/react/components/SelectionDialogReact.tsx:85 | | none |
| text_box_chart | Text Box | TextBoxChartReactWrapper | src/react/components/TextBoxChartReactWrapper.tsx:34 | | none |
| cross_pcf | Cross-PCF matrix | CrossPcfChart | src/spatial_stats/CrossPcfChart.tsx:161 | preview | none |

Registrations that never run or are overwritten:

| Fact | Source |
|---|---|
| TreeDiagram registers under the key `row_summary_box`. RowSummaryBox is imported later and overwrites it, so Tree Diagram is not reachable | src/charts/TreeDiagram.js:90; registerChartModules.ts:16, 24 |
| DensityPlot, MultiPieChart, SingleHeatMapGroup, MultiBoxPlot, VivVolume (viv_volume_scatter_view, viv_volume_view) have registration blocks but no file imports them | src/charts/DensityPlot.js:163; MultiPieChart.js:230; SingleHeatMapGroup.js:307; MultiBoxPlot.js:148; src/webgl/VivVolume.js:195, 219 |
| VivVolume would add "Alter Channels" if loaded | src/webgl/VivVolume.js:11 |

## 4. Callers of ChartManager.addMenuIcon

`ChartManager.addMenuIcon(dataSource, icon, text, func, right)` puts an icon on a data source menu bar (`dsIndex[ds].menuBar`) or, for `"_main"`, on `leftMenuBar` / `rightMenuBar`. It uses `createMenuIcon` and returns the span (src/charts/ChartManager.js:1223-1237; src/utilities/Elements.js:52-84). `leftMenuBar` and `rightMenuBar` are never assigned anywhere in `src`, so a `"_main"` icon is created with no parent and is not shown (Elements.js:80-82). The React `DataSourceMenuBar` mounts in its own `display: contents` span inside `ds.menuBar` so these icons can still be appended after it (ChartManager.js:1485-1497; src/react/components/DataSourceMenuBar.tsx:30-33). The public type promises an `HTMLElement` return (src/charts/charts.d.ts:388).

| Caller | Target | Icon / tooltip | Relies on | Source |
|---|---|---|---|---|
| GridstackManager.getGrid | data source bar | `fa-compress-arrows-alt`, "compact layout" | Keeps the returned element in `grids` and calls `.remove()` on it when the grid is torn down | src/charts/GridstackManager.ts:79-82, 150 |
| hyperion.js | `_main` | `fa-cogs`, "re-cluster selected cells" | Nothing returned is used. Lands nowhere (see above) | src/legacy/hyperion.js:57 |
| hyperion.js | data source `cells` (hard-coded name) | image scatter, upload image, bulk images, abundance chart, Avivator chart | Data source named `cells`; `cm.showCustomDialog` | src/legacy/hyperion.js:385-425 |
| hyperion.js `scAddSSMenuIcon` | each data source with `custom.spatial_stats` | Add Spatial Chart, Adjacency Plot Settings, Single Cell Interaction, Cell Types To Display | `ds.name`, `cm.viewData` | src/legacy/hyperion.js:458-472 |

No file in `src` or `vite.config.mts` imports `src/legacy/hyperion.js`. ADR 0019 lists it as a caller (docs/adr/0019-extract-legacy-js-ui-into-react-in-place.md:15, 46).

Callers of the chart-level `BaseChart.addMenuIcon` from outside the chart class, and what they rely on:

| Caller | Relies on | Source |
|---|---|---|
| ChartManager popout and close icons | Returned span. `popoutIcon` stored on the chart | ChartManager.js:2041-2052 |
| GridstackManager lock icon | Returned span's first child `<i>` to toggle `fa-lock` / `fa-unlock`; stored as `chart.lockButton`; `.remove()` on teardown | GridstackManager.ts:309-329, 396 |
| BaseChart fullscreen | `querySelector("i")` on the span to swap classes and set `aria-label` | BaseChart.ts:302-308, 333-339 |
| TableChart edit toggle | Removes and recreates the child `<i>`, sets `aria-label` on the span | TableChart.js:116-128 |
| BaseChart reset | `style.display` on the span | BaseChart.ts:207-214 |

Charts or helpers that set `resetButton.style.display` directly: CategoryChart.js:70, 122; HistogramChart.js:115, 143; SankeyChart.js:77, 209, 220; DotPlot.js:329, 412; HeatMap.js:159; MultiLineChart.js:131; ViolinPlot.js:43, 187; BoxPlot.js:21; WGLScatterPlot.js:277; WGL3DScatterPlot.js:30; GenomeBrowser.js:143, 280; CategoryHeatMap.js:61, 68; src/react/spatial_context.tsx:259; src/react/components/SelectionDialogComponent.tsx:891; GeneNetworkChart.tsx:32; DotPlotReactWrapper.tsx:76, 92, 99, 115; MultiLineReactWrapper.tsx:133, 141, 149.

## 5. Settings coverage

The settings icon opens `SettingsDialogReactWrapper`, titled "Settings: <title> (<type name>)", which renders `chart.getSettings()` with `SettingsDialogComponent` (src/react/components/SettingsDialogReactWrapper.tsx:27-32). The renderer maps a fixed set of GuiSpec types to components: text, textbox, slider, spinner, dropdown, multidropdown, category_selection, single_category_selection, check, radiobuttons, doubleslider, button, folder, column, multicolumn (src/react/components/SettingsDialogComponent.tsx:937-961; src/charts/charts.d.ts:204-232). There is no GuiSpec type for an arbitrary React component. `BaseChart.getSettings` returns Chart Name, Chart Legend, param settings and a colour folder (BaseChart.ts:930-967). Technical settings from `getDebugSettings()` show in the debug dialog, not the settings dialog (src/react/components/DebugJsonDialogComponent.tsx:397).

| Group | Chart types | Settings UI | Source |
|---|---|---|---|
| GuiSpecs only | bar_chart, row_chart, wordcloud, ring_chart, stacked_row_chart, heat_map, dot_plot, violin_plot, box_plot, sankey_chart, multi_line_chart, single_heat_map, single_series_chart, cell_radial_chart, cell_network_chart, flexible_network_chart, custom_box_plot, density_box_plot, deeptools_heatmap, category_heatmap, row_summary_box, image_table_chart, ImageScatterChart, sv_circos_plot, ucsc_browser, wgl_3d_scatter_plot, density_scatter_plot, image_scatter_plot, compound_dot_plot, split_dot_plot, multi_line_chart_react, DeckSplatter, text_box_chart, highlightedFeature | Settings dialog only | `getSettings` overrides listed in each class, e.g. HistogramChart.js:196, HeatMap.js:333, DotPlotReactWrapper.tsx:357 |
| GuiSpecs with action buttons | flexible_network_chart ("Reset Zoom"), WGLChart family ("Centre Plot"), Deck scatter family ("Center Plot"), igv_browser ("Add Track" with its own name/url text specs), cross_pcf ("Add / remove pair", "Clear pairs") | Still the settings dialog | FlexibleNetworkChart.js:867; WGLChart.js:318; DeckScatterReactWrapper.tsx:385; IGVBrowser.tsx:455-470; CrossPcfChart.tsx:100-111 |
| GuiSpecs plus a context-menu dialog | wgl_scatter_plot and subclasses (`OffsetDialog`) | | WGLScatterPlot.js:615-637 |
| GuiSpecs plus debug settings | Deck scatter family, SpatialDataMdvRegionReact (`getDebugSettings`, `getDebugPanel`) | Debug dialog from the More menu | DeckScatterReactWrapper.tsx:238-241; SpatialDataMDVReact.tsx:309-313 |
| GuiSpecs plus on-canvas tools and gate dialogs | Deck scatter family, VivMdvRegionReact, SpatialDataMdvRegionReact | `SelectionOverlay`, `ManageGateDialogWrapper`, `GateNameDialog` | src/react/components/SelectionOverlay.tsx:219-240 |
| Custom channel dialog | VivMdvRegionReact (`ColorChannelDialogReactWrapper`), viv_scatter_plot (legacy `ColorChannelDialog`) | Opened from the "Alter Channels" title-bar icon | VivMDVReact.tsx:144-149; VivScatterPlot.js:204-209; src/react/components/ColorChannelDialogReactWrapper.tsx:20-45 |
| Custom layer dialog | SpatialDataMdvRegionReact | "Manage Layers" icon and a "Layers" button in settings open `SpatialLayerDialogReactWrapper`. It picks a panel per layer kind: `ImageLayerPanel`, `ShapesLayerPanel`, `PointsLayerPanel`, `LabelsLayerPanel`, `GroupLayerPanel`, `DeckOverlayLayerPanel`. `DeckOverlayLayerPanel` renders GuiSpecs from `getDataLayerSettings` / `getDensityLayerSettings`. The chart's own `getSettings` is "chart-wide settings only" | SpatialDataMDVReact.tsx:390-398, 420-470; src/react/components/SpatialLayerDialogComponent.tsx:205-240; src/react/components/spatialLayers/DeckOverlayLayerPanel.tsx:15, 30 |
| Custom track UI | genome_browser | "Show/Hide Tracks" dialog and location input in the title bar; per-track `getSettings(panel)` inside the browser | GenomeBrowser.js:32-42, 227-238; src/browser/tracks.js:93, 821 |
| Custom column dialogs | table_chart (`EditColumnDialog`, `AddColumnDialog`); table_chart_react (`AddTableColumnDialog`, `FindAndReplaceDialog`, `BulkEditColumnDialog`, `RenameTableColumnDialog`, `ColumnRemovalImpactDialog`) | Opened from title-bar icons or the grid | TableChart.js:90, 261; TableChartReactWrapper.tsx:149; src/react/components/TableChartReactComponent.tsx:4-13 |
| Controls in the chart body | selection_dialog (the chart is the filter controls; settings add one "Filters remain on Reset All" check); geneNetwork (`FeatureSearchBar`) | | SelectionDialogReact.tsx:66-81; src/react/components/GeneNetworkChartComponent.tsx |

Each dialog and each React chart is its own React root through `createMdvPortal`. They share state through MobX-observable `chart.config` and project-scoped singletons, not through a common React tree (docs/design/spatial-tables/03-shared-spatial-contexts.md:15-32).

## 6. Selection tools

`SelectionOverlay` defines five tools: Pan (`TranslateMode`), Rectangle, Polygon, Freehand (lasso icon, `DrawPolygonByDraggingMode`) and Modify (`EditMode`). A Measure tool is commented out (src/react/components/SelectionOverlay.tsx:31-62). It also shows a show/hide-tools toggle saved as `config.selection_toolbar_open`, a density-grid toggle, Manage Gates, Save Gate, and confirm/cancel while editing a gate (SelectionOverlay.tsx:104-140, 260-380). While a gate is edited only Pan and Modify show (SelectionOverlay.tsx:67, 165-167). The chosen tool is React state in `SpatialAnnotationProvider`, not config. Its default comes from the browser drag setting (src/react/spatial_context.tsx:109-115, 171-174).

| Tool set | Charts | Notes | Source |
|---|---|---|---|
| `SelectionOverlay` | DeckScatter, DeckDensity, DeckContourScatter, DeckScatter3D, wgl_scatter_plot_dev, wgl_3d_scatter_plot_dev | In 3D (`config.dimension === "3d"`) gates and density grid are off. The shape filters by screen projection | src/react/components/DeckScatterComponent.tsx:389; SelectionOverlay.tsx:109-111; spatial_context.tsx:205-213 |
| `SelectionOverlay` | VivMdvRegionReact | | src/react/components/VivScatterComponent.tsx:235 |
| `SelectionOverlay` | SpatialDataMdvRegionReact | Not shown in the region grid | src/react/components/SpatialDataMDVReactComponent.tsx:788 |
| Density grid toggle | DeckContourScatter, DeckDensity, SpatialDataMdvRegionReact only | | src/react/components/densityGridUtils.ts:5-9 |
| WGL2DI built-in brush | wgl_scatter_plot, density_scatter_plot, image_scatter_plot, viv_scatter_plot | "Brush Type" setting: Free Draw (`poly`, default) or Rectangle (`default`) | WGLScatterPlot.js:262-283, 815-828 |
| WGL2DI built-in brush | violin_plot, box_plot | Default `poly` | ViolinPlot.js:29-50 |
| WGL2DI built-in brush | wgl_3d_scatter_plot | Default rectangle. Filters on returned indexes | WGL3DScatterPlot.js:19-36 |
| d3 `brushX` | bar_chart | | HistogramChart.js:34-67 |
| React brush per field | multi_line_chart_react | | MultiLineReactWrapper.tsx:99-134 |
| Click on a mark | row_chart, wordcloud, ring_chart, sankey_chart, dot_plot, category_heatmap, compound_dot_plot, split_dot_plot | | CategoryChart.js:57; SankeyChart.js:113, 157; DotPlot.js:506; CategoryHeatMap.js:55-62; DotPlotReactWrapper.tsx:72 |
| Ctrl + drag a range | genome_browser | Without Ctrl the drag zooms | GenomeBrowser.js:134-146 |
| Search terms | geneNetwork | | src/react/components/featureSearchHooks.ts:25-41 |
| Form controls | selection_dialog | | SelectionDialogComponent.tsx |
| Colour legend click or brush | any chart with a colour legend, unless `isColorLegendFilterable()` is false (compound/split dot plot) | | BaseChart.ts:620, 721; DotPlotReactWrapper.tsx:178 |
| Points feature list | SpatialDataMdvRegionReact points layers | Chooses which features a points layer draws. Layer config, not a DataStore filter | src/react/components/spatialLayers/PointsFeatureFilterPanel.tsx; src/react/spatialdata/points_layer_config.ts:31-34 |

## 7. Filter mechanisms

General plumbing:

| Fact | Source |
|---|---|
| BaseChart listens for `filtered` on its DataStore and calls `onDataFiltered(dim)`. On `"all_removed"` it also clears the colour-legend filter | BaseChart.ts:217-224 |
| "Reset All Filters" in the data source bar calls `DataStore.removeAllFilters()`, which skips dimensions marked `noClear` and fires `"all_removed"`. It does not call `chart.removeFilter` | src/react/components/DataSourceMenuBar.tsx:121; src/datastore/DataStore.js:274-290; spatial_context.tsx:242-251 |
| `BaseChart.getFilter()` returns undefined by default. It is used by `updateResetButtonVisibility` and by `ChartManager.getAllFilters`. Only hyperion.js calls `getAllFilters` | BaseChart.ts:448, 663-669; ChartManager.js:2209-2219; hyperion.js:83, 184 |
| Charts that apply a saved filter after construction declare it with `storedFilters.expect` | BaseChart.ts:386; src/charts/viewFilterReadiness.ts:65-72; GeneNetworkChart.tsx:30; SelectionDialogReact.tsx:55 |

Grouped by mechanism:

| Mechanism | Chart types | Where the filter is stored | Set | Cleared | Saved in config? | Enough state to say it in words? |
|---|---|---|---|---|---|---|
| Category list on one column (`category_dimension`, `filterCategories`) | row_chart, wordcloud, ring_chart | `this.filter` (array of category values) | `filterCategories` (CategoryChart.js:57-73) | `removeFilter` (CategoryChart.js:51-55); `"all_removed"` (CategoryChart.js:121-124) | No | Yes. `getFilter()` returns `{column: [values]}` (CategoryChart.js:75-82), e.g. "Tissue: Gut, Skin" |
| Sankey node or link (`category_dimension`) | sankey_chart | `filter`, `nodeFilter` (`"A3"` style), `linkFilter` | `filterCategories` / `filterMultipleCategories` (SankeyChart.js:202-223) | `removeFilter` (SankeyChart.js:58-64) | No | Yes. `getFilter()` decodes node or link to `{column: [value]}` (SankeyChart.js:180-200) |
| Category pair (`filterCategoryPair`) | category_heatmap | `this.cellFilter = [x, y]` | CategoryHeatMap.js:55-62 | `removeFilter` (CategoryHeatMap.js:65-69) | No | Yes. `getFilter()` returns both columns (CategoryHeatMap.js:71-83) |
| Category + field threshold (`catcol_dimension`, `filterCatCol`) | dot_plot | Only in the dimension (`filterArguments`, `filterColumns`). `this.filter` is reset but never filled | DotPlot.js:406-413 | `removeFilter` (DotPlot.js:210-214) | No | Partly. No `getFilter`. Category, field and threshold could be read from the dimension (DotPlot.js:524) |
| Predicate for one dot (`base_dimension`, `filterPredicate`) | compound_dot_plot, split_dot_plot | `this.selectedDot` (MobX observable on the chart) | `filterToDot` (DotPlotReactWrapper.tsx:72-93) | `removeFilter`, `clearDotFilter`, `"all_removed"` (DotPlotReactWrapper.tsx:95-117) | No | Partly. No `getFilter`. `selectedDot` holds row, field and panel; threshold is in config |
| 1D range (`range_dimension`, `filterRange`) | bar_chart | `this.filter = [min, max]` | `setRange` (HistogramChart.js:107-116) | `removeFilter` (HistogramChart.js:82-88) | No | Yes. `getFilter()` returns `{column: [min, max]}` (HistogramChart.js:118-125) |
| Per-field ranges (`range_dimension` per field) | multi_line_chart_react | `config.brush` (field to range), `appliedRanges` | `applyBrushFilters` (MultiLineReactWrapper.tsx:99-134) | `removeFilter` clears all brushes (MultiLineReactWrapper.tsx:137-142) | Yes | Yes, from `config.brush`. No `getFilter` |
| WGL rectangle or polygon on x/y (`range_dimension`, `filterSquare` / `filterPoly`) | wgl_scatter_plot, density_scatter_plot, image_scatter_plot, viv_scatter_plot | `this.range` (rectangle) or `true` (polygon); `this.filter` (ranges or points) | brush handler, `_createFilter`, `_createPolyFilter` (WGLScatterPlot.js:276-283, 533-574) | `WGLChart.removeFilter` clears brush and dimension (WGLChart.js:169-177) | No | Rectangle yes: `getFilter()` gives x and y ranges (WGLScatterPlot.js:503-511). Polygon: `getFilter()` returns null |
| WGL brush on category x value (`catrange_dimension`) | violin_plot, box_plot | `this.range`, `this.filter` | ViolinPlot.js:42-50, 154-170 | `WGLChart.removeFilter` | No | No `getFilter`. The rectangle branch calls `this.dim(...)` rather than `this.dim.filter(...)` (ViolinPlot.js:167) |
| 3D index filter (`filterOnIndex`) | wgl_3d_scatter_plot | Only in the dimension | WGL3DScatterPlot.js:29-36, 142-144 | `WGLChart.removeFilter` | No | No. A list of row indexes. The poly branch calls `_createPolyFilter(range)` with `range` undefined (WGL3DScatterPlot.js:35) |
| Drawn shape on x/y (`range_dimension`, `filterPoly`, or `filterScreenPoly` in 3D) | Deck scatter family, VivMdvRegionReact, SpatialDataMdvRegionReact | `config.selectionFeatureCollection` (GeoJSON, at most one shape). Gate being edited is React state `editingGateId` | effect in `SpatialAnnotationProvider`, region-restricted when a gate has a region (spatial_context.tsx:162-166, 253-296) | `chart.removeFilter` is replaced with `clearSelection` (spatial_context.tsx:233-240); `"all_removed"` listener (spatial_context.tsx:242-251); empty-click on the layer | Yes | Partly. A polygon in data coordinates on two param columns. A gate has a name in `GateManager` (src/react/hooks/useGateActions.ts:38-56) but the saved selection does not record it |
| Linked index filter (`base_dimension` on a linked store) | genome_browser | `cellDim` and a browser highlighted region `_filter` | GenomeBrowser.js:137-145 | `removeFilter` (GenomeBrowser.js:506-510) | No | No `getFilter`. The range is only in the browser highlight |
| Text / term search (`filterPredicate`) | geneNetwork | `config.search` (`query`, `terms`) | featureSearchHooks.ts:25-41 | `removeFilter` calls `clearFeatureSearch` (GeneNetworkChart.tsx:35-37) | Yes | Yes, query and terms |
| Per-column form filters | selection_dialog | `config.filters` keyed by column: `{category: [...]}`, `{category, operand}`, a unique string, or `[min, max]` | column components in SelectionDialogComponent | `removeFilter` nulls every entry (SelectionDialogReact.tsx:58-64). Reset button follows `hasFilter` (SelectionDialogComponent.tsx:888-892). `noClearFilters` keeps them on Reset All | Yes | Yes, typed per column (SelectionDialogReact.tsx:8-18) |
| Colour legend filter (`category_dimension` or `range_dimension`) | Any chart with a filterable colour legend | `chart.colorLegendFilter` and `config.color_legend.filter`: `{kind: "categorical", column, value}` or `{kind: "continuous", column, range}` | `toggleCategoricalColorLegendFilter`, `setContinuousColorLegendFilter`, restored on load by `restoreColorLegendFilter` | Reset icon calls `clearColorLegendFilter` (BaseChart.ts:209-213); `"all_removed"` | Yes | Yes. One column and one value or range (src/react/legend/color_legend/colorLegendFilter.ts:7-19, 84-97, 173-183) |
| Chart-scope display filters | Deck scatter family, VivMdvRegionReact, SpatialDataMdvRegionReact (and legacy scatter `background_filter`) | `config.background_filter`, `config.category_filters` (column plus categories) | Settings dropdowns (VivMDVReact.tsx:172-188; SpatialDataMDVReact.tsx:461-478) | Settings only. The reset icon does not touch them | Yes | Yes. They narrow what this chart draws. They are not a DataStore filter (src/react/hooks.ts:34-70) |
| No filter (dimension used for aggregation only) | heat_map, stacked_row_chart, multi_line_chart, custom_box_plot, density_box_plot, deeptools_heatmap | | | Some hide the reset icon on `"all_removed"` (HeatMap.js:157-160; MultiLineChart.js:130-132) | | Nothing to describe |
| No filter (highlight only or none) | table_chart, table_chart_react, image_table_chart, cell_network_chart, flexible_network_chart, single_heat_map, single_series_chart, cell_radial_chart, row_summary_box, ImageScatterChart, igv_browser, ucsc_browser, sv_circos_plot, DeckSplatter, highlightedFeature, text_box_chart, cross_pcf | | `dataHighlighted` in some (TableChart.js:334; FlexibleNetworkChart.js:320, 438) | | | Nothing to describe |

Pinning: `pinChart` makes the chart ignore other filters until unpinned (CategoryChart.js:84-91; WGLChart.js:216-228). It does not change the chart's own filter.

## Not verified

- Runtime order of icons was read from code, not checked in a browser.
- Whether `ViolinPlot.js:167` and `WGL3DScatterPlot.js:35` throw in practice was not run.
- Each GuiSpecs-only chart's `getSettings` was located but not read line by line for hidden custom behaviour.
