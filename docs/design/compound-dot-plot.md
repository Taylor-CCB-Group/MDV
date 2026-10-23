# Compound Dot Plot — prototype plan

A new React/TypeScript chart that draws the same dot plot once per value of a
**split by** column (e.g. Control | Disease), side by side, on **one shared
colour and size scale**. Without a split it is an ordinary dot plot, so it is
also the first step of moving the legacy `DotPlot.js` to React.

## Why

Biologists compare expression across conditions by putting dot plots side by
side. Two legacy dot plots with background filters can do the layout today, but
each one scales colour (`mean_range`) and dot size (`frac_range`) to its own
data ([DotPlot.js:248](../../src/charts/DotPlot.js), [DotPlot.js:370](../../src/charts/DotPlot.js)),
so the same shade means different values in each panel and the comparison is
misleading.

## Scope of the first prototype

In:

- Rows = category column, columns = genes/fields (same parameters as `dot_plot`)
- **Split by** column + ordered **split values** → one panel per value
- One colour scale and one size scale across all panels
- The existing colour and fraction legends (draggable, resizable), starting in a band above the panels
- Scale: auto (shared) or manual min/max, palette, log
- Y-axis order: data or alphabetical
- Basic axis controls
- Hover tooltip
- Reacts to global filters

Out (later, if asked for; also what the legacy `dot_plot` still has that this lacks):

- Click a dot to filter cells (legacy `filterCategories`)
- Row/column clustering (must run across all panels together)
- Background filter
- Interleaved or stacked layouts
- Independent scales per panel
- Auto-fit / shorten / hide overlapping X labels (reuse `src/charts/dotplot/xAxisLabels.ts`)
- Threshold UI (`threshold_value` exists, no setting today)
- Replacing / migrating legacy `dot_plot` configs

## One chart, not two

This is the React/TS dot plot, with **split by** optional. It is not a separate
chart next to a future React dot plot. It is built in three layers:

1. **Compute** (`dotPlotStats.ts`): pure TS, groups by row × optional split.
2. **`DotPlotPanel`**: a plain React component that draws one dot plot from
   data, scales and size. No settings or data loading.
3. **Chart**: one `BaseReactChart` class that loads data, holds settings, and
   draws the legend strip and N panels (N = 1 with no split).

Two chart classes would duplicate settings, legends, axes and data loading.
Putting whole charts inside a chart doesn't fit MDV either: each chart has its
own title bar, settings and filters, so panels should be components.

Replacing the legacy `dot_plot` later follows the `DeckScatterReact`
precedent: one class registered under several type names, with an
`adaptConfig()` that reads legacy configs. That happens once the features
listed under "Later" exist. It isn't a prerequisite for the workshop.

## Settings

| Folder | Setting | Notes |
|---|---|---|
| General | Chart name, Chart legend | from `BaseChart.getSettings` |
| General | Y-axis order | data / alphabetical |
| General | Show colour legend, Show fraction legend | |
| Parameters | Categories on y-axis, Fields on x axis | from `getParamsGuiSpec`, unchanged |
| Parameters | **Split by** | `dropdown` of text/text16 columns + "None" |
| Parameters | **Split values** | existing `category_selection` control (the chip autocomplete used by "Background Filter Categories"); chip order = panel order |
| Scale | Colour range | Auto / Manual, with Manual min and Manual max |
| Scale | Palette, Log colour scale | |
| Axis controls | X/Y text size, Y axis width, X axis height, rotate X labels | `getAxisGuiSpec()`, see below |

Split by and Split values are appended to the existing "Parameters" folder
returned by `super.getSettings()`.

Dropped on purpose: background filter, cluster rows/columns, right Y axis size
(only used by clustering), trim to percentile. The legacy dot plot's
`color_scale.trim` is saved but never read, so the trim setting does nothing
today.

When a split column is first chosen, preselect its first two values.

## Config

```ts
type DotPlotReactConfig = BaseConfig & {
    type: "compound_dot_plot";     // "split_dot_plot" still accepted; later also the legacy-replacement type
    param: FieldSpecs;            // [rowColumn, ...fields], same shape as dot_plot
    split_by: FieldName | null;   // text / text16 column
    split_values: string[];       // ordered, one panel each
    y_axis_order: "data" | "alphabetical";
    threshold_value?: number;     // "expressed if > x"; no UI, default 0 (as legacy)
    color_scale: {
        palette: "reds" | "blue_white_red";
        log: boolean;
        mode: "auto" | "manual";
        min?: number;
        max?: number;
    };
    color_legend: { display: boolean; pos?: [number, number] };
    fraction_legend: { display: boolean; position?: [number, number] };
    axis: AxisConfig2D;           // from scatter_state, shared with getAxisGuiSpec
};
```

A matching zod schema is registered next to `DotPlotConfigSchema` in
[ChartConfigSchema.ts](../../src/charts/schemas/ChartConfigSchema.ts).

## How it computes

The definitions stay the same as the legacy chart's `addSimpleMean`
([catColWorker.js:241](../../src/datastore/catColWorker.js)), so numbers match
the legacy dot plot:

- **mean** = sum of values `> threshold` / number of cells in the group (zeros included)
- **fraction** = cells with value `> threshold` / cells in the group × 100

The only change is the grouping key: *(row category × split value)* instead of
*row category*. It's one pass over the rows:

```
for each row i that passes the global filter (ds.filterBuffer[i] === 0):
    p = panelIndexOf[splitCodes[i]]  (skip if not a selected split value)
    r = rowCodes[i]
    count[p][r]++
    for each field f: if value > threshold → sum[p][r][f] += value, n[p][r][f]++
```

It also returns the shared `meanRange` (min/max over every panel, row and
field) and `fracMax`. Rows with zero cells in a panel are kept, so every panel
has the same rows.

No `Dimension` is needed: the chart doesn't filter in v1, so the worker only
needs the global filter buffer and the column buffers.

## Files

Reference charts:
- [DeckScatterReactWrapper.tsx](../../src/react/components/DeckScatterReactWrapper.tsx):
  wrapper class, `adaptConfig()` defaults, `getSettings()` adding the axis folder,
  one class registered under several types
- [SplatterPlotReactWrapper.tsx](../../src/react/components/SplatterPlotReactWrapper.tsx):
  a grid of small multiples with row and column labels, `useFieldSpec()` for an
  extra category column, `extra_controls` + `init` in the add-chart dialog

The worker follows the raster-histogram pattern from #565 (pure module + thin
worker + query function).

```
src/react/components/dotplot/
    DotPlotReactWrapper.tsx     class extends BaseReactChart: adaptConfig, getSettings,
                                BaseChart.types registration
    DotPlotComponent.tsx        observer: data fetch, layout, legend strip, N panels
    DotPlotPanel.tsx            one panel: dots, X axis, optional Y axis, tooltip
    dotPlotStats.ts             pure compute (above), no DOM, unit tested
    dotPlotStatsWorker.ts       self.onmessage → computeDotPlotStats
    dotPlotConfig.ts            config type, defaults, colour range choice
    dotPlotLayout.ts            visible rows, panel layout, dot radius
    dotPlotLegends.ts           legend band and fraction legend spec
    dotPlotPalettes.ts          palettes (moved out of DotPlot.js), colour binning
src/tests/charts/dotPlot{Stats,StatsLegacyParity,Layout,Config,Legends}.test.ts
```

Plus one import in [registerChartModules.ts](../../src/charts/registerChartModules.ts)
and the schema registration.

### Axis controls

Reuse `getAxisGuiSpec()` from
[AxisSettingsGui.ts](../../src/charts/dialogs/utils/AxisSettingsGui.ts), as
`DeckScatterReact` does. Two small changes:

- Loosen its parameter from `ScatterPlotConfig2D` to `{ axis: AxisConfig2D }`.
  It only reads `config.axis.x/y`.
- Turn on the commented-out "Rotate labels" control behind an option, so
  scatter plots don't get it. Dot plots need it.

Axes are drawn with `@visx/axis` on band scales, as
[AxisComponent.tsx](../../src/react/components/AxisComponent.tsx) does for scatter plots.

## Rendering

- Defaults come from an `adaptConfig()` (as `DeckScatterReact` does): axis
  sizes, legend display, colour scale, `y_axis_order`, `split_values: []`.
- Columns are loaded with `useParamColumnsExperimental()`, which flattens
  linked/active-link fields. The split column is loaded with `useFieldSpec(config.split_by)`.
- Recompute on filter changes and on relevant config changes. Each request
  creates a fresh worker; stale results are dropped with a cancelled flag
  (same as `CategoryHeatmapComponent`; or `useFilteredIndices()` as the Splatter
  plot does, if a filter-version dependency is enough).
- Plain SVG in React. Band and sqrt scales come from `@visx/scale`, which is
  already a dependency. Size comes from `useChartSize()`.
- Layout: a legend strip on top, then panels in a row. Panel width =
  (width − Y axis width − gaps) / panel count. Y labels are drawn only on the
  first panel; each panel has its split value as its title.
- Legends: the existing `ColorLegend` (through `BaseChart.setColorLegend()`) and
  `FractionLegend` (in a second `LegendWrapper`, with the new `orientation: "horizontal"`).
  The component keeps a band free above the panels and the legends start there; they
  can be dragged and resized, and their positions are saved. On a narrow chart the
  fraction legend starts below the colour legend.
- The colour legend's brush is off (`isColorLegendFilterable()`): it would filter
  cells by the first field's values, which isn't what mean expression shows.
- Auto colour range is fixed while filtering and recalculated when the rows,
  fields, split or threshold change, as in the legacy dot plot.
- Palettes: move the two palettes out of `DotPlot.js` into a shared TS module
  so both charts use the same colours.

## Status

Done: compute + legacy parity test, React dot plot with optional split, shared
colour and size scale, existing legends in a band above the panels, Scale folder
(auto/manual range, palette, log), axis controls, tooltip, Split by and Split
values in the add-chart dialog (values start empty; hidden while Split by is None).
Tried on the synthetic Control/Disease project and on the myeloid project
("Contrasting R vs NR … (split dot plot)" view).

## TODO

- **Scale each gene (0–1).** A colour option alongside mean expression: each
  gene's min/max taken across all panels (scanpy's `standard_scale="var"`). On
  the myeloid R/NR view one gene (S100A8 in Neutrophils, mean ~90) whitens out
  every other gene (means 0–8) on a shared scale.
- **Fix dot size to 0–100%.** Checkbox, so a dot size means the same % in every
  chart, not 0 → largest % in this chart. Consider making it the default.
- **Colour range as a double slider** instead of the Manual min / Manual max
  boxes (`doubleslider` exists). Keep a way to type values or use a log slider:
  with a 0–90 range a linear slider can't pick 0–8 accurately.
- **Limit the number of x-axis fields?** To decide. A hard cap of 5 would stop
  ordinary dot plots with 10–30 marker genes, which is common. Alternatives: a
  limit that shrinks as panels are added, or no limit and a message when the
  columns get too narrow to read. Linked gene queries already have `maxItems`.
- **Keep legends inside the chart when it shrinks.** `LegendWrapper` places a
  legend at a fixed left/top; `makeDraggable` keeps it inside the chart while
  dragging, but nothing moves it back in on resize, so a legend near the right or
  bottom edge disappears when the chart gets smaller. It affects every chart using
  `LegendWrapper` (the legacy dot plot too), so the fix probably belongs there:
  pull each legend back inside the content area when the chart is resized.
- Dropped: capping split values at 4 panels.

Separate task, found along the way (own branch and PR off `main`): the column
dropdown (`ColumnDropdownComponent.tsx`) keys MUI options by display name, so
two columns with the same name (`sub_bucket` and `sub_bucket_caf`, both "Cell
Subgroup") show up as extra, identical rows. Fix: `getOptionKey={(c) => c.field}`,
and show the field id when names clash.

### Parity with the legacy dot plot

The new chart should do everything the legacy `dot_plot` does, then more.

Priority (agreed):

1. **Most important:** X axis label options (done), click a dot to filter cells,
   export as SVG / PNG.
2. **Nice to have:** Rotate Y axis labels, title follows the y-axis column.
3. **Undecided:** cluster rows / columns.
4. **Not needed:** background filter. Split by covers what it was used for (the
   R/NR view used four background-filtered dot plots), and does it on one scale.

What the legacy chart has that this one doesn't yet (from `DotPlot.js`, plus what
it inherits from `SVGChart` and `BaseChart`):

- **X axis label options (done):** Auto-fit X axis labels (on by default),
  Shorten long X axis labels (on), Hide overlapping X axis labels (off), at the
  top of Axis controls in both dot plots, with a full-name tooltip on the X
  labels. Shorten and Hide overlapping now also work with auto-fit off
  (`layoutFixedXLabels`), in both dot plots.
- **Click a dot to filter cells:** cells in that row category (and split value)
  expressing that field above the threshold, with the selected dot outlined and
  the reset button shown. The legacy chart uses `filterCatCol` on a
  `catcol_dimension`.
- **Cluster rows / cluster columns,** with dendrograms. Must cluster across all
  panels together so every panel keeps the same order.
- **Export as SVG / PNG image** from the chart menu (`getImage`, which
  `SVGChart` provides), including the legends.
- **Rotate Y axis labels** (axis controls).
- **Title follows the y-axis column** when it hasn't been renamed (legacy
  `setParams`).

Deliberately not carried over: background filter (split by replaces it), Trim to
percentile (saved but never read by the legacy chart), the colour legend brush
(filters by the first field's per-cell values), row animations.

More than the legacy chart, beyond the split itself: a way to fit long Y labels
(e.g. `annotation_caf`'s cell states) such as an auto-fit Y axis width, and a
threshold setting (`threshold_value` exists in both, with no UI).

Before merging: Playwright smoke test (read `docs/playwright/PLAYWRIGHT_AGENT.md`
first), synthetic demo project generator.
