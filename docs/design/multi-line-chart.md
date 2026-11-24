# Multi-line chart in React, with `values` and a grid — plan

Status: plan (2026-11-20). Spread over several PRs; each one stands on its own. Tracked in #694.

## Why

A marker QC view we have in use has 13 `multi_line_chart`s, one per marker, each with
`param: ["zscore_clipped_<marker>", "sample_id"]` and identical settings
(`scaletrim: "0.05"`, `band_width: 0.1`, `intervals: 100`). That is one chart with 13 values
drawn as small multiples, which is what `grid: { by }` (#683) was made for. The same view also has
a `heat_map` over the same 13 columns, which is the pairing the cross-PCF prototype (#645)
suggests for composite charts.

Today's chart ([MultiLineChart.js](../../src/charts/MultiLineChart.js)) is a D3 `SVGChart` that
takes one value column and one category column through `config.param`.

## Target config

Semantic props, not `param` (as `cross_pcf` in #645 does):

```ts
type MultiLineConfig = BaseConfig & {
    type: "multi_line_chart";
    /** Category column: one line per value of it. */
    categories: FieldName;
    /** Numeric fields: one density curve per category, per field. A live query or column group works too. */
    values: FieldSpecs;
    /** `values`: a cell per field (default when there is more than one). `categories`: a cell per category, a line per field. */
    grid?: { by: "values" | "categories" } | null;
    stacked: boolean;
    fill: boolean;
    /** Trim the x axis to this quantile; null for the column's full range. */
    scaletrim: "0.001" | "0.01" | "0.05" | null;
    band_width: number;
    intervals: number;
    color_legend: { display: boolean; pos?: [number, number] };
};
```

Notes:

- `ChartGridConfigSchema.by` is a single enum (`"densityFields"`) shared by every chart. It
  should become per chart type: each chart says which of its props can make cells.
- The schema says `scaletrim` is a boolean; configs store a quantile string. Fix it with the schema.
- With several values and `grid: null`, overlaying values × categories is unreadable. Proposal:
  treat it as `{ by: "values" }`, so `grid` only chooses the layout.
- y is already normalised per line (each curve is scaled to its own peak), so cells can share a
  y axis for free. x is not: with `scaletrim` the domain is the trimmed range of each column.
  Default to an x scale per cell, with a shared option for comparable units (z-scores, as here).

## PRs

### 1. React chart, one value (parity) — #690

Same three layers as the compound dot plot ([compound-dot-plot.md](compound-dot-plot.md)):

1. **Compute** (`multiLineStats.ts` + worker): kernel density per category for one value column,
   pure TS, testable against the legacy `boxPlotWorker` output. Each request gets its own reply.
   `CatRangeDimension.getKernalDensity` reassigns `worker.onmessage` per call, so concurrent
   requests (one per cell) would drop all but the last. It also overwrites `config.scaletrim`
   with the quantile value it looks up.
2. **`MultiLinePanel`**: plain React SVG, draws curves from data, scales and size. No loading or
   settings, so a composite chart can embed it (see 6).
3. **Chart**: `BaseReactChart` with `categories` / `values`, settings at parity (categories,
   stack, fill, trim, band width, intervals, legend), tooltip, export.

Registered as `multi_line_chart_react`, named "Multi Line Chart". The legacy `multi_line_chart` is
now "Multi Line Chart (Classic)", marked legacy, and converts to it (`legacyCharts.ts`), moving
`param` to `categories` and `values`.

### 2. Several values, and `grid: { by: "values" | "categories" }` — #691

- Cells through `ChartArrayLayout` / `useChartArrayGrid`, which the density grid already uses
  and which already reports `visibleCellIndices`.
- One category legend for the chart, not one per cell.
- Acceptance: the marker QC view's 13 charts as one chart, matching the old ones by eye.

### 3. Take over `multi_line_chart` — #692

The legacy conversion landed with 1. It runs with the "Convert legacy charts" setting or on request,
like the classic table's. Left: whether it should always run (`always: true`, as for the classic 3D
scatter), removing `MultiLineChart.js`, and optionally a view upgrade that merges sibling charts with the same `categories` and
settings into one grid chart; not automatic.

### 4. Lazy `values` — #693

The chart loads a value's data when its cell is about to be seen, not all of them before it
draws. Field names, colours and min/max/quantiles come from column metadata without loading
data, so labels, layout and x scales are ready at once. Cells without data show a placeholder.

- A hook over `useFieldSpecs`: given the specs and the indices wanted now (visible plus some
  overscan), returns `(LoadedDataColumn | undefined)[]` in spec order, and loads the rest as they
  come into view. Loaded columns stay loaded.
- `chartManager.loadColumnSet` takes one callback for a whole set. Scrolling a 1000-gene grid
  needs requests that can be dropped when cells scroll away, and the most recent ones first.
- Then the same hook for the density grid (`useDensityGridCellLayers` already filters by
  `visibleCellIndices`) and for selection dialog entries.
- Anything that needs every value at once (a shared auto scale computed from data, clustering)
  has to come from metadata or wait for all of them; flag it per feature.

### 5. Column groups as a `values` source (#413)

`DataStore.columnGroups` exists (and `ChooseColumnDialog` reads it) but nothing else uses it, and
conversion rarely writes any (the marker QC project above has none).

- A `FieldSpec` variant naming a group, resolved live like `RowsAsColsQuery`, so `values` can be
  "the markers" and pick up columns added to the group.
- Pick a group in the column selector.
- Conversion writes groups: for SpatialData, a table's var axis (e.g. proteins) becomes a group.

### 6. Composite charts (no PR yet)

The constraint on 1–2: panels (`MultiLinePanel`, a future heatmap panel) take data and scales as
props and know nothing about charts, so a chart like `cross_pcf` can hold several. How values
flow between sub-charts belongs to the dataflow-graph sketch in #645
(`docs/design/dataflow-graph.md` on that branch).

## Semantic config: a checklist for chart refactors

Worth recording as an ADR, as this is the SOP and not written down anywhere yet:

- Props named for their role (`categories`, `values`, `split_by`), typed in the chart's zod schema
- `adaptConfig` maps legacy `param` to them, so old views keep opening
- Loading through `useFieldSpecs` (or the lazy hook), not `@loadColumnData` on `param`
- Settings use FieldSpec controls bound to the props, not `getParamsGuiSpec`
- Add-chart dialog builds the props
- Links that write a chart's `param` (`chart_columnval_link` with `set_param`) call the chart's
  `setParamAt(index, column)` when it has one, so implement it to set the matching prop
