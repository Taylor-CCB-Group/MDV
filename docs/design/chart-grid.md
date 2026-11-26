# Grids of small multiples across charts: draft

Status: draft for discussion (2026-11-24). No decision yet. This note collects what four charts
already do, where a generic grid would help, and what could go wrong, before anything is
generalised.

## Why

Several charts already draw small multiples, each in its own way:

| Chart | What makes the cells | Layout | Per-cell state | Shared by all cells |
| --- | --- | --- | --- | --- |
| Deck scatter / SpatialData density grid (#683) | `grid: { by: "densityFields" }`: a cell per numeric field | `ChartArrayLayout`, one deck per cell (`DeckCellGrid`) | the cell's density layer | points, image and shape layers, view state, selection tools |
| Multi-line chart (#691) | `values`: a cell per numeric field once there are several | `ChartArrayLayout`, an SVG per cell | densities, x domain (unless `share_x`), a brush and its range dimension | categories, colours, legend, settings |
| Compound dot plot | `split_by` / `split_values`: a panel per category value | its own panel layout | the rows of that category value | columns, colour scale (`auto` is shared) |
| Splatter plot | categories × `densityFields`: a cell per pair | its own layout | the category's rows, the field's density | view state |
| Deck scatter / SpatialData colour grid | `grid: { by: "color_fields", split }`: a cell per colour column (`columns`), or per category of the first (`categories`) | `ChartArrayLayout`, one deck (`DeckCellGrid`) | the cell's points: all of them coloured by its column, or its category's | grey points, gates, selection, view state; in SpatialData, the render stack but for its points |
| SpatialData layer colour grid | `grid: { by: "layer_fill", layer, split }`: a cell per fill colour column of a shapes or labels layer, or per category of the first | `ChartArrayLayout`, one deck (`DeckCellGrid`) | its copy of the layer: every feature coloured by its column, or its category's | the render stack below and above the layer, view state |
| SpatialData region grid | `grid: { by: "region" }`: a cell per region (a store and coordinate system) | `ChartArrayLayout`, one deck (`DeckCellGrid`) | the region's store, renderer and rows; its view state | the render stack (as a template), colours, filters |

They answer the same questions (what is a cell, what does it share, where does its state live)
differently, and only the first two use `ChartArrayLayout`. A generic grid could make a new
chart's small multiples a matter of naming a prop, and make the layout, virtualisation and
settings the same everywhere. Facets by category (`grid: { by: "categories" }` for the
multi-line chart) are the next case.

## Kinds of cell

A grid can make cells in at least three ways, and they behave differently:

1. **Split by field.** A list of fields (`values`, `densityFields`). Every cell sees the same
   rows; each draws a different column. Cell keys are field names.
2. **Facet by category.** A category column (`split_by`, `categories`). Every cell draws the same
   columns over a different subset of rows. Cell keys are category values, which come from the
   data, can number in the hundreds, and change when the column is edited.

3. **Split by layer.** In spatial charts, cells can differ in what they draw from the
   SpatialData object rather than from MDV's tables: one cell per **image channel** (a marker
   each), or per **segmentation** (nuclei, cell boundaries, another algorithm's shapes) over the
   same image. Every cell shares the view, coordinate system and usually the other layers; each
   draws a different render stack entry, or the same image with a different channel selection.
   Cell keys are element names or channel names, not columns.

`color_fields` (the deck scatter plot and the SpatialData chart) is the first property made
multi-column so that it could be a grid: it replaces the single `color_by` there, and its grid
can split it either way, by column (kind 1) or by the categories of its first column (kind 2).
A SpatialData layer's fill colour is the second, and the first inside the render stack.

Such a property's picker says so (`GuiSpec.gridBy`, and `GridByControl` beside a layer's colour
picker): it takes one column, with a "Grid by this" switch, and turning that on makes the chart a
grid of it. The split follows the first column's type, not a choice: a categorical column (one
category per row) gives a cell per category, and the picker stays single; a numeric one gives a
cell per column, and the picker takes several numeric columns. Choosing a column of the other
kind while it is on changes the split, and one of neither kind (`multitext`, `unique`) turns it
off. Turning it off leaves the columns in config, the first in use. The Grid dropdown offers only
the split the column allows (`colorGridSplitFor`), and the chart's current one.

With a cell per category, each cell can also draw the other categories' points in grey under its
own (`grey_other_categories`, a chart setting rather than part of `grid`, so it outlasts turning
the grid off and on). It is one layer of all the points, drawn in every cell, not a layer per cell.
Layer fill grids do not have it yet.

The selection toolbar's grid button follows whichever grid is shown: with one, it goes back to one
plot; without, it shows the chart's only grid, or a menu of them (`gridChoices` on the chart, the
same list as the Grid dropdown, less the region grid, and density fields only when there are
some).

A **sweep** of a setting (several band widths, several thresholds), is conceivable but has no
use case yet; leave it out.

The list of cells need not be fixed in config. Channels might be driven by another mechanism,
such as a protein picker, as columns might come from a live query or a column group (#413). So a
cell source can be live, and the grid follows it, adding and dropping cells (and their per-cell
state) as it changes.

Two props can make cells at once (the splatter plot's categories × fields). That is a 2D grid,
not a longer list, and wants rows and columns labelled as such.

## What a contract could look like

Not arbitrary props. Each chart type keeps an explicit list of the props its grid can use, as
`chartGridConfigSchema([...])` now does in the config schema, and for each says:

- **Kind**: field split, category facet or layer split, which decides how cells are enumerated
  and keyed, and whether the source is static config or live.
- **How a cell narrows its inputs**: the field, the row subset, or the render stack entry or
  channel selection it draws. A cell gets inputs, not a
  copied config (see the risks below).
- **What is shared and what is per cell**: axes and scales, view state, legend, colour scale,
  filters and selections. Each chart answers this explicitly rather than inheriting a default.

The layout side is already generic: `ChartArrayLayout` / `useChartArrayGrid` with
`visibleCellIndices`. What isn't is everything above, and that's where charts differ.

## Risks

These are the unexpected consequences to check before generalising.

- **Config edited in place.** Charts and the settings dialog write config through MobX, often to
  nested objects (`axis.x.size`). A cell built from a derived copy of the config would not see
  those edits, or would see them a render late. Cells should read the chart's one config and
  receive only their narrowed inputs as props.
- **Settings.** The dialog edits "the" chart. Per-cell settings (a different trim per cell) would
  need a way to address a cell. Avoid them until something needs them.
- **Links that write props.** `chart_columnval_link` and `set_param` links write a chart's prop
  (the multi-line chart maps `param[0]` to `values`, replacing all of them with one). With a
  grid, the link has to say whether it replaces the cells or one of them.
- **Filtering and selection.** A chart leaves out its own filter as a whole
  (`useOwnedFilteredIndices`). The multi-line chart needs a range dimension per brushed value so
  that brushing one cell still shows its effect on the others. A selection tool or gate in a cell
  of a facet grid has the same question: does it filter the other cells? Each answer adds
  dimensions, and a `RangeDimension` starts a worker, so they should be made on first use.
- **Facets and filters interact.** In a category facet, filtering out a category elsewhere could
  empty a cell or remove it. Removing it reflows the grid under the user while they brush;
  keeping it shows empty cells. The compound dot plot keeps an explicit `split_values` list.
- **Persisted per-cell state.** Brushes keyed by field name or category value outlive the cell
  when the field is removed or the category renamed. Prune on change, as the multi-line chart
  does, and decide what a stale key means on load.
- **Scale.** Virtualisation only covers drawing. Loading every cell's column up front (#693),
  one worker per cell's compute, and one dimension per brushed cell all grow with the cell count.
  A grid over 1,000 genes needs lazy loading first.
- **Image export.** `buildChartSvg` copies the chart's SVGs at their on-screen positions, and a
  grid only mounts the cells in view. The multi-line chart's export draws every cell first, 16 at
  a time (each starts a density worker), then exports the grid's whole layout. Once cells load
  their data lazily (#693), export would also load every column: a grid over 1,000 genes needs a
  limit, a warning or a choice of which cells to export before that becomes the default.
- **Legends and axes.** One legend per chart works when every cell has the same categories and
  colours. A field split coloured by field (the multi-line chart's `by: "categories"` mode)
  needs a field legend instead. Sharing an axis only makes sense when the units are comparable,
  so it is an option, not a default (`share_x`).
- **Spatial layers per cell.** The density grid draws the render stack once a frame for all
  cells, with image layers that name no cell, so tiles load wherever they are drawn. A grid over
  channels or segmentations draws something different in each cell, so per-cell layers lose that
  sharing: each cell loads its own tiles or geometry, and memory grows with the cell count.
  Channel settings (contrast limits, colour) live in render stack entry props, which the channel
  controls write in place and immediately; a per-cell channel needs per-cell props there, not a
  second copy of the state, and must not bump `renderStackGeneration`. Shapes geometry must load
  once per element however many cells show it.
- **Selections across segmentations.** A selection in a nuclei cell and one in a cell-boundary
  cell pick different shapes, keyed to different tables. Whether one filters the other depends on
  the relations between elements (#555); until those exist, such a grid may have to keep
  selections per cell.
- **Layout state.** The grid sits inside a chart that gridstack sizes. Anything that sizes the
  chart div directly (`BaseChart.setSize(w, h)`) bypasses gridstack and leaves the chart
  overflowing its slot; cells must take their size from the layout, not the other way round.

## Display options (follow-up, #704)

Every grid takes `ChartArrayLayout`'s defaults today: cells at least 260px wide, square unless
they fill the panel, as many columns as fit, and one scrolling page. Wanted, in the `grid` config
(where `ChartGridConfigSchema` already says layout options belong), shared by every chart:

- **Element sizes**: a cell size (or minimum width) and aspect ratio, and/or a fixed number of
  columns, so a dozen marker cells can match the size of the single charts they replace.
- **Pagination**: a page of cells at a time, as an alternative to scrolling. It bounds how many
  cells are mounted, loaded (#693) and computed, which matters for a grid over hundreds of
  values, and gives export a natural unit (this page, or every page).

Both change how many cells are in view, so `visibleCellIndices`, the lazy loading window and
export should all follow the page rather than the scroll position when paging is on.

## Suggested order

1. Category facets for the multi-line chart (`grid: { by: "categories" }`), within that chart,
   to find out what a facet needs (keys from data, filters emptying cells, a field legend).
2. Move the compound dot plot's `split_by` panels onto `ChartArrayLayout`, keeping its config.
   Two charts sharing the layout for facets would show which parts are common.
3. A grid over image channels in the SpatialData chart, alongside its density grid: the first
   layer split, and the one most likely to be driven live (a protein picker).
4. Only then extract a contract, from several working cases rather than one.

Lazy loading of cells (#693) is independent and helps every case.
