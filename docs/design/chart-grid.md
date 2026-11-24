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

They answer the same questions (what is a cell, what does it share, where does its state live)
differently, and only the first two use `ChartArrayLayout`. A generic grid could make a new
chart's small multiples a matter of naming a prop, and make the layout, virtualisation and
settings the same everywhere. Facets by category (`grid: { by: "categories" }` for the
multi-line chart) are the next case.

## Kinds of cell

A prop can make cells in at least two ways, and they behave differently:

1. **Split by field.** A list of fields (`values`, `densityFields`). Every cell sees the same
   rows; each draws a different column. Cell keys are field names.
2. **Facet by category.** A category column (`split_by`, `categories`). Every cell draws the same
   columns over a different subset of rows. Cell keys are category values, which come from the
   data, can number in the hundreds, and change when the column is edited.

A third, a **sweep** of a setting (several band widths, several thresholds), is conceivable but
has no use case yet; leave it out.

Two props can make cells at once (the splatter plot's categories × fields). That is a 2D grid,
not a longer list, and wants rows and columns labelled as such.

## What a contract could look like

Not arbitrary props. Each chart type keeps an explicit list of the props its grid can use, as
`chartGridConfigSchema([...])` now does in the config schema, and for each says:

- **Kind**: field split or category facet, which decides how cells are enumerated and keyed.
- **How a cell narrows its inputs**: the field, or the row subset. A cell gets inputs, not a
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
- **Image export.** `buildChartSvg` copies the chart's SVGs at their on-screen positions. In a
  grid it exports the mounted cells only, which are the ones in view, at their scrolled
  positions. A grid export needs to lay out every cell, which means drawing the ones out of view.
- **Legends and axes.** One legend per chart works when every cell has the same categories and
  colours. A field split coloured by field (the multi-line chart's `by: "categories"` mode)
  needs a field legend instead. Sharing an axis only makes sense when the units are comparable,
  so it is an option, not a default (`share_x`).
- **Layout state.** The grid sits inside a chart that gridstack sizes. Anything that sizes the
  chart div directly (`BaseChart.setSize(w, h)`) bypasses gridstack and leaves the chart
  overflowing its slot; cells must take their size from the layout, not the other way round.

## Suggested order

1. Category facets for the multi-line chart (`grid: { by: "categories" }`), within that chart,
   to find out what a facet needs (keys from data, filters emptying cells, a field legend).
2. Move the compound dot plot's `split_by` panels onto `ChartArrayLayout`, keeping its config.
   Two charts sharing the layout for facets would show which parts are common.
3. Only then extract a contract, from three working cases rather than one.

Lazy loading of cells (#693) is independent and helps every case.
