# Cell interactions as shared, live state

Status: proposal (2026-11-17). Follows the prototype cross-PCF chart (#645).

## Why

The cross-PCF chart keeps everything it knows — which type pairs matter, at which radii, over which
domain — in its own config. Nothing else can see it. We want that knowledge to drive other views:
the related cell types as density layers on a spatial chart, TCM layers for the chosen pairs, a grid
of TCMs (one per pair) once chart arrays land, and the existing interaction heatmaps / networks.

MDV already has an "interactions" feature. This note says what to keep from it, what to replace, and
the model that replaces it.

## What exists

### Legacy `interactions` (Hyperion era, barely changed since)

- **Data**: a separate datasource, one row per *(type 1, type 2, pivot)*, pivot being an ROI
  (`is_single_region`) or a condition. Statistic columns are arbitrary and precomputed outside MDV;
  the radius is baked into the column name (`gr20`, `"Cross PCF gr20"`). Written by
  `MDVProject.set_interactions` (`python/mdvtools/mdvproject.py:788`), declared by
  `InteractionsConfigSchema` (`src/charts/schemas/DataSourceSchema.ts:142`).
- **Charts**: `single_heat_map` (interaction matrix, one pivot value), `cell_network_chart`,
  `cell_radial_chart`, `single_series_chart` (one pair across ROIs). All d3/SVG.
- **Linking**: the only channel is a *row highlight* in the interactions datasource.
  `ChartManager._addInteractionLinks` (`src/charts/ChartManager.js:411`) turns it into a
  `filterCategories` on the cells datasource (both types, `operand: "and"` — flagged in-code as
  dubious) and opens an `image_scatter_plot` for the ROI. Colours sync via `sync_column_colors`.
- **Dead precedent**: `src/legacy/hyperion.js` once let the highlighted interaction row set
  `category1/category2` on centroid charts to draw density contours — the exact behaviour we want
  back, done by poking another chart's config.
- **Gaps**: no r axis, no domain beyond the pivot, no notion of "these pairs are selected" (a
  highlight is one row index), and the cell filter is destructive.

### Spatial density

- `useLegacyDualContour` (`src/react/contour_state.ts:603`): at most two categorical contours
  (`contourParameter` + `category1` / `category2`, each possibly several categories merged into one
  field) plus one per numeric `densityFields`. Config is flat on the chart.
- On WebGPU these render through `@intraspatial/deck` `DensityLayer` (one field per layer, categories
  by 0/1 weights; `src/webgl/intraspatialDensity.ts`).
- They react to their own chart's config only (a mobx `reaction` on `category1/2`).

### intraspatial (origin/main, 0.3.0)

- `crossPCFMatrixBinned` (all K×K×B at once, domains, edge correction) — what the chart uses.
- TCM: CPU reference, CPU bucket grid, GPU compute, and a fast two-pass WGSL render path
  (`renderTcm`, `src/gpu/spatial/tcmRender.ts`). Inputs are two clouds (A, B), `radius`, `sigma`,
  `alpha`, kernel, `bbox`. **No deck layer; renders on the library's own device; one TCM at a time
  (pooled textures under fixed keys); no domain support** (slice 4 of the domains plan).
- Docs already anticipate an "open-axis tensor field" `[a][b][r]` (`docs/muspan-cell-stats-plan.md`
  §7) and per-pair TCM on hover in < 1 ms (`docs/cell-stats.md` §8).

### Chart array FBO

PR #476 (on #459): one deck, one `OrthographicView` per grid cell, shared static layers drawn once
into an offscreen framebuffer and composited per cell. WebGL / deck 9.1 only, cells are numeric
density fields, and the PR calls itself a prototype to be re-implemented. The useful idea is the
*layer scoping* (`chart-shared` vs `per-viewport`), not the code.

## Decision sketch

**Keep the word, replace the mechanism.** "Interactions" stays the user-facing concept. The legacy
precomputed table becomes one *source* of interaction values among several; the shared state becomes
a first-class, observable object rather than a row highlight.

### The model

Three separable things, each persisted by label (never type id), as the chart does today:

1. **Type axis** — what a "type" is: `{ datasource, column }` (e.g. cells / `leiden`). Everything
   below is relative to one axis. A second axis (e.g. a niche column) is a second set, not a
   generalisation.

2. **Interaction set** — the shared selection: a plain, serialisable value. The model says nothing
   about what holds it; that is a separate question ([Where a set lives](#where-a-set-lives)).

   ```ts
   type TypeTerm = string | { anyOf: string[] };      // a category, or a merged group (as category1 arrays today)
   type Motif = {
       terms: TypeTerm[];                             // arity 1: a type; 2: A→B; n: higher-order
       ordered: boolean;                              // A→B is asymmetric for g and TCM
   };
   type InteractionSet = {
       id: string; name: string;
       axis: { datasource: string; column: string };
       motifs: Motif[];                               // the persisted selection (today's `pairs`)
       measures: MeasureRef[];                        // what is being compared (today's `r_bins`)
       domain: DomainSourceRef | null;                // already exists
   };
   ```

   **Focus** (hover / emphasis) is a second, transient value that travels with a set, not persisted:
   `{ motif?, measure? }`. The chart's current `focus` / `rFocus` become this.

3. **Measures** — how two (or n) types relate, with parameters as data:

   ```ts
   type MeasureRef =
       | { kind: "cross_pcf"; r: number | "mean" }        // reads a computed K×K×B result
       | { kind: "tcm"; r: number; sigma: number; alpha?: number }
       | { kind: "column"; datasource: string; column: string }  // a legacy precomputed stat, e.g. gr20
       | { kind: "contact"; r: number };                  // intraspatial contactNetwork, later
   ```

   Several r's are several measures of one kind. "Multiple measures between the active pairs" is
   just `motifs × measures`, which is exactly the chart's grid of matrices generalised.

**Results** are a cache, not state: keyed by `(axis, domain, measure kind + compute params)`, values
an open-axis tensor `[a][b][param]` (cross-PCF already returns `g[a·K·B + b·B + k]`). A consumer asks
for `(motif, measure)` and gets a value, a curve, or a raster handle. This is where the current chart
compute (worker, latest-wins, abort) moves to, so two charts asking the same question don't compute
twice. The key is the *inputs*, never a set's identity or host: two sets in different places that
ask the same thing share the result.

### Where a set lives

Three concerns that the first implementation will be tempted to merge, kept apart so the set can later
be composed differently:

- **Value** — the `InteractionSet` and its focus, as above. Pure data; no behaviour, no owner.
- **Host** — where a value is stored and persisted. Interchangeable; candidates, roughly in the order
  they're likely to appear:
  - *datastore default*: one set per datastore, persisted in its config like gates
    (`dataStore.config.gates`), shared by every chart that hasn't been told otherwise;
  - *chart-local*: a chart's own unshared set, persisted in its config (what the chart does today, and
    a legitimate choice, not just a legacy state);
  - *view node*: a set belonging to a view's graph, so one view can hold several chart groups on
    different sets;
  - *derived*: not stored at all, computed from other values (a union of sets, "pairs with g > t at
    r", a legacy row highlight read as a focus) — so read-only to consumers.
- **Resolver** — turns a chart's binding into the value and, if the host allows it, a way to write
  it back. The only thing consumers touch.

**First cut:** the datastore-default host, implemented in the GateManager mould
(`getOrCreateInteractionManager(dataStore)`), plus chart-local. The manager is an *implementation of
the default host*, not the model: nothing but the resolver may import it. What keeps later
composition open:

- consumers depend on the resolver and the value types only;
- a binding names *where* to find a set, so adding a host is a new binding kind, not a consumer change;
- writes go through the resolver, which may say "read-only" (derived) — consumers must handle that;
- focus belongs to the binding's channel, not to a global, so two chart groups have two foci;
- the compute cache keys on inputs (above), so it doesn't care which host a set came from;
- a set's axis names its datastore, so a host and a consumer needn't be on the same datastore (a set
  on cells can be read by a chart on the legacy interactions datasource).

### Who reads and writes it

| Consumer | Reads | Writes |
|---|---|---|
| Cross-PCF chart | set, focus, cross-PCF results | motifs (paint), measures (r bands), focus |
| Spatial chart — type density | motifs' types, focus | — |
| Spatial chart — TCM layer | focused / chosen ordered pairs, a `tcm` measure | — |
| TCM grid (chart array) | chosen pairs as cells | focus (hover a cell) |
| Legacy interaction heatmap / network | focus | focus (adapter, below) |

Charts hold a **binding**, not a copy and not a host handle:

```ts
type InteractionsBinding =
    | { channel: "default" }      // the datastore-default host (first cut)
    | { local: InteractionSet }   // chart-local
    | { node: string };           // a view-graph node, source or derived (later)
```

Creating a cross-PCF chart binds it to the default channel, creating the set if none exists; a second
chart on the same channel shares it, and "detach" turns a binding into `local` with a copy. Today's
chart config migrates by lifting `type_column`, `pairs`, `r_bins`, `domain_source` into a set — into
the default host, or `local` if another chart already owns the default and the two differ.

**Spatial layers bind to a set, not to categories.** A render-stack overlay entry gains a binding
such as `{ interactionSet, role: "a" | "b" | "all", follow: "focus" | "chosen" }`. The density
overlay then derives its categories (one `DensityLayer` per type, coloured by the type's own colour)
instead of reading `category1/2`. `category1/2` stay for charts with no binding. This avoids the
hyperion pattern of writing another chart's config.

**Focus is not a filter.** Hovering a pair emphasises its types elsewhere (dim others, draw their
densities); it doesn't filter the cells datasource. Filtering to a motif becomes an explicit action
("filter to these types"), unlike `_addInteractionLinks`, which filters on every highlight.

### Toward a dataflow graph

A longer-standing want: a dataflow graph carrying signals between charts — active interaction sets
among them, but also Dimension filters, region selection, columnar expressions, tagged gene sets — so
one view can hold several groups of charts, each following different values. The general design is
in [dataflow-graph.md](dataflow-graph.md); this section is what it means for interactions. The
interaction design should be the first tenant of that graph, not a parallel mechanism it later has to
absorb.

- **Signals, not stores.** An interaction set and its focus are *values on a channel*. The
  datastore-default host is then just the default graph: one writable source node per datastore,
  which every unbound chart reads; chart-local, view-node and derived hosts are other nodes. Nothing
  about a chart needs to know which it is bound to.
- **Node kinds** this design already implies:
  - *source* (writable state): an interaction set; focus.
  - *derived* (pure, cached): set → its types (for density); set ∪ set; "pairs with g > t at r";
    legacy row highlight → focus.
  - *compute* (async, cancellable): set + domain → cross-PCF tensor; pair + measure → TCM raster.
    The result cache keys on the inputs, so two branches asking the same question share the work.
- **Ports per consumer.** A chart declares typed inputs/outputs (the cross-PCF chart reads `set`,
  `focus`, `pcfResult`; writes `set`, `focus`). A spatial overlay's binding (`role`, `follow`) is an
  input port. Wiring a second chart group = creating a second source node and pointing their ports at
  it.
- **Existing edges this could subsume** (not in scope here, but the graph should be able to express
  them): `viewData.links` `chart_columnval_link` (a highlight → another chart's param), datasource
  `links` (`valueset`, `interactions`, `sync_column_colors`), `view_state` pan/zoom sync, gates and
  domain sources as values. Filters and highlights are the obvious first-class signals after this.
- **Persistence split.** Source *values* default to the datasource (like gates in
  `dataStore.config.gates`) so they survive view switches; *topology* — extra nodes and which chart
  binds to what — is per view, alongside `viewData.links`. A view with no graph behaves exactly as the
  default.
- **Reactivity substrate.** MobX is already how charts observe config; a node is a mobx observable
  (source) or `computed` (derived), with compute nodes as async computeds with abort, which the
  cross-PCF chart effect already is. No new reactive library needed to start.

What this asks of step 1 is small: route every read/write through the resolver
(`resolveInteractions(binding, dataStore)`), never the default host directly, and make focus a value
on the binding's channel rather than per-chart React state. The graph itself — UI to see and wire it,
more hosts and node kinds — can come later without migrating charts again.

### Legacy interop

- A legacy interactions datasource becomes a `column` measure source: given `(a, b, pivot)` it
  returns the row's value. The cross-PCF chart could then show `gr20` alongside live g.
- An adapter maps a legacy row highlight to `focus = {motif: [c1, c2]}` and back (the
  highlighted row for the focused pair and current pivot), so old heatmaps and networks join in
  without being rewritten. `_addInteractionLinks`' filter behaviour stays opt-in for old projects.
- Pivot ≈ domain: a legacy pivot value (ROI) maps to `{kind: "region", region}`. Pivots that are
  conditions (groups of ROIs) have no domain equivalent yet — note, not solve.

### n-types and things not enumerated

The motif's `terms` array is the extension point: arity 1 covers "show these types"; arity 2 covers
pairs; arity 3+ covers co-occurrence / neighbourhood motifs (intraspatial has none yet, only
proposals). `TypeTerm.anyOf` covers merged groups ("all T cells"). Wildcards (A→*, "A's
neighbourhood") are a later `{ any: true }` term. Measures stay a closed, tagged union so each
consumer can say which kinds it understands and ignore the rest.

## What intraspatial needs (to raise there)

For a TCM layer and a TCM grid, in order:

1. **TCM on a caller-supplied device** (as `splatDensityResident` already takes `device`), so it can
   run on deck's.
2. **Per-instance outputs**: no module-level pooled `"mark"` / `"gamma"` textures; the caller (or a
   layer) owns them, so several TCMs coexist.
3. **A `TcmLayer` in `@intraspatial/deck`**, mirroring `DensityLayer` (props: A / B weights over one
   point set, `radius`, `sigma`, `alpha`, kernel, diverging colour range). Weights over a shared point
   set, not two clouds, fits MDV's columns and makes per-pair layers cheap.
4. **Batched pairs**: one pass B-density per distinct B type, reused by every pair sharing it; an
   array-texture output for N pairs (the grid).
5. **TCM on domains** (slice 4): ρ_B·|disk∩W| in place of ρ_B·πr².
6. Already listed for cross-PCF: typed-array `LabelledCells`, `indices?`.

## Plan

Each step lands on its own and is usable.

1. **Set model, resolver, and two hosts** — datastore default and chart-local — migrating the
   cross-PCF chart (no new behaviour; two cross-PCF charts can now share a set, or not). Move the
   compute + cache into a service keyed by inputs. Two hosts from the start keeps the resolver honest:
   with only one, nothing would catch code that bypasses it.
2. **Type-density binding on the SpatialData chart**: an overlay bound to a set draws one density per
   type of the focused / chosen motifs. Uses existing `DensityLayer`; no intraspatial change.
3. **Legacy adapter**: legacy heatmap highlight ↔ focus; legacy stat columns as `column` measures.
4. **intraspatial TCM asks 1–3**, then a TCM overlay bound to a set (focused pair first).
5. **TCM grid**: re-implement the chart-array idea on WebGPU (deck views per cell + shared static
   pass), cells = chosen pairs; needs ask 4 for speed. Coordinate with whoever picks up #476.
6. **Dataflow graph, first cut**: per-view extra source nodes and chart bindings to them (two chart
   groups in one view on different sets), then derived nodes. See
   [dataflow-graph.md](dataflow-graph.md); this note only fixes the seam.

## Open questions

- Persistence is proposed as values at datasource level, topology per view (see the dataflow
  section). Does a view-local set (one that should *not* leak to other views) need to exist, or is a
  per-view node referencing a datasource-level value always enough?
- How a user sees and edits wiring: implicit (chart menu "follow set…", "new set from here") before
  any explicit graph UI?
- Does focus carry a measure (r) as well as a motif? The chart already needs both (rFocus); the
  spatial TCM layer would follow it.
- Do conditions (groups of ROIs) become a domain kind (union of regions, per-ROI averaged g), or stay
  a legacy-only pivot?
- Colour: pairs currently get a palette by position; types have their own colours. A TCM grid cell
  probably wants both (pair swatch, type-coloured labels).
