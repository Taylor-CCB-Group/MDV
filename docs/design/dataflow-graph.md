# A dataflow graph for MDV

Status: sketch (2026-11-17). Prompted by [cell-interactions.md](cell-interactions.md), but broader:
interaction sets are one kind of value among several that should flow between charts.

## Why

Charts in a view talk to each other through a handful of special-purpose mechanisms, each with its
own config shape and its own fixed topology. Two consequences:

- **One of everything per datastore.** One filter, one highlight, one active interaction set. A view
  can't hold two groups of charts that follow different selections, gene sets or pairs.
- **Things get built beside the model.** Where the model has no slot, a parallel mechanism appears
  (chart-scope filters as a main-thread React layer, `useCategoryFilterIndices` bypassing Dimensions,
  hyperion poking another chart's `category1/2`). intraspatial's
  `docs/mdv-dimension-vs-support-facet.md` §6 makes the same diagnosis from the filter side.

A graph whose nodes carry typed values and whose edges are explicit gives each of these a place, and
makes "more than one" the default rather than a special case.

## What already exists

### In MDV — a graph in all but name

| Mechanism | What flows | Implicit topology |
|---|---|---|
| `Dimension` filters (`DataStore`) | row selection | all dimensions AND-ed into one global filter per datastore; background filter as a 4-state byte |
| Highlight (`dataHighlighted`) | row indices | one per datastore, broadcast |
| Gates (`GateManager`, `config.gates`) | named selections | per datastore, persisted |
| Tags (`TagModel`, multitext column) | row labels (e.g. gene tags) | written by the annotation dialog; read as a filterable column |
| `links.rows_as_columns` → `RowsAsColsQuery` | a **live column set** (e.g. highlighted/filtered genes become columns on cells) | serialised *into a chart's config in place of a field name* (`FieldSpec = FieldName \| MultiColumnQuery`) |
| `links.valueset` | value set → filter on another datastore | fixed per link |
| `links.interactions` (`_addInteractionLinks`) | highlighted pair → filter on cells | fixed per link |
| `viewData.links` `chart_columnval_link` | a highlighted row's value → another chart's param / colour | per view, explicit source and targets |
| `view_state` link | pan / zoom | per view |
| `DomainSourceRef` (`src/spatial_stats`) | a region / gate / selection as a domain | per chart |
| SpatialData region selection | which region(s) a chart shows | per chart config |
| Contour `category1/2` | categories for density | per chart config |

`RowsAsColsQuery` is the closest thing to a graph value MDV has: serialised in config, resolved
lazily, live (it re-evaluates when the upstream highlight / filter changes), and usable wherever a
column is. A gene set by tag is already a three-step chain through it: tag column → filter on the
genes datasource → rows-as-columns query → `FieldSpecs` in a chart. Nothing names that chain or lets
a second one coexist in a view.

### In intraspatial — the substrate, on the compute side

- **Operation graph** (`src/gpu/graph/`, docs-site "Operation graphs"): typed ports, edges derived
  from reads/writes, lazy **pull** of a sink, topological execution, content-addressed **memo**
  (`op | params | inputKeys`) so a param change invalidates exactly its dependents. Serialisable; a
  React Flow **composer** and a fluent builder are two front doors to it. Bundled values with named
  parts (ADR-0023); feedback / unit-delay for iteration.
- **`support` facet** (ADR-0005): a filter is a field facet, soft mask ⇄ compact index; a filter is
  "just a subgraph producing a `support`". The reality check against MDV
  (`docs/mdv-dimension-vs-support-facet.md`) agrees on representation (MDV already ships the
  mask⇄index duality via `filteredIndexWorker`), argues for a filter **DAG** (gating trees are what
  users have in their heads), fixes the operator set (`min` / `max` / `1−a`, idempotent so diamonds
  are safe), and leaves open *who owns the graph*, *authored vs inferred*, and `filterSize`.
- **Expression IR** (ADR-0007, proposed): one typed IR with DSL text ⇄ graph ⇄ serialised surfaces,
  CPU-golden and GPU backends. The natural home for **columnar value expressions** (derived columns
  such as `log1p(CD8) − log1p(CD4)`, or a score over a gene set).
- `ParamSpec` metadata on parameters, which the composer already turns into controls.

## Proposal

### Evaluation: change propagates, work is on demand (runtime undecided)

Not decided yet. What is settled is the behaviour, which is neither pure push nor pure pull:

- **Change propagates.** When a source changes (a selection, a set, a parameter), everything derived
  from it learns it is stale, without the consumer polling.
- **Work happens on demand.** Nothing is computed unless something observing it needs it now. A TCM
  grid with twenty chosen pairs computes the few cells on screen, not twenty, and an off-screen chart
  computes nothing. This laziness, together with shared static layers drawn once into an offscreen
  framebuffer, is much of what keeps such views interactive. TCM being cheap per pair does the rest.
  Batching is an optimisation on top, not a requirement.
- **Results are cached by their inputs**, so two consumers asking the same question share the work,
  wherever they are in the graph.

Eager push, which recomputes every dependant on every change, gets the second point wrong. Bare pull,
which recomputes only when asked, gets the first wrong. MobX already has this shape: observables
invalidate, `computed`s are lazy and stay alive only while observed. The intraspatial op graph has the
second and third points (pull a sink, content-addressed memo) but no change notification: re-pulling is
the caller's job.

Options for where the graph lives:

- **A — two layers.** MDV holds a MobX signal layer for the values that change (selection, category
  set, column set, interaction set, region, focus, parameters). Heavy nodes (`support` masks,
  densities, TCM rasters, cross-PCF tensors, GPU expression columns) lower to op-graph pulls when
  observed. This is the least new machinery, but it gives two graphs to serialise and inspect.
- **B — one graph.** The intraspatial op graph gains subscriptions (mutable sources, observed sinks,
  invalidation along edges), and MDV uses it for everything. There is one serialised form and the
  composer can open a view's graph. But host-only values (category labels, column sets, focus) have to
  become op-graph citizens, and MDV's reactivity would depend on intraspatial's runtime.
- **C — no op graph in MDV.** MobX throughout, with compute nodes as async computeds calling
  intraspatial functions directly (as the cross-PCF chart does now). This is the simplest, but it gives
  up the memo's cross-consumer sharing unless MDV rebuilds it.

Whichever is chosen, the shared part is the **vocabulary**: value types, port typing, serialisation
shape, `ParamSpec`. Use intraspatial's where it exists rather than inventing MDV-local ones, and add
MDV value kinds (category set, column set) to the same typing scheme. That keeps A → B open as a later
move.

### Value kinds (port types)

```ts
type Value =
    | { kind: "selection"; datastore: string }                  // support: mask ⇄ index (ADR-0005)
    | { kind: "highlight"; datastore: string }                  // transient selection
    | { kind: "categories"; datastore: string; column: string } // a set of values of one column
    | { kind: "columns"; datastore: string }                    // FieldSpecs; RowsAsColsQuery is one producer
    | { kind: "expression"; datastore: string; dtype: DType }   // a derived column (ADR-0007 IR)
    | { kind: "domain" }                                        // DomainSourceRef: region, gate, selection, mask
    | { kind: "interactions"; datastore: string; column: string } // motifs × measures (cell-interactions.md)
    | { kind: "focus"; of: Value["kind"] }                      // transient emphasis of a value
    | { kind: "param"; spec: ParamSpec }                        // scalar / enum control
    | { kind: "field" };                                        // raster handle (density, TCM) — op-graph resident
```

Examples of the things this should carry, as nodes:

- **Dimension filters** — a `selection` source per filter, combined by `and` / `or` / `not` nodes
  (min / max / 1−a); today's global filter is the default conjunction of every chart's filter node.
- **Region selection** — a `domain` or `categories(region_field)` source that several spatial charts
  follow, instead of each holding its own.
- **Columnar value expressions** — `expression` nodes over columns (and over `columns`, e.g. a mean
  over a gene set), usable wherever a column is, as `RowsAsColsQuery` already is.
- **Sets of genes with a tag** — `categories(tag column)` on the genes datastore → `selection` →
  `columns` on the cells datastore via the rows-as-columns link: today's chain, made explicit and
  multipliable.
- **Interaction sets and their focus** — as in cell-interactions.md.

### Bindings, not copies

Chart config refers to a value by **reference**, never by holding a copy:

```ts
type ValueRef = { channel: "default" } | { node: string };
```

`{ channel: "default" }` resolves to today's implicit per-datastore value (the global filter, the
highlight, the default interaction set), so every existing chart is already "bound to the default
graph" with no migration. `RowsAsColsQuerySerialized` and `DomainSourceRef` are existing special
cases of a ref; they should become producers behind a ref rather than more one-off shapes.

### Persistence

- **Values** default to the datastore, as gates do (`dataStore.config.gates`): they survive view
  switches.
- **Topology** — extra nodes, and which chart port binds to which node — is per view, next to
  `viewData.links`. A view with no graph behaves exactly as today.
- Serialised graph = nodes `{id, op, params}` + port bindings, the same shape the intraspatial
  composer serialises, so the composer can open an MDV view's graph for debugging.

## First steps

No big bang; each step should be useful alone and falsifiable.

1. **Interaction sets through a channel binding** (cell-interactions.md step 1): the first value that
   is born a ref.
2. **A `ValueRef` type and resolver**, with `RowsAsColsQuery` and `DomainSourceRef` moved behind it.
   No behaviour change; proves the indirection on things that already exist.
3. **The background filter as scope → brush**, intraspatial's suggested first filter target: one
   concept that is currently spelled three ways, small enough that if the composition can't
   reproduce DotPlot / WGLScatterPlot behaviour, we find out cheaply.
4. **A second chart group in one view**: a view-level node (say, a second interaction set or a
   region selection) and charts bound to it. The first time "more than one" is real.
5. **Expression nodes** once ADR-0007's IR exists, starting with CPU evaluation into a virtual column.

## Open questions

- **Which runtime** (A / B / C above), and so who owns which graph. Values cross as GPU buffers or
  SAB-backed arrays either way (ADR-0005's zero-copy claim, still to test with one real selection).
  The first steps below don't depend on the answer: refs and resolvers sit in front of any of them.
- **Authored vs inferred.** Most users never see a graph: wiring comes from chart menus ("follow…",
  "new from here", "use this gate as domain"). An explicit editor — the composer, or something
  biologist-shaped for gating trees — is later and is not the same artefact as a developer's node
  editor.
- **Soft vs hard selections.** Does MDV's 4-state byte become two nodes (scope, brush) and does
  anything consume soft weights yet? (intraspatial note §8 q1.)
- **`filterSize` and other per-node reductions** — cheap on the GPU, a readback on the host.
- **Highlights**: row indices today. Do they become a `selection` like any other, or stay a separate,
  lighter channel?
- **Undo / history**: a graph of values is also a natural unit for undo; out of scope, but don't
  design it out.
