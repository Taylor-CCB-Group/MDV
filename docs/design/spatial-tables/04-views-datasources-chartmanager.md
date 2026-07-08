# 04 — Views, DataSources, and ChartManager

> Four related ideas, from cheapest to most radical: (B1) a datasource picker in "Add Chart",
> (B2) cross-datasource shared layout, (A) a new view kind that bypasses much of ChartManager,
> and (C) swapping the ChartManager at runtime onto new data. Plus the "radical vision" that ties
> them together.

## Current orchestration (grounded)

One `ChartManager` ([src/charts/ChartManager.js](../../../src/charts/ChartManager.js)) is a
de-facto singleton (`window.mdv.chartManager`, set at ChartManager.js:155; warns on a second
construction). It is built from **five inputs**:

```js
new ChartManager(holderDiv, datasources[], dataLoader, config, listener)   // project_bootstrap.tsx:232
```

- Constructor creates **one `DataStore` per datasource** (`new DataStore(d.size, d, dataLoader)`,
  ChartManager.js:183) and `this.dsIndex` (name→record).
- `_init(view)` (ChartManager.js:732) splits the content area into **one pane + one toolbar per
  datasource**: `splitPane(contentDiv, {number: dsToView.length, …})` (ChartManager.js:771), where
  `dsToView = Object.keys(viewData.dataSources)` (ChartManager.js:757). Each datasource toolbar's
  "Add Chart" opens `new AddChartDialogReact(dataStore)` with the target **fixed** by which toolbar
  was clicked (ChartManager.js:1626).
- A chart is bound to exactly one `DataStore` at construction:
  `new chartType.class(ds.dataStore, div, config)` (ChartManager.js:2268); registered as
  `charts[id] = {chart, dataSource: ds}` (ChartManager.js:2275); DOM parent is `ds.contentDiv`;
  gridstack is keyed **per DataSource** (`Map<DataSource, GridInstance>`,
  [GridstackManager.ts:39](../../../src/charts/GridstackManager.ts)).

### The serialized view shape

`views.json` is `viewName → View` ([ViewManager.ts:9](../../../src/charts/ViewManager.ts)):

```ts
type View = {
  dataSources:   Record<string, {panelWidth?, layout?, highlight?}>;  // keyed by ds NAME
  initialCharts: Record<string, ChartConfig[]>;                       // ds NAME -> configs
  links?; viewImage?;
};
```

**The critical fact: a chart config never names its datasource.** Its binding is *only* which
`initialCharts[dsName]` array it sits in. `getState()` re-materializes that binding from
`this.charts[id].dataSource` on save (ChartManager.js:1252). `addChart(dsName, config)` already
keys on a name.

### The chart ↔ host coupling surface

A chart needs very little from its host
([BaseChart.ts:118](../../../src/charts/BaseChart.ts)):
- Constructor `(dataStore, div, config)`. It holds `dataStore` + `config` + `div`; it does **not**
  hold a ChartManager reference.
- Real coupling is to the **`DataStore`** (columns, `Dimension`s for filter/highlight, color/legend
  helpers) — not to panes or toolbars.
- ChartManager injects menu icons and drag/resize *after* construction.
- Cross-datasource access is opt-in, pushed in via `setupLinks` / `_giveChartAccess`
  (ChartManager.js:2362) — the chart's *primary* dataStore stays singular.

The `ChartType` descriptor ([ChartTypes.ts:27](../../../src/charts/ChartTypes.ts)):
`{class, name, required?, params?, allow_user_add?}` drives the Add Chart dialog entirely.

## (B1) Datasource picker in Add Chart — low risk

`AddChartDialogReact` ([dialogs/AddChartDialogReact.tsx](../../../src/charts/dialogs/AddChartDialogReact.tsx))
currently takes a fixed `dataStore`; the type list is built from `Object.values(BaseChart.types)`
filtered by `allow_user_add` and `required` columns (AddChartDialogReact.tsx:262). On submit it
calls `window.mdv.chartManager.addChart(dataStore.name, config, true)`.

**Change:** add a datasource `<select>` at the top of the dialog. A ready-made
`DatasourceDropdown.tsx` ([dialogs/DatasourceDropdown.tsx](../../../src/charts/dialogs/DatasourceDropdown.tsx))
**already exists but is unused here**. On change, swap the `DataStoreContext` value and re-derive
the type list + column pickers; on submit call `addChart(selectedDsName, config)`.

**Because `addChart` already takes the ds name, no ChartManager change is needed for the picker
itself.** This is dialog-only and self-contained — the recommended first move here.

## (B2) Cross-datasource shared layout — moderate, structural

Today "one region per datasource" is baked into: `splitPane` (ChartManager.js:771), the pane loop
(776–823), `dsPanes`, `_setUpMenu` (per-ds toolbar), gridstack `Map<DataSource,…>`
(GridstackManager.ts:39,47), and `getState` pane-width readback (1233–1236). It is **structural,
not configurable**.

But nothing about a *chart* requires it: a chart only needs `ds.contentDiv` as a DOM parent and
`ds.dataStore` for data. Two charts with different DataStores can share one container. To share
space:

- Replace the pane loop with a **single content region + single toolbar + single gridstack**.
- Re-key gridstack on the **view/region** instead of per-DataSource.
- Update `getState` pane-width logic (currently requires `dsPanes[ds.name]`).
- The view shape's `dataSources[ds].panelWidth`/`layout` becomes **region-level**.

**Verdict:** structurally clean because chart↔datastore binding is already just a name, independent
of screen space — but it touches pane creation, gridstack keying, `getState`, and the view schema.
Pair it with (B1): one shared toolbar with a datasource picker in Add Chart is exactly the UX you
described (one small chart for a DS without it consuming a whole pane).

> **Schema note (this is a persisted `View` contract change).** Making `panelWidth`/`layout`
> region-scoped changes the shape of `views.json`. The intent is **additive / opt-in via a `kind`
> discriminator, not a migration**: existing per-datasource views keep their current shape *and*
> their current code path unchanged, while new region-scoped views take a discriminated branch —
> the two shapes co-exist (as [@xinaesthete noted on the PR](https://github.com/Taylor-CCB-Group/MDV/pull/522),
> behaviour for existing views doesn't change). So there is no load-time migration to write; the
> real requirement is that `getState()` and `ViewManager.hasUnsavedChanges` **dispatch on the
> discriminator** so they never serialize/compare a new-shape view with the old per-datasource
> assumptions (or vice-versa). Introduce an explicit `schemaVersion`/`kind` field now so the branch
> point is unambiguous rather than inferred. This is the same "audit every reader" caution as (A).

## (A) A new view kind that bypasses ChartManager/DataSources

**Cleanest seam:** branch on a discriminator in `_init(view)` (e.g. `view.kind === "custom"`).
In that branch, skip `splitPane` / per-ds toolbars / the `addChart` loop and render a custom
surface into `this.contentDiv`. `ViewManager.changeView` already tears down and re-inits
([ViewManager.ts:121](../../../src/charts/ViewManager.ts)), so entering/leaving is handled if the
branch mirrors the gridstack-destroy + `removeAllCharts` cleanup.

**What else must branch:**
- `getState()` (ChartManager.js:1224) assumes `dataSources`/`initialCharts`/`dsPanes` — needs a
  custom-view branch or the save / unsaved-changes machinery breaks.
- `ViewManager.hasUnsavedChanges` deep-compares the `view` object (ViewManager.ts:419).

**Hosting legacy charts inside it:** possible, but a legacy chart filters/highlights through a real
`DataStore` + `Dimension`, so you must supply a DataStore-compatible object — or accept that a
custom view only hosts new content. The cheaper path is a view whose `_init` branch renders a
single full-bleed surface **while still using real DataStores** (no pane/toolbar machinery), rather
than reimplementing filtering/linking.

**Verdict:** feasible via a discriminated-union view + branches in `_init`, `getState`,
`hasUnsavedChanges`, and the change/add-view flows. Main risk is the many implicit assumptions that
`viewData.dataSources` / `dsPanes` exist — audit every reader.

## (C) Runtime ChartManager swap onto new data, same views

The pieces a fresh ChartManager needs all funnel through the five constructor inputs, and `_init`
fully rebuilds the DOM from `viewData`. So "load the same view onto new data" is close to
"construct a new `ChartManager(holder, newDatasources, newDataLoader, config, listener)` then load
the same view JSON."

**Two strategies:**
1. **Construct a fresh ChartManager** (cleanest boundary). Must tear down the old one: unmount
   React portals (menu bar ChartManager.js:222, state handler in `project_bootstrap.tsx`), repoint
   `window.mdv.chartManager`, re-bind consumers. `project_bootstrap.tsx` already demonstrates the
   mount/unmount discipline for the state-handler portal.
2. **Keep the ChartManager, swap data behind DataStores** + `dataLoader`, then `_init(sameView)`.
   Lighter DOM churn but needs a reliable DataStore reset and a `dataLoader` swap API — **neither
   exists today** (`dataLoader` is captured per DataStore at construction).

**Blockers to design around:**
- The **`window.mdv.chartManager` singleton** and everything that reads it directly
  (`ViewManager` constructor ViewManager.ts:58, `AddChartDialogReact` AddChartDialogReact.tsx:332,
  gates, the menu-bar and state-handler portals).
- View JSON is **datasource-name + column keyed**, so "same view on new data" only works if the new
  data exposes **matching datasource names and the columns each saved chart references** — else
  `addChart` throws on unknown type / missing columns (ChartManager.js:1828).

**Verdict:** feasible; the "construct a fresh CM" route is least surprising given `_init` already
rebuilds from view JSON. Budget work for singleton/portal teardown and for validating new data
against the saved views' name/column expectations.

## The radical vision: datasources derived from spatial objects at runtime

Composing (A)+(C) with themes 1–2: `project_bootstrap` could synthesize **DataSource configs
directly from a SpatialData object's tables** (names + columns from `getObsColumnNames()` etc.),
build a ChartManager against a JS DataLoader (theme 1) reading those tables, and — when the
underlying spatial object changes — swap the ChartManager while keeping the same chart views (C),
because views bind by **name**. A new "spatial project" view kind (A) could own a lighter runtime
than the full pane/toolbar ChartManager for the charts that don't need it.

The enabling insight is the same one that makes (B1)/(C) cheap: **everything binds by datasource
name and column name.** If synthesized datasources present stable names/columns, existing views
and charts attach to them unchanged. The hard part is not the binding — it's (i) guaranteeing
stable synthesized names/columns across runtime changes, and (ii) the singleton/portal teardown for
the swap.

## Suggested order within this theme

1. **(B1) datasource picker** — dialog-only, immediate UX win.
2. **(B2) shared layout** — removes the per-DS screen segregation; pairs with (B1).
3. **(A) new view kind** — once (B2) has loosened the pane assumptions, a custom view branch is
   less surprising.
4. **(C) runtime swap** — last, gated on the runtime-derived-datasources vision being the
   confirmed direction; tackle the singleton coupling deliberately (consider dependency-injecting
   the ChartManager instead of `window.mdv.chartManager` as a prerequisite).
