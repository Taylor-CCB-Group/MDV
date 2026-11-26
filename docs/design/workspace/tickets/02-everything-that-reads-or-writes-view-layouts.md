Type: research
Status: resolved
Blocked by: none
Map: [Workspace in MDV - Map](../README.md)

## Question

List every place that creates, reads, changes or saves a view's layout or its data source panes, and count how saved projects use them.

- Front end: `ChartManager` (absolute vs gridstack, `changeLayout`, each data source's `contentDiv` and `menuBar`, `splitPane`), `GridstackManager`, `ViewManager`, `BaseChart`, `Popout`, the chart config schema (`gsposition`, `gssize`, `position`, `size`, `layout`).
- Python writers: `mdvproject.py` default views, `template_views/create_default_views.py`, `chart_prototypes.py`, `image_view_prototype.py`, `spatial/spatial_view_template.json`.
- Chat/LLM chart generation, tests and fixtures.
- Across the example projects and the projects on the local dev server: views using absolute vs gridstack, views with more than one data source pane, and typical chart counts per view.

## Context

Feeds [05 Layout engine bake-off on real projects](merged/05-layout-engine-bake-off-on-real-projects.md), [07 Saved layout format and compatibility](merged/07-saved-layout-format-and-compatibility.md), [08 Merging data source panes into one view](merged/08-merging-data-source-panes-into-one-view.md) and [17 Acceptance checks and the test set](merged/17-acceptance-checks-and-the-test-set.md).

## Answer

Layout is set per data source pane (`dataSources[ds].layout`), and a missing value means absolute. Each chart carries its own placement: `size`/`position` in pixels and `gssize`/`gsposition` in 12-column grid units, often both. `ChartManager._init` builds one Split.js pane per `dataSources` key. `GridstackManager` observers and `getState` write placement back, and `set_view` saves it as sent.
Python mostly leaves `layout` out, so those views load as absolute. `create_default_views`, `create_view_with_all_datasources`, `image_view_prototype` and ChatMDV write gridstack. The spatial template is mixed.
Dev server: 280 views in 37 projects (197 distinct). 74% of views have more than one pane. Panes are 379 absolute and 162 gridstack. 29 views mix the two layouts. Charts per view: min 0, median 7, max 23.

Findings: [02 Everything that reads or writes view layouts - findings](../research/02-everything-that-reads-or-writes-view-layouts-findings.md)
