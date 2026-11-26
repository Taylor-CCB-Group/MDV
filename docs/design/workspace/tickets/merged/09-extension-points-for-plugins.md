Type: grilling
Status: merged
Merged into: [21 Chart chrome, inspector, tools and filters](../21-chart-chrome-inspector-tools-and-filters.md)
Blocked by: [03 Chart chrome inventory](../03-chart-chrome-inventory.md), [08 Merging data source panes into one view](08-merging-data-source-panes-into-one-view.md)
Map: [Workspace in MDV - Map](../../README.md)

## Question

What do plugins and extensions get in the new chrome in place of `ChartManager.addMenuIcon` (per data source) and `chart.addMenuIcon` (per chart)? Where do extra actions appear (the ⋯ menu, the toolbar, the menus), does the old API keep working, and what do `hyperion.js` and `GridstackManager`-style callers need?

## Context

ADR 0019 rejected replacing `addMenuIcon` with React state while callers keep the returned element.
