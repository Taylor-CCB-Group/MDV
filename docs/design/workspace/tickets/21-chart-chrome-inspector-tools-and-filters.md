Type: grilling
Status: open
Blocked by: [20 Walking skeleton](20-walking-skeleton.md)
Map: [Workspace in MDV - Map](../README.md)
Merges: [09 Extension points for plugins](merged/09-extension-points-for-plugins.md), [10 What the inspector shows](merged/10-what-the-inspector-shows.md), [11 Global tool rail vs chart-owned tools](merged/11-global-tool-rail-vs-chart-owned-tools.md), [12 How charts describe their filters](merged/12-how-charts-describe-their-filters.md)

## Question

With the walking skeleton running, settle everything around a chart:

**Inspector**
- what does the docked inspector show, and when: the chart's `getSettings()` list, custom panels such as the spatial layer panel, several selected charts?
- popped-out charts have no inspector, so what do they get?
- does the floating settings dialog remain anywhere?
- what are the Chart, Charts and Data tabs each responsible for?

**Tool rail**
- which selection and navigation tools move to one rail (Select, Lasso, Pan, Zoom; Polygon and Modify?), and which chart types respond to them?
- what happens to tools only some charts have, such as spatial layer editing?
- how does someone know which tool is active?

**Filters**
- what must every chart type provide so the filter bar can show and clear its filter: a summary in words, its data source, the rows it keeps, a way to clear it?
- how are filters from linked data sources and colour-legend filters shown?
- what happens for chart types that can't describe their filter yet?

**Plugin extension points**
- what replaces `ChartManager.addMenuIcon` (per data source) and `chart.addMenuIcon` (per chart)?
- where do extra actions appear: the ⋯ menu, the toolbar or the menus?
- does the old API keep working?

**Annotations**
- each pane menu bar had an annotations action; where does it go in one view, given the prototype has no place for it? (from [18 One view, saved layout and screen sizes](18-one-view-saved-layout-and-screen-sizes.md))

**Menus, toolbar, ⌘K, shortcuts and view tabs**
- what goes in each, and which current actions are dropped?

## Context

- Findings: [03 Chart chrome inventory - findings](../research/03-chart-chrome-inventory-findings.md). Six chart types add their own icons; spatial, Viv, genome browser, table and gate charts have custom dialogs the settings renderer can't host; there are about a dozen filter mechanisms; `ChartManager.addMenuIcon`'s only live caller is GridstackManager (`hyperion.js` is never imported).
- ADR 0019 rejected replacing `addMenuIcon` with React state while callers keep the returned element.
- Settled before this map: the inspector replaces the floating settings dialog.
- Four topics in one session: if it runs long, finish the topic in hand, record it, and leave the rest of this ticket open.
