Type: grilling
Status: resolved
Blocked by: none
Map: [Workspace in MDV - Map](../README.md)
Merges: [07 Saved layout format and compatibility](merged/07-saved-layout-format-and-compatibility.md), [08 Merging data source panes into one view](merged/08-merging-data-source-panes-into-one-view.md), [14 Supported screen sizes](merged/14-supported-screen-sizes.md)

## Question

Settle the shape of a view's data, the part that is hardest to change later, in one conversation:

**One view instead of a pane per data source**
- how do existing multi-pane views convert: stack the panes in order, or keep a section per data source?
- what replaces each pane's menu bar actions: add chart, reset filters, colour settings, layout menu, annotations, full screen, row count?
- how does a chart show which data source it uses?
- what happens to linked data sources?

**Saved layout format**
- grid units and row height
- tablet and phone overrides, stored only when someone rearranges them
- a format version
- how existing absolute and gridstack views load
- what the Python writers (default views, template views, chart prototypes, spatial template) must write
- the engine is chosen afterwards in [19 Layout engine bake-off and choice](19-layout-engine-bake-off-and-choice.md), so does the format stay independent of any engine?

**Screen sizes**
- which sizes must the workspace support (desktop, laptop, tablet, phone)?
- what works on each: full editing, viewing and filtering, or read only?
- where do the breakpoints sit?

## Context

- Findings: [02 Everything that reads or writes view layouts - findings](../research/02-everything-that-reads-or-writes-view-layouts-findings.md). Today `gsposition` / `gssize` hold grid units for gridstack, `position` / `size` hold pixels for absolute, and `layout` is set per data source. 74% of the dev server's 197 views have more than one pane, and 14 views mix absolute and gridstack.
- Settled before this map: one view, with the data source chosen in Add chart.
- The prototype uses 12 columns from 960 px, 6 from 600 px and 1 below that.
- Record each decision as an ADR in `docs/adr/` and new terms in `CONTEXT.md`.

## Answer

Recorded as [ADR 0026](../../../adr/0026-preview-one-view-workspace-saved-to-views-v2.md), with four terms added to `CONTEXT.md` under Workspace terms. Decided across two chats on 2026-11-26.

- **An experimental preference.** A `preview` browser setting in `src/lib/browserSettings.ts`, listed while Feature previews is on, like Convert legacy charts. Off: the pane workspace as today. On: every view opens in the new workspace. The pane workspace stays until the team retires it.
- **Saved to `views_v2.json`, never to `views.json`.** Keyed by the same view names. With the preference on, a view opens from v2 if it has an entry, otherwise it is converted from v1 as it opens (an alert says saving keeps it). With it off, the v1 view opens unchanged. The two copies are independent from the first save. *(This reverses an earlier answer in the same grilling, "converted at runtime only, never saved".)*
- **v2 view shape:** `version: 2`; the view's `dataSources` map (keeps highlights); a flat `charts` list where each chart has `id`, `dataSource` and `config`; one `layout` map of chart id to `{x, y, w, h}` in 12 columns, row height fixed by the format version; `links`; `viewImage`. No pixels; positions are not in chart config.
- **Only the 12-column layout is saved.** The 6- and 1-column layouts are always derived from it in reading order. A tablet layout can be added later as a new key.
- **Other code during the preview:** only the new workspace writes v2; Python template writers and ChatMDV keep writing v1; `rename_view` and `delete_view` act on both files; `soft_delete_column` scans both (ADR 0003 fails closed); JSON history and later Icechunk history cover both; a v2-only view is listed only while the preference is on; the catalog thumbnail and static build stay on v1.
- **A "Copy layout JSON" action** so people trying the preview can send arrangements that worked or broke.
- **One flat grid, no sections.** Each chart card shows its data source.
- **Converted views use reading-order flow.** Pane by pane (old left-to-right order), then top to bottom and left to right within a pane; each keeps its width relative to its pane and is packed row by row into 12 columns with no gaps. Deterministic.
- **Pane menu bars are replaced as in the prototype.** Clear all filters in the filter bar; row counts per data source in the status bar; settings in the inspector; Add chart asks for the data source first. The absolute/gridstack layout menu goes. Annotations go to [21 Chart chrome, inspector, tools and filters](21-chart-chrome-inspector-tools-and-filters.md).
- **A view's data sources** are its `dataSources` keys when it opens (including a pane with no charts) plus any data source a chart is added for; none is dropped while the view is open. Saved highlights, filter tracking and the interaction link that adds an image chart keep working on this set.
- **Screen sizes:** 12 columns from 960 px (full editing); 6 from 600 px (everything except moving and resizing; inspector slides over the grid; added charts go at the end); 1 below 600 px (view, filter, select, full screen).
- **Treat the prototype as the reference** and only ask about what it doesn't show (the user's standing preference, saved as a memory by the other chat).
