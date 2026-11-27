# Workspace in MDV: planning map

The plan for bringing the [workspace prototype](prototype.html) (open the file in a browser) into MDV: one view per project instead of a pane per data source, a layout that fits the screen, and one colour token set with a stored light/dark setting. It is worked one decision ticket at a time; each grilling decision becomes an ADR in `docs/adr/`.

New here, person or agent? Read [HANDOVER.md](HANDOVER.md) first. Then start with **Decisions so far** and **Route** below. Open tickets are in [`tickets/`](tickets/), research findings in [`research/`](research/), and tickets folded into bigger ones in [`tickets/merged/`](tickets/merged/).

## Destination

A spec (PRD) for the new MDV workspace, an ADR for every decision that blocks building it, and an ordered set of build-slice issues ready to pick up. It covers the whole [workspace prototype](prototype.html) plus the shared colour tokens and light/dark setting.

## Notes

- **Tracker:** this folder, on the `workspace-showcase` branch. It is the only copy: edit it here, through commits or small PRs into that branch. Tickets are the notes in `tickets/`, numbered. Each has `Type:`, `Status:` (`open`, `claimed`, `resolved`, `merged`) and `Blocked by:` lines. A ticket is takeable when it is `open` and everything it is blocked by is `resolved`; the lowest number wins. Tickets folded into a bigger one sit in `tickets/merged/` with a `Merged into:` link.
- **Tickets are sized to a session, not to a single question.** Closely linked questions share one ticket so they are settled together. If a session runs long, finish the topic in hand, record it, and leave the rest of the ticket open.
- **One build step, the walking skeleton.** This map is planning, except for [20 Walking skeleton](tickets/20-walking-skeleton.md): the thinnest real version, on a branch and behind a flag, so the decisions are checked on real projects before the rest are made. Everything else produces decisions, not code.
- **Grilling tickets** are a conversation between an agent and the humans. The agent asks one question at a time, each with its recommended answer and why, and waits; the humans decide. The agent looks up facts in the code itself instead of asking. Each decision becomes an ADR in `docs/adr/`, and new terms go into `CONTEXT.md`.
- **Research tickets** are answered from primary sources (the code, library docs and source). Findings go in `research/` in this folder, linked from the ticket.
- **Working a ticket:** set its `Status:` to `claimed` before starting so nobody else takes it. When it is settled, add an `## Answer` section, set `Status: resolved`, add a one-line gist under **Decisions so far** linking the ticket, update **Route**, and turn anything newly clear into new tickets. Work one ticket per session.
- **ADR 0019** (`docs/adr/0019-extract-legacy-js-ui-into-react-in-place.md`) applies: legacy UI is extracted in place, and restyling is a separate change.
- **Must work on:** the hosted server (Flask/Docker) and the standalone static build.
- **Settled before this map:** a scrolling grid rather than tiles; one view, with the data source chosen in Add chart; the inspector replaces the floating settings dialog.
- **Proposed next step:** [Showcase plan](showcase-plan.md) (draft): a week-long, never-merged showcase of the prototype's look and feel on real projects, for the team to judge before the rest is planned.
- **References:** the [workspace prototype](prototype.html), the reference for every look-and-feel question: follow it, and only ask about what it doesn't show. Also the [MDV design system](https://claude.ai/code/artifact/a7e9d6e3-f9a7-4fef-b42b-948024ea8649) (a record of today's tokens), the layout lab on the `ui-refresh` branch, and the team's "Grid Stack and ChatMDV issues" notes.
- **Changing a decision.** Never edit a resolved ticket's answer or an accepted ADR to change it. Raise a new ticket (in `tickets/` while this map is open, a GitHub issue after it closes) saying what is changing and what the team found. Grill it into a new ADR that supersedes the old one, mark the old one `Status: superseded by <new ADR>`, and build the change as its own slice. Some decisions are cheaper to change than others:
  - Low cost:
    - the layout engine (ADR 0027's engine seam)
    - breakpoints and the tablet and phone rules (nothing narrower than 12 columns is saved)
    - chart chrome and tools
    - the preview setting
  - Highest cost: the saved format (ADR 0026). Once views are saved in `views_v2.json`, a new format needs a converter, keyed on the `version` field.
  - Hardest to undo: retiring the pane workspace and `views.json`. It is its own decision in [22 Theme, undo and rollout](tickets/22-theme-undo-and-rollout.md), taken only after the team has lived with the preview.
- **When every decision is resolved:** write the PRD, then cut it into build-slice issues.

## Decisions so far

<!-- one line per resolved ticket: ticket - gist of the answer -->

- [04 Where colours and light-dark are set today](tickets/04-where-colours-and-light-dark-are-set-today.md) - four colour layers that don't reference each other (legacy variables with a drifted catalog copy, shadcn, two MUI themes, ~470 hard-coded hex values); light/dark is an `<html>` class taken from the OS and never stored; the data palette has 104 colours. Unblocks [22 Theme, undo and rollout](tickets/22-theme-undo-and-rollout.md).
- [02 Everything that reads or writes view layouts](tickets/02-everything-that-reads-or-writes-view-layouts.md) - layout is per data source pane (missing = absolute); charts store pixel `size`/`position` and often grid `gssize`/`gsposition` too; saved via `getState` to `views.json` as sent. Dev server: 197 distinct views, 74% with more than one pane, 379 absolute vs 162 gridstack panes, 14 distinct views mix both, median 7 charts (max 23). Python writers split between absolute and gridstack. Unblocks [18 One view, saved layout and screen sizes](tickets/18-one-view-saved-layout-and-screen-sizes.md) and [20 Walking skeleton](tickets/20-walking-skeleton.md).
- [03 Chart chrome inventory](tickets/03-chart-chrome-inventory.md) - every chart shares one BaseChart title bar; only six types add icons; most of the 51 types are fully described by `getSettings()`, but spatial, Viv, genome browser, table and gate charts have custom dialogs the settings renderer can't host; SelectionOverlay serves only Deck scatter and the Viv/spatial region charts; about a dozen filter mechanisms, most describable in words, few saved in config; `ChartManager.addMenuIcon`'s only live caller is GridstackManager (hyperion.js is never imported). Unblocks [21 Chart chrome, inspector, tools and filters](tickets/21-chart-chrome-inspector-tools-and-filters.md).
- [01 Layout engine candidates](tickets/01-layout-engine-candidates.md) - gridstack 14.1 (MDV runs 7.3) has breakpoints, touch and locking built in but no keyboard moves and no documented save of narrow layouts; react-grid-layout 2.3 stores a layout per breakpoint but copies widths unscaled and needs a React wrapper for legacy charts; an in-house engine (~2 KB) already fits the design but lacks touch, west/north resize, screen-reader support, chart hosting and tests. "Unsaved Grid Stack View" is MDV's integration, not gridstack (early baseline, pixel sizes saved, pushed charts never written back). Unblocks [19 Layout engine bake-off and choice](tickets/19-layout-engine-bake-off-and-choice.md).
- [18 One view, saved layout and screen sizes](tickets/18-one-view-saved-layout-and-screen-sizes.md) - ADR 0026: the new workspace is a `preview` browser preference; pane views convert on open into one flat grid in reading order and save to a separate `views_v2.json` (flat chart list naming each data source, one 12-column `{x, y, w, h}` layout), leaving `views.json` untouched; a view's data sources are its `dataSources` keys plus any added; 12/6/1 columns at 960/600 px, tablet derived and no move/resize, phone view and filter. Unblocks [19 Layout engine bake-off and choice](tickets/19-layout-engine-bake-off-and-choice.md).
- [19 Layout engine bake-off and choice](tickets/19-layout-engine-bake-off-and-choice.md) - ADR 0027: the in-house (prototype) engine, behind a seam (charts + 12-column layout in, whole layout out per gesture; `flow`, `compact`, keyboard moves, chart host and saving outside it) so gridstack 14 can replace it if the team prefers it; gridstack felt the same by hand but can't grow a chart past the window bottom and needed its handles, column layouts and scrolling worked around; react-grid-layout was least smooth. Unblocks [20 Walking skeleton](tickets/20-walking-skeleton.md).

## Not yet specified

- Accessibility targets: keyboard use of the tools, screen-reader announcements, contrast in both themes. The grid's keyboard moves and 24 px handles are settled in [19 Layout engine bake-off and choice](tickets/19-layout-engine-bake-off-and-choice.md).
- The PRD and the build slices, once everything above is resolved.

## Out of scope

- The catalog page redesign: goes ahead separately.
- The desktop app (`mdv_desktop`) and embedding MDV in other pages: not required at first.

## Route

Restructured on 2026-11-26 from 13 open tickets into five sessions. Each session runs `/wayfinder` with this map and takes the lowest open, unblocked ticket.

1. ~~[18 One view, saved layout and screen sizes](tickets/18-one-view-saved-layout-and-screen-sizes.md)~~: resolved
2. ~~[19 Layout engine bake-off and choice](tickets/19-layout-engine-bake-off-and-choice.md)~~: resolved, in-house engine with gridstack 14 as the fallback
3. [20 Walking skeleton](tickets/20-walking-skeleton.md): the first real build slice, behind a flag, plus the acceptance checks
4. [21 Chart chrome, inspector, tools and filters](tickets/21-chart-chrome-inspector-tools-and-filters.md) and [22 Theme, undo and rollout](tickets/22-theme-undo-and-rollout.md): grilling with the skeleton running; either order
5. The PRD, then build-slice issues
