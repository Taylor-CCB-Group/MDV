# Workspace redesign: handover

Read this first. It applies to a person or an agent picking up the workspace redesign. It says what we're building, what's decided, how we work, and what's next. Everything it points to is on this branch, `docs/workspace-redesign-plan` (PR #719).

## What we're building

MDV today shows a project as one split pane per data source, each with its own menu bar and layout. The [workspace prototype](prototype.html) (open the file in a browser) shows one view per project instead:
- one 12-column grid of charts from every data source
- a docked inspector in place of the floating settings dialog
- a filter bar, a tool rail, view tabs
- light and dark themes

The team liked it. The goal now is a **showcase**: real MDV, with real projects, in the prototype's look and feel, so the team can judge whether this is the right direction. It is not a full conversion of MDV. Plan: [showcase-plan.md](showcase-plan.md).

**The prototype is the reference.** For any look-and-feel question, do what the prototype does, and only ask about what it doesn't show.

## Where things are

| What | Where |
|---|---|
| The plan: map, decisions, route | [README.md](README.md) |
| Tickets, with their questions and answers | [tickets/](tickets/), and folded-in ones in [tickets/merged/](tickets/merged/) |
| Research behind the decisions | [research/](research/) |
| The showcase plan (proposed, not yet agreed by the team) | [showcase-plan.md](showcase-plan.md) |
| The prototype | [prototype.html](prototype.html) |
| Decisions | [ADR 0026](../../adr/0026-preview-one-view-workspace-saved-to-views-v2.md): the preview setting, one view, `views_v2.json`, screen sizes. [ADR 0027](../../adr/0027-build-the-workspace-grid-on-the-prototypes-own-engine.md): the in-house layout engine behind a seam, gridstack 14 as the fallback |
| Terms | `CONTEXT.md`, under "Workspace terms" |
| Repo rules for agents | `AGENTS.md` at the repo root |
| The layout lab (reference code) | `src/react/prototypes/layoutLab/` on this branch |

## Where it stands

- **Resolved:**
  - 01: layout engine candidates
  - 02: everything that reads or writes view layouts
  - 03: chart chrome inventory
  - 04: where colours and light/dark are set
  - 18: one view, saved layout and screen sizes (ADR 0026)
  - 19: layout engine choice (ADR 0027)
- **Proposed:** the showcase plan. It replaces the old route (walking skeleton, then long grilling on tickets 21 and 22, then a PRD) with a one-week showcase. If the team agrees, apply its last section, "Changes to the map if this plan is accepted", to the map and tickets.
- **Open questions for the team, best answered with the showcase in front of them:**
  - the saved format
  - the inspector for charts that have their own dialogs (spatial, Viv, genome browser, table, gates)
  - how filters are shown
  - MDV-wide theming
  - undo
  - when, and whether, the preview merges and later becomes the default

## How we work

### Branches

- **`docs/workspace-redesign-plan` (PR #719)** is the one working branch, for both decisions and code. Nothing merges to main until the team has tried the showcase and we have iterated on it.
- **Day 1 of the showcase**, the grid with real charts, is built directly on this branch, because everything else builds on it.
- **After that, parallel work** goes on short-lived branches off it, each with a small PR back into it. That's where colleagues review each piece. People editing together live can commit straight to it.
- **The layout lab** (`src/react/prototypes/layoutLab/`, opened with `?layoutLab=`) is throwaway code on this branch. Build the showcase from it, then delete the folder, its `project_bootstrap.tsx` hook and the `gridstack14` package before anything merges.
- **When the team says yes,** this branch merges to main behind the Workspace preview, with its docs and ADRs.

### Decisions

- **The tracker** is the README's map plus `tickets/`, and its rules are in the README's Notes. In short:
  - a ticket is `open`, `claimed`, `resolved` or `merged`
  - claim it before starting
  - take the lowest open ticket whose blockers are resolved
  - record an `## Answer`, then update **Decisions so far** and **Route**
- **Grilling:** the agent asks one question at a time, each with its recommended answer and the reason. The humans decide, and the agent looks facts up in the code instead of asking.
- **Every decision becomes an ADR** in `docs/adr/` (next free number), and new terms go into `CONTEXT.md`.
- **To change a decision,** don't edit the old ADR or answer. Write a new ADR that supersedes it, mark the old one `Status: superseded by …`, and build the change as its own piece. The README's Notes list which decisions are cheap to change and which aren't.
- **ADR 0019 applies:** legacy UI is extracted into React in place, and restyling is a separate change.

### Code rules that matter most here

- **The engine seam (ADR 0027).** The grid takes the charts and the 12-column layout, and reports the whole new layout after each drag or resize. Keep everything else outside it, and keep engine class names out of other code and tests (find cards by `data-chart-id`), so gridstack 14 could replace it.
- **Hosting real charts.** The lab found these, and the showcase must keep each one (details in ADR 0027):
  - move the chart's own div into the card
  - remove the pane workspace's drag and resize from hosted charts
  - give the host `isolation: isolate`
  - mount the grid only after the view's layout has loaded
  - never redraw a chart while its card is animating, only once after a gesture
- **`AGENTS.md`:**
  - `pnpm exec tsgo` when TypeScript changes
  - the Python tests and pyright
  - no `as` casts where avoidable
  - plain, concise commit messages
  - read `docs/playwright/PLAYWRIGHT_AGENT.md` before any Playwright work
  - the spatial rules: measure with `perf.ts` and switch it off afterwards, and don't debounce image channel controls
- **Public repo.** Project names and data from anyone's dev server stay out of commits, including test scripts and docs.

## Running it

- **Front end:** `pnpm install`, then `pnpm run dev` (Vite on port 5170, or set `PORT=…`).
- **Backend:** the Flask server is expected at `http://127.0.0.1:5055` (`flaskURL` in `vite.config.mts`). Run it the way you normally do, for example in Docker.
- **The layout lab:** open a project with `?layoutLab=inhouse` or `?layoutLab=gridstack` on the end of its URL. It works in dev builds only.
- **Lab checks** (in `src/react/prototypes/layoutLab/checks/`):
  - `bakeoff.mjs`: layout rules at four widths, a drag, a resize, save and reload
  - `glitches.mjs`: the four glitches found by hand

  Both take project ids as arguments rather than in the file. Day 5 of the showcase points them at the workspace.

### Test projects

Three kinds of project caught the most problems:
- one with spatial or WebGL charts
- one view with many charts (23)
- one view over three data sources mixing pane layouts

Ask the plan's owner for the actual projects on the shared dev server. Their names stay out of the repo.

## Starting an agent on this

> Read `AGENTS.md`, `docs/design/workspace/HANDOVER.md`, the map at `docs/design/workspace/README.md`, ADRs 0026 and 0027, and `docs/design/workspace/showcase-plan.md`. Open `docs/design/workspace/prototype.html`: it is the reference. Then explain in your own words what we're building, what's decided and what's next, before changing anything.

If the explanation it gives back is wrong, fix this file before going on.
