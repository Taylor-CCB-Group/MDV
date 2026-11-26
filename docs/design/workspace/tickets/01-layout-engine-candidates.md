Type: research
Status: resolved
Blocked by: none
Map: [Workspace in MDV - Map](../README.md)

## Question

For gridstack (its current stable major), react-grid-layout 2.x and an in-house grid engine, what are the facts on:

- columns that change with the container's width (breakpoints), and how a layout per breakpoint is stored
- hosting DOM that legacy code builds: moving an existing chart `div` into a grid cell, and getting a callback when its size settles
- drag by a handle, resize from edges, collision and compaction, locking one item in place
- touch and keyboard support
- bundle size and dependencies
- maintenance: release cadence, open issues, breaking changes since gridstack 7.3
- known bugs that matter to MDV, including "Unsaved Grid Stack View" in Grid Stack and ChatMDV issues: is that gridstack's fault, or MDV's (`GridstackManager` writing `gsposition` from a MutationObserver)?

## Context

- MDV uses gridstack 7.3.0 through `src/charts/GridstackManager.ts`.
- The layout lab (worktree `MDV-ui-refresh`, `src/react/prototypes/layoutLab`) measured react-grid-layout and React Mosaic on a large real project: both held ~120 fps; react-grid-layout's own narrower layouts left gaps until we scaled them ourselves. Results: `MDV-ui-refresh/artifacts/playwright-cli/layout-lab/results.json`.
- The workspace prototype has an in-house engine of about 150 lines (`compact`, `flow`, `findSpot`).
- Facts and trade-offs only; the choice is made in [06 Choose the layout engine](merged/06-choose-the-layout-engine.md).

## Answer

- gridstack 14.1.0 (MDV runs 7.3.0, seven majors behind): no dependencies, built-in breakpoints, touch, locking and auto-scroll, hosts existing DOM via `makeWidget`. No keyboard move. Upgrading touches `oneColumnSize` and `float`. Per-breakpoint layouts are cached in memory, not saved as a set.
- react-grid-layout 2.3.0: stores one layout per breakpoint, but generated narrower layouts do not scale widths (the lab's gaps). Legacy chart DOM needs a React wrapper. No keyboard support. 3.0 removals already scheduled.
- In-house (~110 dense lines, ~2 KB gzip): matches the workspace design and has Alt+arrow moves, but lacks west/north resize, robust touch, continuous auto-scroll, Esc cancel, screen reader support, DOM hosting and tests.
- "Unsaved Grid Stack View" is MDV's integration, not gridstack: the baseline is taken before MDV's observers finish writing sizes, pixel sizes depend on the window, and charts pushed by collisions are never saved. High confidence on MDV's side; not reproduced in a browser.

Findings: [01 Layout engine candidates - findings](../research/01-layout-engine-candidates-findings.md)
