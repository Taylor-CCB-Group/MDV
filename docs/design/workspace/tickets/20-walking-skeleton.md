Type: prototype
Status: open
Blocked by: [19 Layout engine bake-off and choice](19-layout-engine-bake-off-and-choice.md)
Map: [Workspace in MDV - Map](../README.md)
Merges: [17 Acceptance checks and the test set](merged/17-acceptance-checks-and-the-test-set.md)

## Question

Does the workspace hold up on real projects once it's built on the decisions from [18 One view, saved layout and screen sizes](18-one-view-saved-layout-and-screen-sizes.md) and [19 Layout engine bake-off and choice](19-layout-engine-bake-off-and-choice.md)? And what must every later slice pass before it ships?

Build the thinnest real version, on a branch and behind the Workspace preview preference (ADR 0026): one view, charts from every data source in the chosen engine, pane views converted on open, saving to and reloading from `views_v2.json`, and working at every supported width. Keep the existing chart chrome. The inspector, tool rail and filter bar wait for [21 Chart chrome, inspector, tools and filters](21-chart-chrome-inspector-tools-and-filters.md).

While it runs:

- the human opens real projects and says what feels wrong
- settle the test set: which real views and projects, and which automatic checks run on them (no overlap, nothing past the edge at each width, layout saved and reloaded unchanged, limits on redraws and frame rate, a ceiling on chart count for large projects)
- write those checks against the skeleton so later slices inherit them

The answer records what held, what didn't, and the agreed acceptance checks. The branch is the first build slice, ready to continue after `/to-issues`.

## Context

- This is the one build step this map carries; see the map's Notes.
- The last attempt went wrong at exactly this point: the prototype looked fine, then the real charts ended up all over the place. The skeleton exists to catch that before more is built on it.
- The prototype's own check (8 charts after a drag and a resize: none outside the grid, none overlapping) is the starting shape for the checks.
- ADR 0019 applies to anything extracted from legacy JS.
- Engine (ADR 0027, from [19 Layout engine bake-off and choice](19-layout-engine-bake-off-and-choice.md)): the in-house engine, built behind the engine seam so gridstack 14 can replace it. Start from the layout lab on this branch (`src/react/prototypes/layoutLab`: `engine.ts`, `InHouseMode.tsx`, `ChartHost.tsx`); build from what helps, delete the lab before merging, and carry over the hosting rules the ADR lists. Write the acceptance checks against `data-chart-id` and roles, never engine class names, so they can judge a gridstack 14 engine unchanged. The team's feel for the skeleton is the first test of the fallback.
