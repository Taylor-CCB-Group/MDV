Ticket: [01 Layout engine candidates](../tickets/01-layout-engine-candidates.md)
Map: [Workspace in MDV - Map](../README.md)
Sources: MDV repo at `main` `b3c5a868` (2026-11-25). Layout lab worktree `MDV-ui-refresh` at `ca488030` (2026-11-20). Workspace prototype `mdv-workspace.html` (scratchpad copy, 1,477 lines). Library facts were read on 2026-11-26 from npm (`npm view`), the GitHub repos and the GitHub API. Repo paths are relative to the repo root. Line numbers are for those commits.

Facts and trade-offs only. The choice is made in [06 Choose the layout engine](../tickets/merged/06-choose-the-layout-engine.md).

## Summary

- **gridstack** is now at 14.x. MDV is on 7.3.0 (April 2023), seven majors behind. The latest release at least two weeks old is 14.1.0 (2026-11-04). It has no dependencies, built-in touch, breakpoints, locking, auto-scroll and a `change` event. It has no keyboard move or resize. Upgrading MDV touches at least two options it uses today (`oneColumnSize`, `float`).
- **react-grid-layout** 2.3.0 (2026-10-05) is a TypeScript rewrite released in December 2025. It stores one layout per named breakpoint. It must host chart DOM through a React wrapper, as the layout lab did. Its generated narrower layouts do not scale widths, which is why the lab saw gaps. It has no keyboard support. A 3.0 with removals is already planned.
- **In-house** is about 110 dense lines in the prototype: `compact`, `flow`, `findSpot`, plus pointer handling. It does what the workspace design asks (12/6/1 columns, reading-order reflow, locked items that others flow around) and has Alt+arrow keyboard moves. It lacks resize from the west and north edges, robust touch, continuous auto-scroll, Esc to cancel, screen reader announcements, legacy DOM hosting and any tests.
- **"Unsaved Grid Stack View"** is explained by MDV's integration, not by gridstack. The view's baseline state is taken before MDV's observers finish writing chart sizes and positions, and MDV saves pixel sizes that depend on the window. Confidence: high that the cause is on MDV's side, medium on which of three mechanisms dominates. Not reproduced in a browser.

## Comparison

Versions compared: gridstack 14.1.0, react-grid-layout 2.3.0, in-house engine as it stands in the workspace prototype. MDV today runs gridstack 7.3.0.

| Criterion | gridstack 14.1.0 | react-grid-layout 2.3.0 | In-house (prototype) |
|---|---|---|---|
| Columns change with width | Yes. `columnOpts.breakpoints: [{w, c, layout}]`, or automatic `columnWidth`. Grid or window width. | Yes. `breakpoints` and `cols` maps on `Responsive`. Width from `useContainerWidth`. | Yes. Fixed rule: 12 from 960 px, 6 from 600 px, else 1. |
| Layout per breakpoint | Cached in memory per column count. `save()` writes the widest layout unless given a column count. No documented one-call save and load of all breakpoints. | First-class. `layouts` is `{breakpoint: Layout}`. `onLayoutChange(layout, allLayouts)` returns every stored one. | `l12` plus `custom[cols]`, stored only when someone rearranges at that width. Otherwise derived by `flow`. |
| How narrower layouts are made | Per breakpoint mode: `list`, `compact`, `moveScale`, `move`, `scale`, `none` or a callback. | Clone the nearest wider layout, clamp to bounds, compact. Widths are not scaled. | `flow`: reading order, widths scaled by `cols/12` with a per-type minimum. |
| Host legacy chart DOM | Yes, directly. `makeWidget(el)` on existing DOM. MDV does this today. | Only through React children. Legacy DOM is moved into a slot by a wrapper component (lab `ChartHost`). | Engine positions `.card` elements. A host slot would have to be written. |
| Size-settled callback | `resizestop`, `dragstop`, `change` (includes pushed items), `resizecontent`. Per-item DOM size needs a ResizeObserver. | `onResizeStop`, `onDragStop`, `onLayoutChange` (grid units). Per-item DOM size needs a ResizeObserver. | Redraw on pointer-up via `renderGrid` and `maybeRedraw`. No public callback. |
| Drag by handle | `handle` selector. | `dragConfig.handle` and `cancel`. 3 px threshold. | `.card-head` only, 4 px threshold. |
| Resize from edges | All eight via `resizable.handles`. Docs warn about top-side handles. | All eight via `resizeConfig.handles`. | East, south and south-east only. |
| Collision and compaction | `mode`: `top` (compact up), `float`, `list`, `compact`. | Pluggable compactor: vertical (default), horizontal, none, wrap, fast variants. `preventCollision`, `allowOverlap`. | Vertical compaction around pinned items. Drop onto a locked item is refused. |
| Lock one item | `locked` (not pushed by others) is separate from `noMove` and `noResize` (user cannot move or resize). | `static: true` (fixed, others flow around), or per-item `isDraggable` and `isResizable`. | `locked` flag. Not draggable or resizable, others flow around it. Ignored by `flow` on narrow layouts. |
| Touch | Built in since v6. Long-press to drag vs scroll since 13.3.0. Several touch fixes in 2025-2026. | Via react-draggable. `allowMobileScroll`, touch drop and edge auto-scroll since 2.3.0. | Pointer events with `touch-action: none` on the title bar and 9 px handles. Untested on devices. |
| Keyboard | Esc cancels a drag or resize, `r` rotates while dragging. No keyboard move or resize. PR #3092 open since 2025-07. | None documented. Issue #1907 open since 2023. | Alt+arrows move, Shift+Alt+arrows resize. Cards are focusable with `aria-label`. |
| Auto-scroll while dragging | Built in. Reworked in 12.5.0 and 14.1.0. | Built in since 2.3.0. | Only on pointer move near an edge. Stops when the pointer is still. |
| Bundle (min / gzip) | 97.3 KB / 26.6 KB (`gridstack-all.js`). CSS 6.1 KB. 7.3.0 was 77.7 KB / 20.5 KB. | 65.7 KB / 20.9 KB for `Responsive` with React external. CSS 2.7 KB. | About 5.1 KB / 2.2 KB for the extracted engine and pointer code. |
| Dependencies | None. | `react-draggable`, `react-resizable`, `fast-equals`, `clsx`, `prop-types`, `resize-observer-polyfill`. React 18+. | None. |
| Release cadence | 69 releases since 7.3.0. 27 of them in majors 12, 13 and 14, between 2025-04 and 2026-11. | 9 releases of 2.x since 2025-12-09. 2.4.0 unreleased. 3.0 planned. | Ours. |
| Open issues / PRs (2026-11-26) | 17 / 2. | 9 / 11. | n/a |
| Maintainers | One main author: 77 of the last 100 commits. | One main author: 92 of the last 100 commits. | MDV team. |
| Breaking changes that hit MDV | `oneColumnSize` removed (v10, shim removed 12.1). `float` replaced by `mode` (v14). ES2020 target (v8). CSS moved to variables (v12). | Whole new API in 2.0. Deprecations in 2.4 for removal in 3.0. | None external. |
| Known bugs that matter to MDV | Resize rounding drift with fractional `cellHeight` (fixed 12.5.0). Touch detection broke desktop drag in Chrome 144 and Firefox 147 (fixed 13.0.0). Several locked-item fixes. | Generated narrower layouts leave gaps (lab). Many 2.x fixes for drag-from-outside loops and double `onLayoutChange`. | `findSpot` never returns if `w > cols`. Callers clamp today. |

## gridstack

### Version and history

- Latest on npm: 14.1.1 (2026-11-17). Latest at least two weeks old: 14.1.0 (2026-11-04). Source: `npm view gridstack time`.
- MDV pins 7.3.0 in the lockfile (`pnpm-lock.yaml:4483`). `package.json:145` says `^7.1.1`.
- Releases since 7.3.0: 69, across majors 8 to 14. 27 of them are 12.x to 14.x, from 2025-04-13 on (`npm view gridstack time`).
- The README asks for donations to one person. The GitHub contributor list puts `adumesny` at 2,235 commits, and 77 of the last 100 commits are theirs.
- 9,265 stars. 17 open issues and 2 open PRs on 2026-11-26 (GitHub API search).

### Columns, breakpoints and per-breakpoint layouts

- v10 replaced `oneColumnSize` with `columnOpts`. Breakpoints are `{w, c, layout}` entries, measured against the grid or the window (`breakpointForWindow`). Source: [API.md, Breakpoint and Responsive](https://github.com/gridstack/gridstack.js/blob/master/doc/API.md), [README, Migrating to v10](https://github.com/gridstack/gridstack.js/blob/master/README.md#migrating-to-v10).
- Re-layout on column change is chosen per breakpoint: `list`, `compact`, `moveScale` (default), `move`, `scale`, `none`, or a function. `list` keeps items in sorted order, like the prototype's `flow`. Source: API.md, `ColumnOptions`.
- The engine caches one layout per column count in memory (`engine.cacheLayout`). `save()` writes the widest layout unless a column count is passed, and then only if that layout is cached. Source: API.md, `save()` and `cacheLayout()`. `save(columnCount)` arrived in 12.3.1 ([CHANGES.md](https://github.com/gridstack/gridstack.js/blob/master/doc/CHANGES.md)).
- The v10 migration notes say per-column layouts were "planned" and "TBD". I found no documented single call that saves and restores all breakpoint layouts. An app would call `save(…, column)` per breakpoint and `engine.cacheLayout()` on load. Not verified in code.
- MDV today stores one layout per chart: `gsposition` and `gssize` in the chart config (`src/charts/schemas/ChartConfigSchema.ts:51-61`). It sets `oneColumnSize: 400` (`src/charts/GridstackManager.ts:62`) and leaves `column` at the default 12. Line 52 computes a column count that is never used.

### Hosting legacy DOM and size callbacks

- `makeWidget(el)` turns existing DOM into a widget. API.md: "If you add elements to your grid by hand (or have some framework creating DOM), you have to tell gridstack afterwards". MDV wraps each chart div in two divs and calls `makeWidget` (`GridstackManager.ts:218-227`, `:298`).
- Events: `change` fires when widgets "change their position/size due to constraints or direct changes". Also `dragstop`, `resizestop`, `resizecontent`. Source: API.md, `on()`. MDV listens only to `resizestart` and `resizestop`, for a blur effect (`GridstackManager.ts:72-77`).
- In 7.3.0, widget positions are written as `gs-x`/`gs-y`/`gs-w`/`gs-h` attributes, not inline style (`node_modules/gridstack/dist/gridstack.js:1291-1305`). In 14.1.0 they are written as inline `top`/`left`/`width`/`height` with CSS variables, plus the attributes ([src/gridstack.ts:2024-2049 at v14.1.0](https://github.com/gridstack/gridstack.js/blob/v14.1.0/src/gridstack.ts#L2024)). This matters for MDV's MutationObserver (see the bug section).
- The grid watches its own size with a ResizeObserver in 14.1.0 (`src/gridstack.ts:2241`).
- v13 ships new React and Vue wrappers. 14.1.1 fixed the React wrapper unmounting item content when `change` fires mid-drag (#3459). MDV would not need the wrapper to host legacy DOM.

### Drag, resize, collision, locking

- Drag handle: `handle` selector. 14.0.0 finds handles inside open shadow roots. 14.1.1 makes the drop target follow the grabbed handle rather than the whole widget (#3463).
- Resize: `resizable.handles` takes any of the eight sides. API.md: "it is not recommended to resize from the top sides as weird side effect may occur". 13.2.0 shipped "better north and west resize handles" (#3360). MDV uses all eight (`GridstackManager.ts:63-65`).
- Packing: v14 replaced `float` with `mode: 'top' | 'float' | 'list' | 'compact'`. Old `float` values are auto-migrated at runtime, but the TypeScript field is gone ([README, Migrating to v14](https://github.com/gridstack/gridstack.js/blob/master/README.md#migrating-to-v14)). MDV passes `float: true` (`GridstackManager.ts:61`).
- Locking: `locked` "prevents being pushed by other widgets or api ... which is different from `noMove` (user action only)". Source: API.md, `GridStackWidget`. MDV's lock button sets only `locked` (`GridstackManager.ts:325`), so by the docs a "locked" chart can still be dragged and resized by the user. The lock state lives in a closure and is not saved in the chart config (`:318`). Not tested in a browser.
- Locked-item fixes since 7.3: #2503 drop onto a locked widget (10.1.2), #2981 docs clarified (11.5.1), #1959 `update()` refuses a move onto a locked item (14.1.0).

### Touch, keyboard, accessibility

- Touch is native since v6 ([README, Touch devices support](https://github.com/gridstack/gridstack.js/blob/master/README.md#touch-devices-support)). 13.3.0 added pause-to-drag versus scroll on mobile. 12.4.1, 12.4.2 and 14.1.0 fixed touch issues on Linux, hybrid devices and iOS auto-scroll.
- 13.0.0 fixed "drag broken in Firefox 147.0.4+ and Chrome 144+ due to `navigator.maxTouchPoints > 0` now returning true on desktop (macOS trackpad)" (#3231). MDV's 7.3.0 has the same check (`node_modules/gridstack/dist/dd-touch.js:14-22`). The issue report mentions drag-in from a sidebar. I did not check whether in-grid drag in MDV is affected.
- Keyboard: Esc cancels a move or resize and `r` rotates during a drag (10.2.0). There is no keyboard move or resize. PR [#3092 "Move items with keyboard controls"](https://github.com/gridstack/gridstack.js/pull/3092) has been open since 2025-07-10, last updated 2025-10-12.
- 13.1.0: "DOM matches the visual order when widgets are moved around (accessibility tabbing, regular printing)" (#3329).
- Auto-scroll: reworked in 12.5.0 (#2188, #2074) and 14.1.0 (#3356, #3355, #2666).

### Size and dependencies

- `dependencies: {}` (`npm view gridstack@14.1.0`).
- `dist/gridstack-all.js`: 97,265 bytes min, 26,607 bytes gzip -9 for 14.1.0. 77,717 / 20,487 for 7.3.0. File sizes from the jsDelivr package API, gzip measured locally on the CDN file.
- `dist/gridstack.min.css`: 6,138 bytes for 14.1.0, 8,586 for 7.3.0. v12 dropped `gridstack-extra.css`.

### Breaking changes from 7.3.0 to 14.1.0 that touch MDV's code

From [README, Migrating](https://github.com/gridstack/gridstack.js/blob/master/README.md#migrating) and CHANGES.md:

- v8: ES2020 package target. CSS rules renamed (`.gs-12`). `minWidth` removed.
- v9: `onParentResize()` renamed `onResize()`. MDV does not call it.
- v10: `oneColumnSize` and `disableOneColumnMode` removed. One-column switching is no longer on by default. 12.1 removed the conversion shim (#3027). MDV uses `oneColumnSize: 400`.
- v11: `addWidget` takes only a widget object. `renderCB` added. MDV uses `makeWidget`, which is unaffected.
- v12: dynamic stylesheet replaced by CSS variables. ES5 build removed. Nested-grid events go to the top grid.
- v13: breaking only for the Angular wrapper.
- v14: `float` option and `float()`/`getFloat()` replaced by `mode`.
- The methods MDV calls still exist in 14.1.0: `init`, `cellWidth`, `compact`, `destroy`, `makeWidget`, `removeWidget`, `update` (`src/gridstack.ts:92, 1002, 1059, 1169, 1311, 1445, 1610`). `compact()` now takes a layout mode argument. `el.gridstackNode` is still documented.

### Known bugs relevant to MDV

- #3241 resize from S/SE/E drifts the position through rounding, present since v4.2.1. Fixed in 12.5.0, which the changelog describes as "drifting rounding issue when cellHeight is small/fraction". MDV sets a fractional `cellHeight` (`GridstackManager.ts:55`).
- #2669 `load()` restores Y wrongly in some cases (fixed 10.2.1). MDV does not use `load()`.
- #3231 desktop drag broken by touch detection in new browsers (fixed 13.0.0). See above.

## react-grid-layout

### Version and history

- Latest: 2.3.0 (2026-10-05). A `legacy` dist-tag points at 1.5.4 (2026-07-29). Source: `npm view react-grid-layout`.
- 2.0.0 (2025-12-09) is "a complete TypeScript rewrite with a modernized hooks-based API". The v1 API lives on under `react-grid-layout/legacy` ([CHANGELOG 2.0.0](https://github.com/react-grid-layout/react-grid-layout/blob/master/CHANGELOG.md)).
- 2.x releases: 2.0.0, 2.1.0, 2.1.1, 2.2.0, 2.2.1, 2.2.2, 2.2.3, 2.2.4, 2.3.0. Six of them fall in December 2025.
- 2.4.0 is unreleased. It deprecates `droppingItem`, `onDropDragOver`, `data-grid`, `validateLayout` and several compactors "removed in 3.0". RFC 0003 "plans the 3.0 engine".
- The repo moved from `STRML/react-grid-layout` to `react-grid-layout/react-grid-layout`. npm metadata still names the old URL. 22,501 stars. 9 open issues and 11 open PRs on 2026-11-26. `STRML` made 92 of the last 100 commits.
- The layout lab already uses 2.3.0 (`MDV-ui-refresh/package.json:167`).

### Columns, breakpoints and per-breakpoint layouts

- `Responsive` takes `breakpoints` (`{lg: 1200, …}`), `cols` per breakpoint and `layouts` per breakpoint. `onLayoutChange(layout, layouts)` returns all of them. Source: [README, ResponsiveGridLayout Props](https://github.com/react-grid-layout/react-grid-layout/blob/master/README.md#responsivegridlayout-props). Features list: "Separate layouts per responsive breakpoint".
- Width is not measured for you since 2.0. Use `useContainerWidth` (README, Providing Grid Width). 2.2.4 fixed an update loop at fractional zoom. 2.3.0 added `debounceTimeout`.
- A breakpoint with no stored layout is generated by `findOrGenerateResponsiveLayout`: clone the nearest wider layout, `correctBounds`, then compact ([src/core/responsive.ts:110 at 2.3.0](https://github.com/react-grid-layout/react-grid-layout/blob/2.3.0/src/core/responsive.ts#L110), same code in the lab's `node_modules/react-grid-layout/dist/chunk-VBU6B7II.mjs:534-557`). Widths are not scaled. This is why the lab's narrower layouts had gaps until `scaleRects` scaled them (`MDV-ui-refresh/src/react/prototypes/layoutLab/GridMode.tsx:14-25`, `:38-44`).
- 2.4.0 (unreleased) bounds a generated copy to the target breakpoint's `cols` and moves items saved outside `cols` back inside, reporting it through `onLayoutChange`.

### Hosting legacy DOM and size callbacks

- Grid children must be React elements that forward `ref`, `style`, `className` and mouse and touch handlers (README, Custom Child Components). Legacy chart DOM must be moved into a slot by a wrapper.
- The lab does this in `ChartHost.tsx:12-46`: it appends `chart.getDiv()` into a slot, resizes through a ResizeObserver debounced 30 ms, and puts the div back on unmount.
- 2.3.0 memoizes `GridItem`. The README says children are compared by reference and must be memoized. A wrapper that re-creates children on every render will re-render every item.
- Callbacks are in grid units: `onDragStop`, `onResizeStop`, `onLayoutChange`. 2.4.0 fixes `onLayoutChange` firing twice per gesture (#1984). There is no per-item DOM size callback.

### Drag, resize, collision, locking

- `dragConfig`: `handle`, `cancel`, `threshold` (3 px), `bounded`. `resizeConfig.handles`: any of `s w e n sw nw se ne` (README, DragConfig and ResizeConfig). The lab used the chart title bar as the handle and seven resize handles (`GridMode.tsx:63-67`).
- Compaction is pluggable: `verticalCompactor` (default), `horizontalCompactor`, `noCompactor`, and in `/extras` `wrapCompactor`, `fastVerticalCompactor`, `fastHorizontalCompactor` (README, Compactor; CHANGELOG 2.1.0).
- Locking: `static: true` ("not draggable or resizable"), or per-item `isDraggable` and `isResizable` (README, Layout Item). Compaction treats static items as fixed (CHANGELOG 2.1.0, #2194).

### Touch, keyboard, accessibility

- Touch runs through react-draggable. 2.3.0 added `dragConfig.allowMobileScroll`, touch drop from outside, and `requestAnimationFrame` edge auto-scroll that "covers fixed-height grids and touch".
- No keyboard move or resize is documented. [#1795 "Accessibility support"](https://github.com/react-grid-layout/react-grid-layout/issues/1795) was closed on 2025-12-30 as a duplicate of [#1907](https://github.com/react-grid-layout/react-grid-layout/issues/1907), which is open and tracks keyboard access to resize handles.

### Size and dependencies

- Dependencies: `clsx`, `prop-types`, `fast-equals`, `react-draggable`, `react-resizable`, `resize-observer-polyfill`. Peer: React and React DOM. The README says 2.x needs React 18+. MDV uses React 19.2.6 (`package.json:163`). 2.4.0 drops `prop-types` and `resize-observer-polyfill`.
- Measured with esbuild 0.27.7 in the lab worktree, React external, importing `Responsive` and `useContainerWidth`: 65,667 bytes min, 20,895 gzip -9. Split: react-grid-layout 28.6 KB, react-draggable 19.1 KB, react-resizable 12.0 KB, fast-equals 3.6 KB, prop-types and clsx under 1 KB each.
- `react-grid-layout/core` alone (algorithms, no React): 12,766 bytes min, 4,608 gzip.
- `css/styles.css`: 2,709 bytes.

### Layout lab results

`MDV-ui-refresh/artifacts/playwright-cli/layout-lab/results.json` (2026-11-22, Apple M4, project 26 by default in `measure.mjs:9`):

- react-grid-layout grid at `lg`: 118-119 fps average, worst frame 25-33 ms, no long frames, 0-4 `setSize` calls per gesture.
- React Mosaic tiles at 1440 px: 119-123 fps, worst frame 9-32 ms. Divider drags cost 116-168 `setSize` calls each.

## In-house engine

### What the prototype does

Source: `mdv-workspace.html`, section "layout manager" (lines 662-709), geometry (729-735), pointer handling (958-1011).

- Model: every chart is `{i, x, y, w, h}` on a 12-column grid. The engine "is the only thing that positions or sizes a chart" (line 448-449).
- `compact(layout, cols, pinned)` (669-681): pinned items (locked, or the one being dragged) stay put. Every other item, in top-to-bottom then left-to-right order, is clamped to the column count and moved up until it hits something.
- `flow(l12, cols)` (685-693): builds a narrower layout from the 12-column one in reading order. Widths scale by `cols/12`, with a per-type minimum (`minW`, 682). Rows wrap when full. Order and rows stay intact.
- `findSpot(layout, w, h, cols)` (694): first free slot scanning rows top to bottom. Used for new and duplicated charts (1025, 1036) and by `tidy` (705-709) to close gaps.
- Breakpoints: `colsFor` picks 12 from 960 px, 6 from 600 px, else 1 (666). A narrower layout is stored in `custom[cols]` only after someone rearranges at that width (704). Adding or duplicating a chart clears all custom layouts (1026, 1037). "Reset to desktop order" deletes one (1369).
- Drag and resize (960-1004): pointer events on the title bar or a resize handle, 4 px start threshold. During a move the dragged item is pinned, others are compacted around it, and the preview is applied by CSS transforms. A target that overlaps a locked item is refused (990-991). The layout is stored on pointer-up.
- Redraw: charts redraw once on drop, not on every frame (`maybeRedraw`, 800; design note at 1328, citing the lab's 116-168 redraws per gesture).
- Locking: a locked card hides its resize handles (CSS line 192) and ignores drags (967). Others flow around it in `compact`.
- Keyboard: Alt+arrows move the selected chart by one cell, Shift+Alt+arrows resize it (`nudge`, 1005-1011, key handler 1453-1454). Cards have `tabIndex=0` and an `aria-label` (777).
- Undo and redo snapshot the whole state, including layouts (712-724).
- Size: the extracted engine, geometry and pointer code is 109 lines of dense code, about 5.1 KB min and 2.2 KB gzip (esbuild minify). No dependencies.

### What it would still need

- **Resize from every edge.** Only east, south and south-east exist (784, 985-987). MDV offers all eight today.
- **Hosting legacy DOM.** The prototype draws its own canvases. It would need a slot that adopts `chart.getDiv()` and calls `setSize` when the size settles, like the lab's `ChartHost`.
- **Touch.** Pointer events work in principle. There is no long-press to tell drag from scroll, and `touch-action: none` on the title bar blocks page scrolling from it. Edge handles are 9 px and the corner 18 px (CSS 188-190), below the WCAG 2.2 minimum target size of 24 px. No pointer capture. Never tested on a device.
- **Auto-scroll.** It scrolls 14 px per pointer move near an edge (979), so it stops when the pointer stops. It needs a timer or `requestAnimationFrame` loop.
- **Cancel.** No Esc to cancel a drag. `pointercancel` commits like `pointerup` (998-1003).
- **Accessibility.** Keyboard moves have no live-region announcement, no grab-and-drop mode and no instructions for screen readers. DOM order is insertion order, not visual order, so tab order does not follow the layout.
- **Locks on narrow layouts.** `flow` ignores `locked`, so a locked chart reflows at 6 and 1 columns.
- **Edge cases.** `findSpot` loops forever if `w > cols`. Callers clamp `w` today (699-700). No `maxW`/`maxH`, no RTL, no nested grids, no drag-in from outside.
- **Tests.** None. Compaction, flow, lock handling and breakpoint switching would each need unit tests, and drag and resize would need Playwright coverage.

## The "Unsaved Grid Stack View" bug

Report (Grid Stack and ChatMDV issues): after saving a gridstack view and coming back to it, the view shows unsaved, and the console logs "Unexpected unsaved changes shortly after changing view".

### Where the message comes from

- `src/charts/ViewManager.ts:253-267`. After `_init` resolves, `changeView` stores `getState()` as `lastSavedState`. A check 500 ms later logs the message if `hasUnsavedChanges()` is true. The comment there says "this happens when changing to a view with gridstack layout".
- `hasUnsavedChanges` deep-compares `currentState.view` with `lastSavedState.view` (`ViewManager.ts:493-538`). Each chart's config is in `view.initialCharts`, including `size`, `gsposition` and `gssize`.
- On `main` this check and the 1 s polling are commented out "until we fix state issues" (`ViewManager.ts:75-78`, `215-233`, `257-267`; `src/react/components/ViewSelectorComponent.tsx:71-92`). `hasUnsavedChanges` is still called by `RenameDatasourceDialog.tsx:170,180` and `useSlickGridReact.ts:807,1141`, so those paths can still see a false "unsaved" on gridstack views.

### How MDV writes gridstack layout

- `manageChart` (`GridstackManager.ts:164-422`) adds two observers per chart:
  - a ResizeObserver on the outer div, debounced 20 ms, that calls `chart.setSize()` and then saves `gsposition`/`gssize` from `outer.gridstackNode` (`:252-273`);
  - a MutationObserver on the outer div's `style` attribute only, debounced 20 ms, that saves `gsposition`/`gssize` (`:278-287`).
- It also calls `chart.deferredInit(() => chart.setSize(), 100)` (`:293-295`). The comment above it says the absolute-layout size "interferes with detecting unsaved changes" and "would still also have a false positive if window is resized".
- `setSize` writes the chart's pixel size into `config.size` (`BaseChart.ts:1273-1291`). `getState` keeps `config.size` for gridstack charts and only skips `position` (`ChartManager.js:1152-1159`).

### Three mechanisms, all in MDV

1. **The baseline is taken before the layout settles.** `addChart` calls `_makeChartRD` (which runs `manageChart`) and then awaits `chart.deferredInitsReady()` (`ChartManager.js:2039`, `:2082`). `deferredInit` runs its callback once synchronously inside a promise that resolves at once, and schedules the real call with `setTimeout` (`BaseChart.ts:396-407`). So `deferredInitsReady` resolves before the 100 ms `setSize` runs. The synchronous call in `manageChart` also runs before `makeWidget` and `update` have sized the widget (`:294` runs before `:298-307`). The ResizeObserver's first callback comes after observation starts ([Resize Observer spec](https://www.w3.org/TR/resize-observer/): observations start from a recorded size of 0×0) and is debounced 20 ms. So `size`, and possibly `gsposition`/`gssize`, change after `lastSavedState` is stored. `ChartManager.js:816-818` notes that the code relies on config being settled when these promises resolve.
2. **Pixel sizes that depend on the window and the previous view.** `config.size` is in pixels. `cellHeight` is one field on the manager, shared by every data source and kept across view changes. `getGrid` recomputes it from the current container and its own previous value (`GridstackManager.ts:37`, `:51-55`, comment "hack"). The same saved view can therefore produce different pixel sizes in a different window or after a different view.
3. **Pushed charts are not saved.** In 7.3.0 gridstack moves widgets by writing `gs-x`/`gs-y` attributes (`node_modules/gridstack/dist/gridstack.js:1291-1305`). MDV's MutationObserver only watches `style` (`GridstackManager.ts:285`). Only the dragged chart gets inline style during a drag (`dd-draggable.js:238`, `:278-279`). Charts pushed aside by a collision keep their old `gsposition` in config. MDV's own `compact()` saves positions by hand for this reason (`:88-106`). On reload, the stale positions overlap, gridstack resolves the overlap, and the first ResizeObserver callback writes the new positions into config after the baseline.

### Verdict

- The evidence points at MDV's integration. Nothing in it needs a gridstack bug. Each mechanism is visible in MDV's code, and two are acknowledged in MDV's own comments.
- Confidence: high that the cause is on MDV's side. Medium on which mechanism dominates in a given project. I did not reproduce it in a browser.
- Gridstack does offer what the integration lacks: a `change` event that covers pushed items, and `save()` in grid units. MDV uses neither. In 14.x, positions are also written as inline style, so the same MutationObserver would see pushed items. That would not fix mechanisms 1 and 2.
- A related risk, not verified: with `oneColumnSize: 400`, a data source pane narrower than 400 px switches gridstack to one column. The ResizeObserver would then save one-column `x`/`w` into each chart's config.

## Open questions for the bake-off

These are for [05 Layout engine bake-off on real projects](../tickets/merged/05-layout-engine-bake-off-on-real-projects.md) and [06 Choose the layout engine](../tickets/merged/06-choose-the-layout-engine.md).

1. Does gridstack 14 `list` mode, or `columnOpts` with `layout: 'list'`, reproduce the prototype's reading-order reflow at 6 and 1 columns without gaps?
2. Can gridstack save and restore a layout per breakpoint reliably with `save(…, column)` and `engine.cacheLayout()`? Is the result stable across reloads?
3. With react-grid-layout, is our own `scaleRects` (or a custom compactor) enough to avoid gaps at every breakpoint, including after a user edits a narrow layout?
4. Does hosting legacy chart DOM in react-grid-layout survive drag, resize, breakpoint changes and React re-renders without detaching or redrawing charts? Measure with `perf.ts` on spatial charts.
5. Does gridstack 14 `makeWidget` on the existing chart DOM keep WebGL contexts and spatial charts intact through drag and column changes?
6. Lock semantics: should a locked chart be immovable by the user (gridstack `noMove` plus `noResize`, RGL `static`) or only unpushable (gridstack `locked`)? Should locks hold on narrow layouts?
7. Touch on a real tablet: can a user scroll the page past charts and still drag by the title bar? Compare long-press behaviour.
8. Keyboard: none of the libraries moves or resizes by keyboard. How much of the prototype's Alt+arrow `nudge` carries over to each library, and does focus and tab order follow the visual order?
9. Save and reload: does each candidate round-trip a layout with no spurious "unsaved" state, once the baseline is taken after the layout settles?
10. Does saving store only grid units (no pixel `size`), so the layout is independent of window size? This overlaps [07 Saved layout format and compatibility](../tickets/merged/07-saved-layout-format-and-compatibility.md).
11. Upgrade cost if gridstack stays: replacing `oneColumnSize` and `float`, and checking the touch-detection fix in 13.0.0 against desktop drag in current Chrome and Firefox.
12. Maintenance risk: both libraries depend mostly on one maintainer. react-grid-layout has removals scheduled for 3.0. gridstack has shipped three majors in 2026.

## Could not verify

- The "Unsaved Grid Stack View" mechanisms were traced in code, not reproduced in a browser.
- Whether gridstack 7.3's touch detection breaks in-grid drag in current desktop Chrome or Firefox for MDV.
- Whether MDV's lock button lets a user drag a "locked" chart. The docs say it does.
- Per-breakpoint save and restore in gridstack 14 was read from the API docs only.
- The in-house engine size is an estimate from the extracted prototype code, not a production build.
- The workspace prototype's touch behaviour was not tried on a device.
