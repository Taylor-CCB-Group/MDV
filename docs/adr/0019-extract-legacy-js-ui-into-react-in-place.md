# Extract legacy JS UI into React one piece at a time, mounted in place

**Status:** proposed. Applied once, to the datasource menu bar: `ChartManager._setUpMenu` now
mounts `src/react/components/DataSourceMenuBar.tsx`. The main menu bar
(`MenuBarComponent.tsx`) was extracted the same way earlier. The change that extracted the
datasource menu bar also restyled it (`DataSourceMenuBar.css`), against step 6 below; the next
extraction should keep the two apart.

## Context

A lot of MDV's UI is built imperatively in JS: a method on `ChartManager`, `BaseChart` or a
dialog calls `createEl` / `createMenuIcon`, keeps the elements on `this` or on a datasource
record, and updates them from listeners. We want that code in TypeScript and React, but
`ChartManager.js` is too large and too widely depended on to convert in one go, and plugin-style
code (`GridstackManager`, `src/legacy/hyperion.js`) calls its imperative methods directly.

So the conversion happens a piece at a time, with the legacy owner still in charge of where and
when the piece appears. This records how to do one such extraction, using the datasource menu
bar as the worked example, so that the next one does not have to rediscover the same problems.

## Decision

Extract one DOM subtree at a time into a React component, mount it from the legacy code that
used to build that subtree, and delete the legacy code in the same change. The extraction
preserves behaviour and appearance; improvements come in separate changes.

### 1. Choose a piece that one method builds

A good piece is a subtree that a single method creates into a container it was handed, such as
`_setUpMenu(ds)` filling `ds.menuBar`. If building it is spread over several methods, extract
them together or not at all.

### 2. List what the legacy code shares before writing anything

Search for every field the method assigns on `this` or on a shared record, every listener it
registers, and every element it hands to something else. Each one is either state that only
this UI used, which moves into the component, or something other code depends on, which has to
keep working. For the datasource menu bar:

| Legacy code | Used by | What happened to it |
| --- | --- | --- |
| `cm.progressBars[name]`, updated in `_addDSListeners` | only the bar | moved into the component (`useFilterSize`) |
| `cm.isFullscreen` and a `fullscreenchange` listener | only the bar | moved into the component (`useIsFullscreen`) |
| `ds.filterBar` | nothing: it was created with a mistyped parent and never attached | deleted |
| `cm.layoutMenus[name]` | built elsewhere in `ChartManager` | still read from `ChartManager` |
| `ds.menuBar` | `cm.addMenuIcon(dsName, ...)`, called by `GridstackManager` and `hyperion.js` | kept as a legacy element (see 4) |

The search also turns up behaviour that depends on the shared state by accident. Here, the
`"filtered"` listener returned early when a datasource had no progress bar, which also meant
`ChartManager` listeners were never told about filters on a datasource that had never been in a
view. That is kept, under an honest name (`cm.dataSourcesShown`), and not silently changed.

### 3. Mount with `createMdvPortal`, and unmount where the legacy DOM is thrown away

```js
ds.menuBarRoot = createMdvPortal(createElement(DataSourceMenuBar, { dataSource: ds }), container, undefined, this);
```

`createMdvPortal` supplies the theme, the query client and the `ChartManager` context that
hooks such as `useChartManager()` need. Keep the root it returns.

Legacy code discards DOM by clearing it (`contentDiv.innerHTML = ""` on a view change), which
does not unmount a React root: its effects never clean up and its listeners stay registered.
Find where the piece's container is discarded and call `root.unmount()` there. For the menu bar
that is the loop at the top of `ChartManager._init` that resets `ds.menuBar`.

### 4. If other code appends to the container, give React a child of it

`addMenuIcon` must return an element synchronously and append it to `ds.menuBar`, and a React
root renders later. So the legacy code still creates `ds.menuBar`, and React renders into a
`<span style="display: contents">` inside it. The span has no box of its own, so React's
elements lay out as if they were direct children, and icons appended afterwards follow them as
they always did. This relies on `createMdvPortal` adding no wrapper elements, which is true
today (it adds context providers only).

### 5. Subscribe to legacy event sources with `useSyncExternalStore`

`DataStore` and `ChartManager` are not observable to React. Read them with
`useSyncExternalStore`: subscribe by adding a listener, return its removal, and read the current
value in the snapshot function. `DataStore.addListener(id, fn)` replaces any listener with the
same id, so generate a unique id per subscription. Reading the live value on mount, not a value
captured when the legacy code ran, also removes a class of stale-display bugs.

Reach `ChartManager` through `useChartManager()`. Pass the component only the thing it is about
(here, the datasource record).

### 6. Keep the markup and CSS classes

The component renders the same elements, class names and attributes as the legacy helper it
replaces (`MenuIcon` mirrors `createMenuIcon`, including the `role="tooltip"` that the microtip
tooltips need). This keeps the look identical, keeps React-rendered icons consistent with the
ones still added imperatively beside them, and makes "nothing changed" something that can be
checked by comparing the DOM before and after. Moving a piece to MUI components is a separate,
visible change.

### 7. Fix types where the legacy code is untyped, not with casts

TypeScript will object to things JS never declared. Fix them at the source:

- an untyped field gets a JSDoc type in the JS file
  (`/** @type {Record<string, ContextMenu>} */ this.layoutMenus`);
- a class only reachable through `BaseDialog.experiment[...]`, which is typed `unknown`, gets a
  named export and is imported directly;
- a legacy method that takes a DOM `Event` is given `event.nativeEvent`, not the React event;
- a new field on a shared record is added to its type in `src/charts/charts.d.ts`.

### 8. Test the component, then check the app

- A Vitest file under `src/tests/react/` renders the component with `@/react/hooks` mocked to
  return a small fake `ChartManager`, and the dialogs it opens mocked out. Cover what moved into
  the component: here the row count following `"filtered"`, the listener being removed on
  unmount, and the full screen state.
- In the running app, compare the container's children with the legacy build, use every
  control, and switch views back and forth while counting the listeners the component adds
  (`Object.keys(dataStore.listeners)`), to confirm that unmounting works.

## Considered and rejected

- **Converting `ChartManager.js` to TypeScript first.** A much larger change with no visible
  result, and it blocks every extraction behind it. Extractions add types at the edges they
  touch, which is where the types are needed.
- **A single React root for the whole app, with pieces rendered through portals.** This is the
  direction `createMdvPortal`'s own comment points to, but it needs the layout (`splitPane`, the
  chart containers) to be React-owned first. Separate roots work now and can be folded into one
  later without changing the components.
- **Restyling with MUI components in the same change**, which is how the main menu bar is
  built. It makes the change impossible to verify by comparison and leaves imperatively added
  icons looking different from their neighbours.
- **Replacing `addMenuIcon` with React state** (a registry of extra icons the component
  renders). Its callers keep the returned element and manipulate it directly, so this changes a
  plugin-facing API. It can follow once those callers are React.
- **Leaving both implementations behind a flag.** Two code paths to keep in step, for a change
  that is easy to revert as one commit.

## Consequences

- Each extraction removes a slice of `ChartManager.js` and replaces state spread across
  `ChartManager` with state owned by a component, which is testable on its own.
- There is one more React root per extracted piece: for the menu bar, one per datasource pane
  per view. This matches what charts and dialogs already do.
- A component updates on React's schedule, so a value driven by a legacy event appears a tick
  after the event, not within it. Code that reads the DOM straight after firing an event would
  notice; none in this repository does.
- Fields that only the extracted UI used disappear from `ChartManager` (`progressBars`,
  `isFullscreen`), and the container's direct children change (`ds.menuBar` now holds one
  `display: contents` span and any imperatively added icons). Nothing in this repository read
  either; code outside it that did would break.
- Behaviour that changed for the menu bar, both as fixes: the `fullscreenchange` listener is no
  longer added again on every view change, and a rebuilt bar shows the current filtered count
  where the legacy one showed the total until the next filter event.
- Not addressed: the menu icons are still `<span>`s with click handlers, so they are not
  reachable by keyboard, exactly as before.
