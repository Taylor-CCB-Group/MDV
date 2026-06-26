# React Table Performance Notes

This note records the React table performance investigation from the 6.2M-row / 100-column stress test. It is intended as context for a follow-up PR, not as a finished performance plan.

## Baseline

The comparison was between:

- legacy table: `src/charts/TableChart.js`
- React table: `src/react/components/TableChartReactComponent.tsx` and `src/react/hooks/useSlickGridReact.ts`

The committed improvement before this note was the stable `SlickGridDataProvider` change:

- keep one provider instance in `useSlickGridReact`
- update it in place with new columns / indices / `include_index`
- call `grid.setData(...)` only on grid creation
- call `grid.invalidate()` after provider updates

That change is still considered valid because it removes an unnecessary provider recreation and grid `setData` cycle.

## What Was Measured

Temporary logs were added and later removed for:

- legacy table filter update
- legacy table sort
- React table total filter update
- React table total sort update
- React table indices copy/sort time
- React table provider update and `grid.invalidate()` time

Typical no-sort filter logs looked like this:

```text
[Legacy Table] filter update: 6.9ms { rows: 1033333, sorted: false }
methodfilterCategories: 70.9ms
[React Table] indices copy: 0.7ms { rows: 1033333, sorted: false }
[React Table] grid update: 7.8ms { rows: 1033333, providerUpdateMs: 0, invalidateMs: 7.8 }
[React Table] filter update: 87.7ms { rows: 1033333, sorted: false }
```

The important observation is that the React table's own measured work was small:

- index copy: around 1ms
- provider update: approximately 0ms
- grid invalidate: usually a few ms

But the total React filter update remained around 50-100ms or more.

## Conclusions

### No-sort filtering

The React table latency during selection-dialog filtering is not primarily caused by:

- `SlickGridDataProvider` cell extraction
- the `Uint32Array` copy in `useSortedFilteredIndices`
- provider update cost
- `grid.invalidate()` cost
- sorting, when `sorted: false`

The gap appears to be before the table update effect runs:

```text
DataStore filtered event
-> React scheduling / effects
-> selection dialog updates
-> filtered indices arrive
-> React table effect runs
-> provider update + grid invalidate
```

The legacy table is faster in this benchmark because it updates imperatively during its datastore listener path:

```text
TableChart.onDataFiltered()
-> DataModel.updateModel()
-> optional sort
-> grid update
```

So the broad React table filter timing includes React scheduling and surrounding React component work, not just table work.

### Sorting

Both tables currently block the main thread while sorting table rows.

The legacy table does not use `SortableDimension` / `sortWorker.js` for table header sorting. Its path is:

```text
TableChart._sort()
-> DataModel.sort()
-> typed array sort on main thread
-> grid.invalidateAllRows()
-> grid.render()
```

The React table sort path is also synchronous on the main thread inside `useSortedFilteredIndices`.

Therefore, a sort worker is still the most meaningful change for long sort operations and for UI responsiveness during sorting.

## Experiments Tried And Reverted

### Lazy cell extraction

Idea:

- make `SlickGridDataProvider.getItem()` avoid reading every column for each row
- resolve values lazily when SlickGrid asks for a specific cell

Result:

- did not improve the no-sort filter benchmark
- added behavior complexity around row object shape and edit paths

Status: reverted / stashed, not recommended for this PR.

### Remove no-sort `Uint32Array` copy

Idea:

- when no sort is active, return filtered indices directly instead of copying them

Result:

- did not materially improve the no-sort filter benchmark
- later split logs showed the copy itself was around 1ms

Status: reverted / stashed, not recommended as a standalone fix.

### Table-specific synchronous filtered indices

Idea:

- bypass `DataStore.getFilteredIndices()` worker for the React table
- build indices directly from `dataStore.filterArray`, similar to `DataModel.updateModel()`

Result:

- did not materially improve the no-sort filter benchmark
- later split logs showed the remaining delay was mostly before the table update effect ran
- this also weakens the shared React chart data architecture

Status: reverted, not recommended.

## Recommended Next Steps

1. Keep the stable data-provider commit.

2. Do not continue with micro-optimizations for no-sort filtering unless a profiler shows a new local bottleneck.

3. If filter latency is revisited, instrument the shared React filtered-indices pipeline and the selection dialog, not the table provider:

```text
DataStore filtered event timestamp
getFilteredIndices request timestamp
getFilteredIndices resolve timestamp
useSimplerFilteredIndices setState timestamp
useSortedFilteredIndices effect timestamp
useSlickGridReact provider/grid update timestamp
```

4. For a high-impact table performance change, prioritize a React table sort worker:

- move the expensive `Uint32Array.sort(...)` work off the main thread
- preserve current sort behavior for numeric, text, unique, multitext, and `__index__`
- ignore stale worker results when filters/sorts change quickly
- keep the current table visible while the worker sorts

5. If parity with legacy no-sort filtering is mandatory, the likely architectural option is an imperative table update path. That would trade React chart consistency for speed and should be discussed explicitly before implementation.
