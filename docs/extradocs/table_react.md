# Table Chart React

Documentation for the React-based Table Chart component using SlickGrid-React.

## Overview

The React Table Chart provides:
- Large dataset rendering via virtual scrolling
- Column sorting and reordering (persisted to config)
- Cell editing for editable columns
- Find & Replace functionality
- Cross-chart highlighting/selection

## File Structure

| File | Purpose |
|------|---------|
| `TableChartReactWrapper.tsx` | Chart class, config adaptation, getConfig() serialization |
| `TableChartReactComponent.tsx` | Main React component, composes hooks |
| `useSlickGridReact.ts` | Grid lifecycle, event handlers, core state |
| `useSortedFilteredIndices.ts` | Sorting logic based on config.sort |
| `useFindReplace.ts` | Find/Replace dialog logic |
| `useEditCell.ts` | Cell editing validation and persistence |
| `SlickGridDataProvider.ts` | Adapter between MDV data and SlickGrid API |
| `valueReplacementUtil.ts` | Read/write values across MDV data types |

## Data Flow

DataStore (raw data) ->  useOrderedParamColumns (column definitions from config.param + config.order) -> useSortedFilteredIndices (applies config.sort to filteredIndices) -> SlickGridDataProvider (adapts MDV data model to SlickGrid DataView API) -> SlickGrid (renders visible rows via virtual scrolling)

## State Management

**Config (MobX Observable)** - Persisted state:
- `config.sort` - Sort state: `{ columnId, ascending }` (syncs config.sort and internal grid state when sort changes)
- `config.order` - Column order: `Record<field, position>`
- `config.column_widths` - Initial widths (read once on mount)
- `config.include_index` - Show index column

**React State** - Volatile UI state:
- Find/Replace dialog open state
- Search results and navigation

**SlickGrid** - Runtime state:
- Actual column widths during session
- Cell selection, scroll position

## Key Patterns

### 1. Refs vs Values

- `useEditCell` and `useSlickGridReact` uses refs (`orderedParamColumnsRef`, `sortedFilteredIndicesRef`) because SlickGrid stores callback references and we need to avoid stale closures.

- `useFindReplace` uses direct values because the dialog must close before sort/filter changes, so callbacks recreate with fresh data.

### 2. Grid Row vs Data Index

- **Grid Row**: Position in displayed grid (0 to length-1)
- **Data Index**: Position in underlying data arrays

The `indices` array in SlickGridDataProvider maps grid rows to data indices.

### 3. External Sorting

We manage sorting externally via `useSortedFilteredIndices` instead of SlickGrid's built-in sorting because:
- Sort state persists via `config.sort` (MobX)
- MDV data is in typed arrays, not SlickGrid format
- `SlickGridDataProvider.sort()` is a no-op

### 4. Sort Synchronization

Bidirectional sync between `config.sort` and SlickGrid's visual sort indicators:
- `suppressSortSyncRef` prevents feedback loops
- Grid onSort → updates config.sort
- config.sort changes → updates grid visual state

### 5. Selection/Highlighting

- `selectionSourceRef` tracks whether selection is 'user' or 'programmatic'
- Prevents echo effects when cross-chart highlighting triggers grid selection

### 6. Editing Data

- `setCellValueFromString()` modifies typed array
- `dataStore.dataChanged([column])` notifies listeners
- `grid.invalidate()` + `grid.render()` refreshes display
- Note: dataProvider is NOT recreated since array references don't change.

## Supported Data Types

| Type | Storage | Notes |
|------|---------|-------|
| text | Uint8Array index to values[] | Max 256 values |
| text16 | Uint16Array index to values[] | Max 65536 values |
| double | Float64Array | IEEE 754 |
| int32/integer | Int32Array | 32-bit signed |
| multitext | Uint16Array[len] to values[] | Multiple per cell |
| unique | Uint8Array[len] UTF-8 | Fixed-length string |


