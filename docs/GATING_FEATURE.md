# Gating Feature

Documentation for the gating feature on 2D scatter plots (Viv and Deck).

## Overview

Gating lets users:

- Draw named regions (gates) on scatter plots using rectangle, polygon, or freehand tools
- Save a selection as a gate with a name and color
- Edit, rename, delete gates, or export gate geometry as GeoJSON
- Have per-cell gate membership stored in a `__gates__` column for filtering

Each gate is tied to a specific X/Y column pair. Only gates for the current chart axes are shown.

## File Structure


| File                                                   | Purpose                                                                            |
| ------------------------------------------------------ | ---------------------------------------------------------------------------------- |
| `src/react/gates/GateManager.ts`                       | Gate lifecycle (add/update/delete), `__gates__` column   |
| `src/react/gates/types.ts`                             | `Gate` interface               |
| `src/react/gates/useGateManager.ts`                    | Hook: one GateManager per DataStore (WeakMap)                                      |
| `src/react/gates/gateUtils.ts`                         | Geometry helpers |
| `src/react/hooks/useGateLayers.ts`                     | Deck layers: gate polygons, draggable labels              |
| `src/react/hooks/useGateActions.ts`                    | Actions: save, delete, rename, edit, etc            |
| `src/react/components/SelectionOverlay.tsx`            | Toolbar: draw tools, Save as Gate, Manage Gates, Confirm/Cancel when editing       |
| `src/react/components/ManageGateDialog.tsx`            | Gate list and all related operations                           |
| `src/react/spatial_context.tsx`                        | Selection state, editable layer       |
| `VivScatterComponent.tsx` / `DeckScatterComponent.tsx` | Add gate layers from useGateLayers, tooltips, dragPan control                      |


## How It Works

**Data:** Gates live in `GateManager.gates` (Map) and in `dataStore.config.gates` for persistence. A multitext column `__gates__` stores which gate name(s) each cell belongs to (value `"N/A"` when none). Membership is computed with point-in-polygon on the gate’s X/Y columns.

**Create:** User draws a shape → selection is in `selectionFeatureCollection`. “Save selection as Gate” → name (and color) → GateManager adds the gate, updates the `__gates__` column for cells inside the polygon, then persists to config.

**Edit:** User clicks a gate label or uses Manage Gates → Edit. That gate’s geometry is loaded into the selection layer and `editingGateId` is set. User reshapes with Pan/Modify, then Confirm updates the gate and clears the selection, or Cancel discards.

**Display:** `useGateLayers` builds a polygon layer (all gates for current axes except the one being edited) and a text layer (draggable gate label). The scatter views add these layers to Deck and pass through `controllerOptions.dragPan` so pan is disabled while dragging a label.

## Logic Summary

- **GateManager** owns the `__gates__` column (create if missing, load from config on init). On add/update/delete it runs point-in-polygon for affected gates and writes gate names into the column; it also keeps `dataStore.config.gates` and dirty flags in sync.
- **Relevant gates** are those whose `columns` match the chart’s current X and Y; only those are rendered and listed in Manage Gates associated with that particular chart.
- **Editing:** When `editingGateId` is set, the selection layer shows that gate’s geometry (editable) and uses its color; the gate display layer hides that gate so it isn’t drawn twice.

## GateManager Details

- **When constructed** (`new GateManager(dataStore)`):
  - Looks for an existing multitext column `__gates__`; if missing, creates it with default value `"N/A"` for every row.
  - Loads any saved gates from `dataStore.config.gates` / `dataStore.gates` into an in-memory `Map<string, Gate>`.
  - Rebuilds the `__gates__` column when the underlying X/Y data columns are loaded.
- **Gate membership**:
  - For each gate, only the configured X/Y columns are used for point-in-polygon checks.
  - For each row, if `(x, y)` is inside the gate polygon, the gate’s name is added to that row’s `__gates__` multitext cell (supporting multiple gates per row).
  - If a gate is renamed or deleted, GateManager rebuilds the `values` array of the multitext column so the dropdowns and filters only show current gate names.
- **Persistence**:
  - After any change, GateManager writes the gates array to `dataStore.config.gates` and `dataStore.gates`, and marks both metadata and the `__gates__` column as dirty so it is saved.

## Gate Layers and Labels

- **Gate display layer** (`useGateLayers`):
  - Renders polygons for all gates that match the current X/Y axes, except the one currently being edited.
  - Each polygon uses its gate color (with a light fill and stronger outline).
  - Tooltips show the gate name and a hint that labels are clickable for editing.
- **Gate label layer** (`useGateLayers`):
  - Renders one text label per gate (truncated name), positioned at `labelPosition` or the centroid of the gate geometry.
  - Labels are **clickable** (to start editing that gate) and **draggable**; on drag end, GateManager updates `labelPosition`.
  - While a label is being dragged or hovered, `dragPan` is disabled so the chart does not pan instead of moving the label.

## Editing Gates

- **Starting an edit**:
  - From the toolbar: open **Manage Gates**, choose **Edit Geometry**.
  - Or directly: click the gate’s label in the scatter plot.
  - In both cases, the gate’s geometry is copied into `selectionFeatureCollection` and `editingGateId` is set.
- **During edit**:
  - The selection layer (editable GeoJSON) shows the gate shape with its color and edit handles.
  - Only Pan and Modify tools are shown while editing to keep interactions simple.
  - Filters are driven from the current selection polygon, so changing the shape immediately affects which points are inside.
- **Finishing or cancelling**:
  - **Confirm**: writes the edited geometry and a fresh centroid-based `labelPosition` back via GateManager, recomputes membership in `__gates__`, then clears `editingGateId` and selection.
  - **Cancel**: clears `editingGateId` and selection without changing the stored gate.


## Future Work:

- Make the gates column as a derived column
- The coloring of gates need to work with the color by of the chart, but the multitext color by needs to be implemented properly.
- Combine density fields and gates as suggested by a user
- UI/UX enhancements based on user feedback