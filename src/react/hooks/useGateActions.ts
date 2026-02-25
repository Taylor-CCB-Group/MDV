import { useCallback } from "react";
import { useGateManager } from "../gates/useGateManager";
import { useParamColumns } from "../hooks";
import type { Gate } from "../gates/types";
import { computeCentroid, generateGateId } from "../gates/gateUtils";
import { useSpatialLayers } from "../spatial_context";
import { action } from "mobx";
import { getEmptyFeatureCollection } from "../deck_state";
import { useChart } from "../context";
import type { DeckScatterConfig } from "../components/DeckScatterReactWrapper";
import { Tools } from "../components/SelectionOverlay";

const useGateActions = () => {
    const gateManager = useGateManager();
    const paramColumns = useParamColumns();
    const { selectionProps } = useSpatialLayers();
    const {
        selectionFeatureCollection,
        setEditingGateId,
        editingGateId,
        setSelectedTool,
        setSelectionFeatureCollection,
        setSelectionMode,
    } = selectionProps;
    const chart = useChart<DeckScatterConfig>();

    const clearSelection = useCallback(() => {
        action(() => {
            chart.config.selectionFeatureCollection = getEmptyFeatureCollection();
        })();
    }, [chart.config]);

    const onSaveGate = useCallback(
        (gateName: string) => {
            if (gateManager.hasGateName(gateName)) {
                throw new Error("A gate with this name already exists");
            }

            // Get X and Y columns
            const [xCol, yCol] = paramColumns.slice(0, 2);
            if (!xCol || !yCol) {
                throw new Error("Chart must have X and Y axes");
            }

            const centroid = computeCentroid(selectionFeatureCollection);

            // Create gate
            const gate: Gate = {
                id: generateGateId(),
                name: gateName.trim(),
                geometry: selectionFeatureCollection,
                columns: [xCol.field, yCol.field],
                createdAt: Date.now(),
                labelPosition: centroid,
            };

            // Add to gate store
            gateManager.addGate(gate);

            // Clear selection
            clearSelection();
        },
        [gateManager, paramColumns, selectionFeatureCollection, clearSelection],
    );

    const onDeleteGate = useCallback(
        (gateId: string) => {
            gateManager.deleteGate(gateId);
        },
        [gateManager],
    );

    const onRenameGate = useCallback(
        (gateId: string, newName: string) => {
            if (gateManager.hasGateName(newName)) {
                throw new Error("A gate with this name already exists");
            }
            gateManager.updateGate(gateId, { name: newName });
        },
        [gateManager],
    );

    /**
     * Export the geometry of the gate in GeoJSON format
     */
    const onExportClick = useCallback(
        (gateId: string) => {
            const gate = gateManager.gatesArray.find((g) => g.id === gateId);
            if (!gate) return;

            const geojson = gate.geometry;
            const jsonStr = JSON.stringify(geojson, null, 2);
            const blob = new Blob([jsonStr], { type: "application/geo+json" });
            const url = URL.createObjectURL(blob);
            const link = document.createElement("a");
            link.href = url;
            link.download = `${gate.name || "gate"}.geojson`;
            link.click();
            URL.revokeObjectURL(url);
        },
        [gateManager],
    );

    /**
     * Set editing gate id and change selection mode to pan mode
     */
    const onEditGate = useCallback(
        (gateId: string) => {
            if (gateId) {
                const gate = gateManager.gates.get(gateId);
                if (gate) {
                    setSelectionFeatureCollection(JSON.parse(JSON.stringify(gate.geometry)));
                    setEditingGateId(gateId);
                    const panMode = Tools["pan"].mode;
                    setSelectionMode(new panMode());
                    setSelectedTool("Pan");
                }
            }
        },
        [setSelectionFeatureCollection, setEditingGateId, setSelectionMode, setSelectedTool, gateManager],
    );

    /**
     * Update the gate with the new geometry and label position of the gate
     * Compute centroid of label position of editing gate
     */
    const onConfirmEditGate = useCallback(() => {
        if (!editingGateId) return;
        const currentGeometry = chart.config.selectionFeatureCollection;

        const newLabelPosition = computeCentroid(currentGeometry);
        gateManager.updateGate(editingGateId, { geometry: currentGeometry, labelPosition: newLabelPosition });
        setEditingGateId(null);
        clearSelection();
    }, [gateManager, clearSelection, chart.config, setEditingGateId, editingGateId]);

    /**
     * Clear the editing id and clear selection
     */
    const onCancelEditGate = useCallback(() => {
        setEditingGateId(null);
        clearSelection();
    }, [setEditingGateId, clearSelection]);

    return {
        onSaveGate,
        onDeleteGate,
        onRenameGate,
        onExportClick,
        onEditGate,
        onConfirmEditGate,
        onCancelEditGate,
    };
};

export default useGateActions;
