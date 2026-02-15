import { useCallback } from "react";
import { useGateManager } from "../gates/useGateManager";
import { useParamColumns } from "../hooks";
import type { Gate } from "../gates/types";
import { generateGateId } from "../gates/gateUtils";
import { useSpatialLayers } from "../spatial_context";
import { action } from "mobx";
import { getEmptyFeatureCollection } from "../deck_state";
import { useChart } from "../context";
import type { DeckScatterConfig } from "../components/DeckScatterReactWrapper";

const useGateActions = () => {
    const gateManager = useGateManager();
    const paramColumns = useParamColumns();
    const { selectionProps } = useSpatialLayers();
    const { selectionFeatureCollection } = selectionProps;
    const chart = useChart<DeckScatterConfig>();

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

            // Create gate
            const gate: Gate = {
                id: generateGateId(),
                name: gateName.trim(),
                geometry: selectionFeatureCollection,
                columns: [xCol.field, yCol.field],
                createdAt: Date.now(),
            };

            // Add to gate store
            gateManager.addGate(gate);

            // Clear selection
            action(() => {
                chart.config.selectionFeatureCollection = getEmptyFeatureCollection();
            })();
        },
        [gateManager, paramColumns, selectionFeatureCollection, chart.config],
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

    return {
        onSaveGate,
        onDeleteGate,
        onRenameGate,
        onExportClick,
    };
};

export default useGateActions;
