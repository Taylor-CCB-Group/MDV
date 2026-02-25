import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useGateManager } from "../gates/useGateManager";
import { useChartID, useConfig, useParamColumns } from "../hooks";
import type { LoadedDataColumn } from "@/charts/charts";
import { computeCentroid } from "../gates/gateUtils";
import type { DeckScatterConfig } from "../components/DeckScatterReactWrapper";
import { GeoJsonLayer, TextLayer } from "deck.gl";
import { getVivId } from "../components/avivatorish/MDVivViewer";
import { useSpatialLayers } from "../spatial_context";
import useGateActions from "./useGateActions";

const MAX_LABEL_LENGTH = 18;

function truncateGateLabel(name: string, maxLen: number = MAX_LABEL_LENGTH) {
    if (name.length <= maxLen) return name;
    return `${name.slice(0, maxLen - 1)}...`;
}
//todo: Be able to drag and drop a gate overlay layer
const useGateLayers = () => {
    const gateManager = useGateManager();
    const [cx, cy] = useParamColumns() as LoadedDataColumn<"double">[];
    const cz = useParamColumns()[2] as LoadedDataColumn<"double">;
    const config = useConfig<DeckScatterConfig>();
    const { dimension } = config;
    const is2d = dimension === "2d";
    const [draggingId, setDraggingId] = useState<string | null>(null);
    const [isHoveringLabel, setIsHoveringLabel] = useState(false);
    const draggingIdRef = useRef<string | null>(null);
    const [dragPos, setDragPos] = useState<[number, number] | null>(null);
    const gates = gateManager.gatesArray;
    const hasRebuiltGateColumnRef = useRef(false);
    const chartId = useChartID();
    const { selectionProps } = useSpatialLayers();
    const {
        onEditGate
    } = useGateActions();
    const {
        editingGateId,
        selectionFeatureCollection,
    } = selectionProps;

    const relevantGates = useMemo(() => {
        if (!cx || !cy) return [];
        return gateManager.gatesArray.filter(
            (gate) => gate.columns[0] === cx.field && gate.columns[1] === cy.field
        );
    }, [gateManager.gatesArray, cx, cy]);

    useEffect(() => {
        // Wait for columns to be loaded to update the gates column
        if (!cx || !cy || hasRebuiltGateColumnRef.current) return;
        hasRebuiltGateColumnRef.current = true;
        gateManager.rebuildGatesColumnWhenReady().catch(console.error);
    }, [cx, cy, gateManager]);

    const gateLabelLayer = useMemo(() => {
        if (!cx || !cy || gates.length === 0) return null;

        if (relevantGates.length === 0) return null;

        const layerData = relevantGates.map((gate) => {
            let position;

            if (gate.id === editingGateId) {
                // If gate is editing, we compute the centroid of the gate
                position = computeCentroid(selectionFeatureCollection);
            } else if (gate.id === draggingId && dragPos) {
                // If the label is getting dragged, we use the dragging position
                position = dragPos;
            } else {
                // Get the label position of the gate for normal gates
                position = gate.labelPosition;
            }
            // Get the average of z if it exists
            const z = is2d ? 0 : cz?.minMax ? (cz.minMax[0] + cz.minMax[1]) / 2 : 0;
            return {
                position: [position[0], position[1], z] as [number, number, number],
                text: truncateGateLabel(gate.name),
                gateId: gate.id,
            };
        });

        return new TextLayer({
            id: `text-layer-${getVivId(`${chartId}detail-react`)}`,
            data: layerData,
            getPosition: (d: { position: [number, number, number] }) => d.position,
            getText: (d: { text: string }) => d.text,
            getSize: 12,
            getColor: [255, 255, 255, 255],
            getBackgroundColor: [0, 0, 0, 100],
            getTextAnchor: "middle",
            getAlignmentBaseline: "center",
            padding: [3, 8],
            pickable: true,
            background: true,
            backgroundPadding: [4, 2],
            backgroundBorderRadius: 4,
            updateTriggers: {
                getPosition: [gates.length, dragPos, selectionFeatureCollection],
                getColor: [draggingId],
                getBackgroundColor: [draggingId],
            },
            onHover(pickingInfo) {
                if (pickingInfo.index !== -1) {
                    setIsHoveringLabel(true);
                } else {
                    setIsHoveringLabel(false);
                }
            },
            onDrag(pickingInfo) {
                const currentId = draggingIdRef.current;
                if (!currentId || !dragPos || !pickingInfo.object || !pickingInfo.coordinate) return;

                const currentPosition: [number, number] = [pickingInfo.coordinate[0], pickingInfo.coordinate[1]];

                // Update the drag position of the label with the current dragged position
                setDragPos(currentPosition);
            },
            onDragStart(pickingInfo) {
                // Not allowing dragging when editing a gate to avoid issues
                if (!pickingInfo.object || !pickingInfo.coordinate || editingGateId) return false;

                // Set teh dragging id and drag position state
                draggingIdRef.current = pickingInfo.object.gateId;
                setDraggingId(pickingInfo.object.gateId);
                setDragPos([pickingInfo.coordinate[0], pickingInfo.coordinate[1]]);
                return true;
            },
            onDragEnd() {
                const currentId = draggingIdRef.current;
                if (!currentId) return;

                const position = dragPos;
                if (position) {
                    // Update the label position of the gate, when the dragging ends
                    gateManager.updateGate(currentId, { labelPosition: position });
                }
                draggingIdRef.current = null;
                setDraggingId(null);
                setDragPos(null);
            },
        });
    }, [
        cx,
        cy,
        cz,
        is2d,
        gates,
        draggingId,
        dragPos,
        gateManager,
        chartId,
        editingGateId,
        relevantGates,
        selectionFeatureCollection,
    ]);

    const gateDisplayLayer = useMemo(() => {
        if (!cx || !cy || gateManager.gatesArray.length === 0) return null;

        // Filter the editing gate
        const filteredGates = relevantGates.filter((gate) => gate.id !== editingGateId);

        // Create the objects of features to supply as the feature collection
        // Store gate id and gate name as feature properties
        const features = filteredGates.flatMap((gate) =>
            gate.geometry.features.map((feature) => ({
                ...feature,
                properties: {
                    ...feature.properties,
                    gateId: gate.id,
                    gateName: gate.name,
                },
            })),
        );

        return new GeoJsonLayer({
            id: `gate_${getVivId(`${chartId}detail-react`)}`,
            data: {
                type: "FeatureCollection",
                features,
            } as any,
            filled: true,
            getFillColor: [76, 175, 80, 30], // Semi-transparent green fill
            getLineColor: [76, 175, 80, 200], // Green border
            getLineWidth: 2,
            lineWidthMinPixels: 4,
            pickable: true,
            onClick(pickingInfo) {
                const gateId = pickingInfo.object?.properties?.gateId;
                onEditGate(gateId);
            },
        });
    }, [
        relevantGates,
        chartId,
        editingGateId,
        gateManager,
        onEditGate,
        cx,
        cy,
    ]);

    const getCursor = useCallback(
        ({ isDragging }: { isDragging: boolean; isHovering: boolean }) => {
            if (draggingId || isDragging) return "grabbing";

            return "grab";
        },
        [draggingId],
    );

    const dragPan = !(draggingId || isHoveringLabel);

    return {
        gateLabelLayer,
        gateDisplayLayer,
        controllerOptions: { dragPan },
        getCursor,
    };
};

export default useGateLayers;
