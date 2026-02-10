import { useEffect, useMemo, useRef, useState } from "react";
import { useGateManager } from "../gates/useGateManager";
import { useChartID, useConfig, useParamColumns } from "../hooks";
import type { LoadedDataColumn } from "@/charts/charts";
import { computeCentroid } from "../gates/gateUtils";
import type { DeckScatterConfig } from "../components/DeckScatterReactWrapper";
import { TextLayer } from "deck.gl";
import { MonkeyPatchEditableGeoJsonLayer } from "@/lib/deckMonkeypatch";
import { getVivId } from "../components/avivatorish/MDVivViewer";

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
    const [labelPositions, setLabelPositions] = useState<Map<string, [number, number]>>(new Map());
    const [draggingId, setDraggingId] = useState<string | null>(null);
    const [isHoveringLabel, setIsHoveringLabel] = useState(false);
    const draggingIdRef = useRef<string | null>(null);
    const [dragPos, setDragPos] = useState<[number, number] | null>(null);
    const gates = gateManager.gatesArray;
    const hasRebuiltGateColumnRef = useRef(false);
    const chartId = useChartID();
    const vivLayerId = getVivId(`${chartId}detail-react`);

    useEffect(() => {
        const positions = new Map<string, [number, number]>();
        gates.forEach((gate) => {
            if (gate.labelPosition) {
                positions.set(gate.id, gate.labelPosition);
            } else {
                const centroid = computeCentroid(gate.geometry);
                positions.set(gate.id, [centroid[0], centroid[1]]);
            }
        });
        setLabelPositions(positions);
    }, [gates]);

    useEffect(() => {
        // Wait for columns to be loaded to update the gates column
        if (!cx || !cy || hasRebuiltGateColumnRef.current) return;
        hasRebuiltGateColumnRef.current = true;
        gateManager.rebuildGatesColumnWhenReady().catch(console.error);
    }, [cx, cy, gateManager]);

    const gateLabelLayer = useMemo(() => {
        if (!cx || !cy || gates.length === 0) return null;

        const relevantGates = gates.filter((gate) => gate.columns[0] === cx.field && gate.columns[1] === cy.field);

        if (relevantGates.length === 0) return null;

        const layerData = relevantGates.map((gate) => {
            let position = labelPositions?.get(gate.id);
            if (!position) {
                const centroid = computeCentroid(gate.geometry);
                position = [centroid[0], centroid[1]];
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
            id: `text-layer-${vivLayerId}`,
            data: layerData,
            getPosition: (d: { position: [number, number, number] }) => d.position,
            getText: (d: { text: string }) => d.text,
            getSize: 12,
            // todo: fix the color based on the theme
            getColor: [255, 255, 255, 255],
            getTextAnchor: "middle",
            getAlignmentBaseline: "center",
            padding: [3, 8],
            // billboard: true,
            pickable: true,
            updateTriggers: {
                getPosition: [gates.length, labelPositions.size],
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

                setLabelPositions((prev) => {
                    const newPositions = new Map(prev);
                    newPositions.set(currentId, currentPosition);
                    return newPositions;
                });

                setDragPos(currentPosition);
            },
            onDragStart(pickingInfo) {
                if (!pickingInfo.object || !pickingInfo.coordinate) return false;

                draggingIdRef.current = pickingInfo.object.gateId;
                setDraggingId(pickingInfo.object.gateId);
                setDragPos([pickingInfo.coordinate[0], pickingInfo.coordinate[1]]);
                return true;
            },
            onDragEnd() {
                const currentId = draggingIdRef.current;
                if (!currentId) return;

                const position = labelPositions?.get(currentId);
                if (position) {
                    gateManager.updateGate(currentId, { labelPosition: position });
                }
                draggingIdRef.current = null;
                setDraggingId(null);
                setDragPos(null);
            },
        });
    }, [cx, cy, cz, is2d, gates, draggingId, dragPos, gateManager, labelPositions, vivLayerId]);

    const gateOverlayLayer = useMemo(() => {
        if (!cx || !cy || gates.length === 0) return null;

        const relevantGates = gates.filter((gate) => gate.columns[0] === cx.field && gate.columns[1] === cy.field);

        if (relevantGates.length === 0) return null;

        //* Why do we need to do this?
        const layerData = relevantGates.flatMap((gate) => {
            return gate.geometry.features.map((feature) => {
                return {
                    ...feature,
                    properties: {
                        ...feature.properties,
                        gateId: gate.id,
                        gateName: gate.name,
                    },
                };
            });
        });

        return new MonkeyPatchEditableGeoJsonLayer({
            id: `gate-layer-${vivLayerId}`,
            // todo: Fix typing here
            data: {
                type: "FeatureCollection",
                features: layerData,
            } as any,
            filled: true,
            getFillColor: [76, 175, 80, 30], // Semi-transparent green fill
            getLineColor: [76, 175, 80, 200], // Green border
            getLineWidth: 2,
            lineWidthMinPixels: 1,
            pickable: false,
        });
    }, [cx, cy, gates, vivLayerId]);

    return {
        gateLabelLayer,
        gateOverlayLayer,
        draggingId,
        isHoveringLabel,
    };
};

export default useGateLayers;
