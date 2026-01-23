import { useEffect, useId, useMemo, useState } from "react";
import { useGateManager } from "../gates/useGateManager";
import { useConfig, useParamColumns } from "../hooks";
import type { LoadedDataColumn } from "@/charts/charts";
import { computeCentroid } from "../gates/gateUtils";
import type { DeckScatterConfig } from "../components/DeckScatterReactWrapper";
import { TextLayer } from "deck.gl";
import { MonkeyPatchEditableGeoJsonLayer } from "@/lib/deckMonkeypatch";

//todo: Be able to drag and drop a gate overlay layer and fix issues related to the whole thing moving when the gate layer or text layer is moved
const useGateLayers = () => {
    const gateManager = useGateManager();
    const id = useId();
    const [cx, cy] = useParamColumns() as LoadedDataColumn<"double">[];
    const cz = useParamColumns()[2] as LoadedDataColumn<"double">;
    const config = useConfig<DeckScatterConfig>();
    const { dimension } = config;
    const is2d = dimension === "2d";
    const [labelPositions, setLabelPositions] = useState<Map<string, [number, number]>>(new Map());
    const [draggingId, setDraggingId] = useState<string | null>(null);
    const [dragPos, setDragPos] = useState<[number, number] | null>(null);
    const gates = gateManager.gatesArray;

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
        if (cx && cy) {
            gateManager.rebuildGatesColumnWhenReady().catch(console.error);
        }
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
                text: gate.name,
                gateId: gate.id,
            };
        });

        return new TextLayer({
            id: `text-layer-${id}`,
            data: layerData,
            getPosition: (d: { position: [number, number, number] }) => d.position,
            getText: (d: { text: string }) => d.text,
            getSize: 12,
            getColor: [255, 255, 255, 255],
            getTextAnchor: "middle",
            getAlignmentBaseline: "center",
            padding: [3, 8],
            billboard: true,
            pickable: true,
            updateTriggers: {
                getPosition: [gates.length, labelPositions.size],
                getColor: [draggingId],
                getBackgroundColor: [draggingId],
            },
            onDrag(pickingInfo) {
                if (!draggingId || !dragPos || !pickingInfo.object || !pickingInfo.coordinate) return;
                console.log("onDrag.....");

                const currentPosition: [number, number] = [pickingInfo.coordinate[0], pickingInfo.coordinate[1]];

                setLabelPositions((prev) => {
                    const newPositions = new Map(prev);
                    newPositions.set(draggingId, currentPosition);
                    return newPositions;
                });

                setDragPos(currentPosition);
            },
            onDragStart(pickingInfo) {
                if (!pickingInfo.object || !pickingInfo.coordinate) return false;
                console.log("onDragStart.....");

                setDraggingId(pickingInfo.object.gateId);
                setDragPos([pickingInfo.coordinate[0], pickingInfo.coordinate[1]]);
                console.log("drag position", [pickingInfo.coordinate[0], pickingInfo.coordinate[1]]);
                return true;
            },
            onDragEnd(_pickingInfo, _event) {
                if (!draggingId) return;
                console.log("onDragEnd.....");

                const position = labelPositions?.get(draggingId);
                if (position) {
                    gateManager.updateGate(draggingId, { labelPosition: position });
                }
                setDraggingId(null);
                setDragPos(null);
            },
        });
    }, [cx, cy, cz, is2d, id, gates, draggingId, dragPos, gateManager, labelPositions]);

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
        }
        );

        return new MonkeyPatchEditableGeoJsonLayer({
            id: `gate-layer-${id}`,
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
            onDrag(pickingInfo, event) {
                console.log("onDrag.....")
                console.log("pickingInfo: ", pickingInfo);
            },
            onDragStart(pickingInfo, event) {
                console.log("onDragStart.....")
                console.log("pickingInfo: ", pickingInfo);
            },
            onDragEnd(pickingInfo, event) {
                console.log("onDragEnd.....")
                console.log("pickingInfo: ", pickingInfo);
            },
            
        });
    }, [cx, cy, id, gates]);

    return {
        gateLabelLayer,
        gateOverlayLayer,
        draggingId,
    };
};

export default useGateLayers;
