// deck.gl-community
// SPDX-License-Identifier: MIT
// Copyright (c) vis.gl contributors / Peter Todd

import {
    type DraggingEvent,
    DrawPolygonByDraggingMode,
    type Feature,
    type ModeProps,
    type Polygon,
    type Position,
    type SimpleFeatureCollection,
    type StopDraggingEvent,
} from "@deck.gl-community/editable-layers";

function rectangleFeature(p1: Position, p2: Position, properties: Record<string, unknown>): Feature<Polygon> {
    const minX = Math.min(p1[0], p2[0]);
    const minY = Math.min(p1[1], p2[1]);
    const maxX = Math.max(p1[0], p2[0]);
    const maxY = Math.max(p1[1], p2[1]);
    return {
        type: "Feature",
        properties,
        geometry: {
            type: "Polygon",
            coordinates: [
                [
                    [minX, minY],
                    [maxX, minY],
                    [maxX, maxY],
                    [minX, maxY],
                    [minX, minY],
                ],
            ],
        },
    };
}

function withoutTentativeRectangle(data: SimpleFeatureCollection): SimpleFeatureCollection["features"] {
    return data.features.filter((feature) => !feature.properties?.tentativeRectangle);
}

/**
 * Draw a rectangle by dragging
 *
 * !this is already a feature, we just needed to pass in `dragToDraw` in props...
 * leaving this for now, for review later, as the other version isn't an immediate drop-in replacement.
 */
export class DrawRectangleByDraggingMode extends DrawPolygonByDraggingMode {
    /**
     * Built from this event rather than from the tentative rectangle in `props.data`: the layer is
     * only given that once React re-renders, so a quick drag could otherwise finish with the
     * previous selection.
     */
    handleStopDragging(event: StopDraggingEvent, props: ModeProps<SimpleFeatureCollection>) {
        const rectangle = rectangleFeature(event.pointerDownMapCoords, event.mapCoords, {});
        props.onEdit({
            updatedData: { ...props.data, features: [...withoutTentativeRectangle(props.data), rectangle] },
            editType: "finishDrawingRectangle",
            editContext: {
                position: event.mapCoords,
            },
        });
    }

    handleDraggingAux(event: DraggingEvent, props: ModeProps<SimpleFeatureCollection>) {
        const feature = rectangleFeature(event.pointerDownMapCoords, event.mapCoords, { tentativeRectangle: true });
        props.onEdit({
            updatedData: { ...props.data, features: [...withoutTentativeRectangle(props.data), feature] },
            editType: "setTentativeRectangle",
            editContext: {
                position: event.mapCoords,
            },
        });
    }
}
