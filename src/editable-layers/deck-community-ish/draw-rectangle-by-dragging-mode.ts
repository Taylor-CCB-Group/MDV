// deck.gl-community
// SPDX-License-Identifier: MIT
// Copyright (c) vis.gl contributors / Peter Todd

import { DrawPolygonByDraggingMode, type DraggingEvent, type FeatureCollection, type ModeProps, type StopDraggingEvent } from "@deck.gl-community/editable-layers";
import type { Feature } from "@turf/helpers";

/**
 * Draw a rectangle by dragging
 * 
 * This is a feature that should be in deck.gl-community/editable-layers I think...
 * Will make a PR for it soon.
 */
export class DrawRectangleByDraggingMode extends DrawPolygonByDraggingMode {
  handleStopDragging(event: StopDraggingEvent, props: ModeProps<FeatureCollection>) {
    // I don't think we want the clickSequence stuff from DrawPolygonByDraggingMode
    // we should be safe to assume here that the last feature is the one we're working on
    props.data.features[props.data.features.length - 1].properties = {};
    props.onEdit({
      updatedData: props.data,
      editType: 'finishDrawingRectangle',
      editContext: {
        position: event.mapCoords
      }
    });
  }

  handleDraggingAux(event: DraggingEvent, props: ModeProps<FeatureCollection>) {
    const p1 = event.pointerDownMapCoords;
    const p2 = event.mapCoords;
    const minX = Math.min(p1[0], p2[0]);
    const minY = Math.min(p1[1], p2[1]);
    const maxX = Math.max(p1[0], p2[0]);
    const maxY = Math.max(p1[1], p2[1]);
    const start = [minX, minY];
    const end = [maxX, maxY];
    const feature: Feature = {
      type: "Feature",
      properties: {
        tentativeRectangle: true
      },
      geometry: {
        type: "Polygon",
        coordinates: [
          [start, [end[0], start[1]], end, [start[0], end[1]], start]
        ]
      }
    };
    const features = [...props.data.features];
    if (features.length > 0 && features[features.length - 1].properties?.tentativeRectangle) {
      features[features.length - 1] = feature as any; //todo improve typing
    } else {
      features.push(feature as any);
    }
    const updatedData: FeatureCollection = {
      ...props.data,
      features
    };
    props.onEdit({
      updatedData,
      editType: 'setTentativeRectangle',
      editContext: {
        position: event.mapCoords
      }
    });
  }
}
