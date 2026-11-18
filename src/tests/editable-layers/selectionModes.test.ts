import {
    type DraggingEvent,
    DrawPolygonByDraggingMode,
    type EditAction,
    type FeatureCollection,
    fromDeckCoordinateSystem,
    type GeoJsonEditMode,
    type ModeProps,
    type Pick,
    type Position,
    type SimpleFeatureCollection,
} from "@deck.gl-community/editable-layers";
import { describe, expect, it, vi } from "vitest";
import {
    EditMode,
    PolygonMode,
    RectangleMode,
    SELECTION_COORDINATE_SYSTEM,
    selectionFromEdit,
} from "@/editable-layers/selectionModes";

// `EditableGeoJsonLayer` forwards each mjolnir event to `getActiveMode().handleX(event, props)`,
// so driving the mode directly is the same path a real click / double-click takes.

type XY = [number, number];

function emptyCollection(): FeatureCollection {
    return { type: "FeatureCollection", features: [] };
}

function squareCollection(size = 10): SimpleFeatureCollection {
    return {
        type: "FeatureCollection",
        features: [
            {
                type: "Feature",
                properties: {},
                geometry: {
                    type: "Polygon",
                    coordinates: [
                        [
                            [0, 0],
                            [size, 0],
                            [size, size],
                            [0, size],
                            [0, 0],
                        ],
                    ],
                },
            },
        ],
    };
}

function modeProps(data: FeatureCollection, selectedIndexes: number[] = []) {
    const onEdit = vi.fn<(action: EditAction<FeatureCollection>) => void>();
    const props: ModeProps<FeatureCollection> = {
        data,
        modeConfig: {},
        selectedIndexes,
        cursor: null,
        lastPointerMoveEvent: pointerEvent([0, 0]),
        // what EditableGeoJsonLayer.getModeProps derives from the selection layer's coordinateSystem
        coordinateSystem: fromDeckCoordinateSystem(SELECTION_COORDINATE_SYSTEM),
        onEdit,
        onUpdateCursor: () => {},
    };
    return { props, onEdit };
}

/** screen and map coordinates coincide here; `downAt` / `downPicks` only matter for drags */
function pointerEvent(xy: XY, picks: Pick[] = [], downAt: XY = xy, downPicks: Pick[] | null = null) {
    return {
        picks,
        screenCoords: xy,
        mapCoords: xy,
        sourceEvent: null,
        pointerDownPicks: downPicks,
        pointerDownScreenCoords: downAt,
        pointerDownMapCoords: downAt,
        cancelPan: () => {},
    } satisfies DraggingEvent;
}

function outerRing(data: FeatureCollection | undefined): Position[] {
    const geometry = data?.features[0]?.geometry;
    if (geometry?.type !== "Polygon") throw new Error(`expected a Polygon, got ${geometry?.type}`);
    return geometry.coordinates[0];
}

function featurePick(index: number): Pick {
    return { index, isGuide: false, object: squareCollection().features[index] };
}

function existingHandlePick(position: Position, positionIndexes: number[]): Pick {
    return {
        index: 0,
        isGuide: true,
        object: {
            type: "Feature",
            properties: {
                guideType: "editHandle",
                editHandleType: "existing",
                featureIndex: 0,
                positionIndexes,
            },
            geometry: { type: "Point", coordinates: position },
        },
    };
}

/** the drawing-in-progress vertex handles DrawPolygonMode renders as guides */
function tentativeHandlePick(position: Position, positionIndex: number): Pick {
    return {
        index: positionIndex,
        isGuide: true,
        object: {
            type: "Feature",
            properties: {
                guideType: "editHandle",
                editHandleType: "existing",
                featureIndex: -1,
                positionIndexes: [positionIndex],
            },
            geometry: { type: "Point", coordinates: position },
        },
    };
}

function editTypes(onEdit: ReturnType<typeof modeProps>["onEdit"]) {
    return onEdit.mock.calls.map(([action]) => action.editType);
}

describe("selection edit modes", () => {
    // A double-click is recognised as one `dblclick`, not as clicks, so its position is the last vertex.
    it("Polygon tool adds the double-clicked point as the last vertex and completes the polygon", () => {
        const mode: GeoJsonEditMode = new PolygonMode();
        const { props, onEdit } = modeProps(emptyCollection());
        for (const p of [
            [0, 0],
            [10, 0],
            [10, 10],
        ] satisfies XY[]) {
            mode.handleClick(pointerEvent(p), props);
        }
        mode.handleDoubleClick(pointerEvent([0, 10]), props);

        const added = onEdit.mock.calls.find(([a]) => a.editType === "addFeature")?.[0];
        expect(added?.updatedData.features).toHaveLength(1);
        expect(outerRing(added?.updatedData)).toEqual([
            [0, 0],
            [10, 0],
            [10, 10],
            [0, 10],
            [0, 0],
        ]);
    });

    it("Polygon tool completes a triangle when the third vertex is double-clicked", () => {
        const mode: GeoJsonEditMode = new PolygonMode();
        const { props, onEdit } = modeProps(emptyCollection());
        for (const p of [
            [0, 0],
            [10, 0],
        ] satisfies XY[]) {
            mode.handleClick(pointerEvent(p), props);
        }
        mode.handleDoubleClick(pointerEvent([10, 10]), props);

        const added = onEdit.mock.calls.find(([a]) => a.editType === "addFeature")?.[0];
        expect(outerRing(added?.updatedData)).toEqual([
            [0, 0],
            [10, 0],
            [10, 10],
            [0, 0],
        ]);
    });

    it("Polygon tool doesn't duplicate the last vertex when its handle is double-clicked", () => {
        const mode: GeoJsonEditMode = new PolygonMode();
        const { props, onEdit } = modeProps(emptyCollection());
        for (const p of [
            [0, 0],
            [10, 0],
            [10, 10],
        ] satisfies XY[]) {
            mode.handleClick(pointerEvent(p), props);
        }
        // nearby rather than exact: the handle is picked within a pixel radius
        mode.handleDoubleClick(pointerEvent([10.01, 10], [tentativeHandlePick([10, 10], 2)]), props);

        const added = onEdit.mock.calls.find(([a]) => a.editType === "addFeature")?.[0];
        expect(outerRing(added?.updatedData)).toEqual([
            [0, 0],
            [10, 0],
            [10, 10],
            [0, 0],
        ]);
    });

    it("Polygon tool keeps drawing if a double-click would leave fewer than three vertices", () => {
        const mode: GeoJsonEditMode = new PolygonMode();
        const { props, onEdit } = modeProps(emptyCollection());
        mode.handleClick(pointerEvent([0, 0]), props);
        mode.handleDoubleClick(pointerEvent([10, 0]), props);
        expect(editTypes(onEdit)).not.toContain("addFeature");

        mode.handleDoubleClick(pointerEvent([10, 10]), props);
        const added = onEdit.mock.calls.find(([a]) => a.editType === "addFeature")?.[0];
        expect(outerRing(added?.updatedData)).toEqual([
            [0, 0],
            [10, 0],
            [10, 10],
            [0, 0],
        ]);
    });

    it("Modify tool removes a vertex on double-click", () => {
        const mode: GeoJsonEditMode = new EditMode();
        const { props, onEdit } = modeProps(squareCollection(), [0]);
        mode.handleDoubleClick(pointerEvent([10, 0], [existingHandlePick([10, 0], [0, 1])]), props);

        expect(editTypes(onEdit)).toEqual(["removePosition"]);
        const [{ updatedData }] = onEdit.mock.calls[0];
        expect(outerRing(updatedData)).not.toContainEqual([10, 0]);
    });

    it("Modify tool still removes a vertex on single click", () => {
        const mode: GeoJsonEditMode = new EditMode();
        const { props, onEdit } = modeProps(squareCollection(), [0]);
        mode.handleClick(pointerEvent([10, 0], [existingHandlePick([10, 0], [0, 1])]), props);

        expect(editTypes(onEdit)).toEqual(["removePosition"]);
    });

    it("Modify tool translates a selection by exactly the drag delta in data coordinates", () => {
        // scatter / spatial coordinates are cartesian, and often well outside lng/lat ranges
        const mode: GeoJsonEditMode = new EditMode();
        const { props, onEdit } = modeProps(squareCollection(100), [0]);
        const from: XY = [50, 50];
        const to: XY = [80, 60];
        const picks = [featurePick(0)];
        mode.handlePointerMove(pointerEvent(from, picks), props);
        mode.handleStartDragging(pointerEvent(from, picks, from, picks), props);
        mode.handleDragging(pointerEvent(to, picks, from, picks), props);
        mode.handleStopDragging(pointerEvent(to, picks, from, picks), props);

        const last = onEdit.mock.calls.at(-1)?.[0];
        expect(last?.editType).toBe("translated");
        const ring = outerRing(last?.updatedData);
        const expected = [
            [30, 10],
            [130, 10],
            [130, 110],
            [30, 110],
            [30, 10],
        ];
        ring.forEach((p, i) => {
            expect(p[0]).toBeCloseTo(expected[i][0], 6);
            expect(p[1]).toBeCloseTo(expected[i][1], 6);
        });
    });
});

/**
 * The selection layer only receives a committed shape as its `data` once React re-renders, so
 * events in between are handled with the data it already holds. Model that: `props.data` stays
 * at what the layer was last given, while edits go through the chart's handler.
 */
function staleSelectionLayer<D extends FeatureCollection>(data: D) {
    const onEdit = (edit: EditAction<FeatureCollection>) => {
        const next = selectionFromEdit(edit);
        if (next) layer.selection = next;
    };
    const props: ModeProps<D> = {
        data,
        modeConfig: {},
        selectedIndexes: [],
        cursor: null,
        lastPointerMoveEvent: pointerEvent([0, 0]),
        coordinateSystem: fromDeckCoordinateSystem(SELECTION_COORDINATE_SYSTEM),
        onEdit,
        onUpdateCursor: () => {},
    };
    const layer: {
        selection: FeatureCollection;
        props: ModeProps<D>;
        /** editable-layers types some handlers' props as any `FeatureCollection` */
        anyProps: ModeProps<FeatureCollection>;
    } = { selection: data, props, anyProps: { ...props, onEdit } };
    return layer;
}

describe("selection layer edits", () => {
    it("a hover move after finishing a freehand shape keeps the new shape", () => {
        // the Freehand tool's mode
        const mode = new DrawPolygonByDraggingMode();
        const previous = squareCollection();
        const layer = staleSelectionLayer(previous);
        const down: XY = [20, 20];
        mode.handleStartDragging(pointerEvent(down), layer.props);
        for (const p of [
            [30, 20],
            [30, 30],
            [20, 30],
        ] satisfies XY[]) {
            mode.handleDragging(pointerEvent(p, [], down), layer.props);
        }
        // drawing in progress is the layer's own guides, not the selection
        expect(layer.selection.features).toEqual(previous.features);

        mode.handleStopDragging(pointerEvent([21, 29], [], down), layer.props);
        const drawn = layer.selection;
        expect(outerRing(drawn)).toEqual([
            [30, 20],
            [30, 30],
            [20, 30],
            [21, 29],
            [30, 20],
        ]);

        // the mouse moves on before the layer has been given `drawn`
        mode.handlePointerMove(pointerEvent([22, 28]), layer.anyProps);
        expect(layer.selection).toBe(drawn);
    });

    it("a hover move after finishing a polygon keeps the new shape", () => {
        const mode: GeoJsonEditMode = new PolygonMode();
        const previous = squareCollection();
        const layer = staleSelectionLayer<FeatureCollection>(previous);
        for (const p of [
            [20, 20],
            [30, 20],
            [30, 30],
        ] satisfies XY[]) {
            mode.handleClick(pointerEvent(p), layer.props);
        }
        mode.handleDoubleClick(pointerEvent([20, 30]), layer.props);
        const drawn = layer.selection;
        expect(drawn).not.toBe(previous);

        mode.handlePointerMove(pointerEvent([22, 28]), layer.props);
        expect(layer.selection).toBe(drawn);
    });

    it("a rectangle dragged before the layer re-renders is the one committed", () => {
        // a quick drag can end before the layer has been given any tentative rectangle
        const mode: GeoJsonEditMode = new RectangleMode();
        const previous = squareCollection();
        const layer = staleSelectionLayer<FeatureCollection>(previous);
        const down: XY = [20, 20];
        mode.handleStartDragging(pointerEvent(down), layer.props);
        mode.handleDragging(pointerEvent([25, 25], [], down), layer.props);
        mode.handleStopDragging(pointerEvent([40, 30], [], down), layer.props);

        const drawn = layer.selection;
        expect(drawn.features).toHaveLength(1);
        expect(drawn.features[0].properties).toEqual({});
        expect(outerRing(drawn)).toEqual([
            [20, 20],
            [40, 20],
            [40, 30],
            [20, 30],
            [20, 20],
        ]);
        // and the previous selection isn't modified in place
        expect(previous).toEqual(squareCollection());

        mode.handlePointerMove(pointerEvent([41, 31]), layer.props);
        expect(layer.selection).toBe(drawn);
    });
});
