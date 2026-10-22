import {
    type DraggingEvent,
    type EditAction,
    type FeatureCollection,
    type GeoJsonEditMode,
    type ModeProps,
    type Pick,
    type Position,
    fromDeckCoordinateSystem,
} from "@deck.gl-community/editable-layers";
import { describe, expect, it, vi } from "vitest";
import { EditMode, PolygonMode, SELECTION_COORDINATE_SYSTEM } from "@/editable-layers/selectionModes";

// `EditableGeoJsonLayer` forwards each mjolnir event to `getActiveMode().handleX(event, props)`,
// so driving the mode directly is the same path a real click / double-click takes.

type XY = [number, number];

function emptyCollection(): FeatureCollection {
    return { type: "FeatureCollection", features: [] };
}

function squareCollection(size = 10): FeatureCollection {
    return {
        type: "FeatureCollection",
        features: [
            {
                type: "Feature",
                properties: {},
                geometry: {
                    type: "Polygon",
                    coordinates: [[[0, 0], [size, 0], [size, size], [0, size], [0, 0]]],
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

function editTypes(onEdit: ReturnType<typeof modeProps>["onEdit"]) {
    return onEdit.mock.calls.map(([action]) => action.editType);
}

describe("selection edit modes", () => {
    it("Polygon tool completes the polygon on double-click", () => {
        const mode: GeoJsonEditMode = new PolygonMode();
        const { props, onEdit } = modeProps(emptyCollection());
        for (const p of [[0, 0], [10, 0], [10, 10], [0, 10]] satisfies XY[]) {
            mode.handleClick(pointerEvent(p), props);
        }
        mode.handleDoubleClick(pointerEvent([0, 10]), props);

        expect(editTypes(onEdit)).toContain("addFeature");
        const added = onEdit.mock.calls.find(([a]) => a.editType === "addFeature")?.[0];
        expect(added?.updatedData.features).toHaveLength(1);
        expect(outerRing(added?.updatedData)).toHaveLength(5);
    });

    it("Polygon tool keeps drawing if double-clicked before it has three vertices", () => {
        const mode: GeoJsonEditMode = new PolygonMode();
        const { props, onEdit } = modeProps(emptyCollection());
        for (const p of [[0, 0], [10, 0]] satisfies XY[]) {
            mode.handleClick(pointerEvent(p), props);
        }
        mode.handleDoubleClick(pointerEvent([10, 0]), props);
        expect(editTypes(onEdit)).not.toContain("addFeature");

        mode.handleClick(pointerEvent([10, 10]), props);
        mode.handleDoubleClick(pointerEvent([10, 10]), props);
        const added = onEdit.mock.calls.find(([a]) => a.editType === "addFeature")?.[0];
        expect(outerRing(added?.updatedData)).toHaveLength(4);
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
        const expected = [[30, 10], [130, 10], [130, 110], [30, 110], [30, 10]];
        ring.forEach((p, i) => {
            expect(p[0]).toBeCloseTo(expected[i][0], 6);
            expect(p[1]).toBeCloseTo(expected[i][1], 6);
        });
    });
});
