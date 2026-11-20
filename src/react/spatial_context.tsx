import type { Viewport } from "@deck.gl/core";
import type {
    GeoJsonEditMode,
    Position,
    SimpleFeatureCollection,
    SimpleGeometry,
} from "@deck.gl-community/editable-layers";
import { Matrix4 } from "@math.gl/core";
import { action, reaction, toJS } from "mobx";
import { observer } from "mobx-react-lite";
import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useState } from "react";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import { loadColumn } from "@/dataloaders/DataLoaderUtil";
import {
    panMode,
    RectangleMode,
    SELECTION_COORDINATE_SYSTEM,
    selectionFromEdit,
} from "@/editable-layers/selectionModes";
import { chartDragSetting } from "@/lib/browserSettings";
import { MonkeyPatchEditableGeoJsonLayer } from "@/lib/deckMonkeypatch";
import type RangeDimension from "../datastore/RangeDimension";
import { getVivId } from "./components/avivatorish/MDVivViewer";
import type { Tool } from "./components/SelectionOverlay";
import type { HoveredContour } from "./contour_state";
import { getEmptyFeatureCollection } from "./deck_state";
import { DEFAULT_GATE_COLOR, SELECTION_FILL_COLOR, SELECTION_LINE_COLOR } from "./gates/gateUtils";
import { useGateManager } from "./gates/useGateManager";
import { type FilterPolyRegionOpts, useChartID, useConfig, useParamColumns, useRangeDimension2D } from "./hooks";
import {
    DEFAULT_AXIS_SCALES,
    type OrbitSelectionView,
    orbitSelectionViewOf,
    screenProjectionMatrix,
    validAxisScales,
} from "./orbitView";
import { type ScatterPlotConfig, type ScatterPlotConfig3D, useRegionScale, useScatterplotLayer } from "./scatter_state";

/*****
 * Persisting some properties related to SelectionOverlay in "SpatialAnnotationProvider"... >>subject to change<<.
 * Not every type of chart will have a range dimension, and not every chart will have a selection overlay etc.
 * Needs will also get more complex, and now we have a somewhat convoluted way of doing something simple.
 * Probably going to be a zustand store in not too long.
 */

export type P = [number, number];
export type RangeState = {
    /** `undefined` until the owning effect has created it. */
    rangeDimension: RangeDimension | undefined;
    selectionFeatureCollection: SimpleFeatureCollection;
    editableLayer: MonkeyPatchEditableGeoJsonLayer;
    editingGateId: string | null;
    setEditingGateId: (id: string | null) => void;
    setSelectionFeatureCollection: (newSelection: SimpleFeatureCollection) => void;
    selectionMode: GeoJsonEditMode;
    setSelectionMode: (mode: GeoJsonEditMode) => void;
    selectedTool: Tool;
    setSelectedTool: (tool: Tool) => void;
    modelMatrix: Matrix4;
    /**
     * Whether a drag pressed at a point should go to another layer rather than draw a selection: a
     * gate label there, which drags to move (see `useGateLayers`). Null for none.
     */
    yieldDragRef: { current: ((point: [number, number], viewport: Viewport) => boolean) | null };
};
export type MeasureState = {
    startPixels: P;
    setStart: (p: P) => void;
    endPixels: P;
    setEnd: (p: P) => void;
};
export type SpatialAnnotationState = {
    rectRange: RangeState;
    measure: MeasureState;
    scatterProps: ReturnType<typeof useScatterplotLayer>;
    hoveredContour?: HoveredContour | null;
    setHoveredContour?: (contour: HoveredContour | null) => void;
};

// Could more usefully be thought of as SpatialContext?
const SpatialAnnotationState = createContext<SpatialAnnotationState>(undefined as any);

function useSelectionCoords(selection: SimpleFeatureCollection) {
    // where should we keep this in config for persisting?
    const feature = selection?.features[0];
    const coords = useMemo(() => {
        if (!feature) return [];
        //these casts are unsafe in a general sense, but should be ok in our editor.
        //?we could set a property in the feature to say when it's simple AABB?
        //^^ need to be careful about managing that property.
        const geometry = feature.geometry as SimpleGeometry;
        const raw = geometry.coordinates as Position[][];
        //! without toJS, this can be orders of magnitude slower than before - careful with that mobx, eugene...
        // still adds significant overhead - may need a different strategy for critical/hot paths.
        return toJS(raw[0]);
    }, [feature]);
    return coords as [number, number][];
}

/** The tool a chart starts with (see `chartDragSetting`): Pan moves the view, Rectangle draws. */
function defaultTool(): Tool {
    return chartDragSetting.get() === "select" ? "Rectangle" : "Pan";
}

function defaultToolMode(): GeoJsonEditMode {
    return chartDragSetting.get() === "select" ? new RectangleMode() : panMode();
}

/** A chart that draws in 3D says what its camera is, for a selection drawn on it to record. */
function getOrbitSelectionView(chart: BaseChart<any>): OrbitSelectionView | null {
    if (!("getOrbitSelectionView" in chart) || typeof chart.getOrbitSelectionView !== "function") return null;
    return orbitSelectionViewOf({ view: chart.getOrbitSelectionView() }) ?? null;
}

/** In 3D, the selection with the camera it was drawn in (see `projectRowsToScreen`). */
function withSelectionView(chart: BaseChart<any>, selection: SimpleFeatureCollection): SimpleFeatureCollection {
    if (chart.config.dimension !== "3d") return selection;
    const view = getOrbitSelectionView(chart);
    return {
        ...selection,
        features: selection.features.map((feature) => ({
            ...feature,
            properties: { ...feature.properties, view: view ?? undefined },
        })),
    };
}

/**
 * for this to be more useful as a hook will depend on state/context...
 * and with proper spatialdata support we should review this so we have proper coordinate system.
 * hopefully we can move some things around for better HMR editing DX as well.
 */
function useScatterModelMatrix() {
    const scale = useRegionScale();
    const config = useConfig<ScatterPlotConfig & Partial<ScatterPlotConfig3D>>();
    const is3d = config.dimension === "3d";
    // a 3D chart can stretch each axis (`axis_scales`), and its axes need not be the region's x / y
    const [ax, ay, az] = is3d ? validAxisScales(config.axis_scales) : DEFAULT_AXIS_SCALES;
    const s = is3d ? 1 : 1 / scale;
    const modelMatrix = useMemo(() => new Matrix4().scale([s * ax, s * ay, s * az]), [s, ax, ay, az]);
    return { modelMatrix };
}

function useCreateRange(chart: BaseChart<ScatterPlotConfig & BaseConfig>) {
    const id = useChartID();
    const { modelMatrix } = useScatterModelMatrix();
    const gateManager = useGateManager();
    // making selectionFeatureCollection part of config, so it can be persisted
    // !nb as of this writing, the scale of these features will be wrong if there is useRegionScale() / modelMatrix that compensates for image being different to 'regions'
    // so when we are persisting editable-geojson in a way that will be used elsewhere we need to address that later.
    //const [selectionFeatureCollection, setSelectionFeatureCollection] = useState<SimpleFeatureCollection>(getEmptyFeatureCollection());
    // config marks this `observable.ref` (chartConfigUtils): deep observability converted the whole
    // collection on every edit and made deck.gl read its coordinates through proxies
    const { selectionFeatureCollection } = chart.config;
    const setSelectionFeatureCollection = useCallback(
        action((newSelection: SimpleFeatureCollection) => {
            chart.config.selectionFeatureCollection = newSelection;
        }),
        [],
    );
    // what a plain drag does until a tool is picked: move the view (Pan) or draw a rectangle, as
    // chosen in this browser (`chartDragSetting`); ⌘-drag moves the view either way (`viewGestureGate`)
    const [selectionMode, setSelectionMode] = useState<GeoJsonEditMode>(() => defaultToolMode());
    const [selectedTool, setSelectedToolState] = useState<Tool>(() => defaultTool());
    const toolPicked = useRef(false);
    const setSelectedTool = useCallback((tool: Tool) => {
        toolPicked.current = true;
        setSelectedToolState(tool);
    }, []);
    // a chart whose tool nobody has picked follows the setting when it changes (as when it is first chosen)
    useEffect(
        () =>
            reaction(
                () => chartDragSetting.get(),
                () => {
                    if (toolPicked.current) return;
                    setSelectionMode(defaultToolMode());
                    setSelectedToolState(defaultTool());
                },
            ),
        [],
    );
    const [editingGateId, setEditingGateId] = useState<string | null>(null);
    const yieldDragRef = useRef<RangeState["yieldDragRef"]["current"]>(null);
    const { filterPoly, removeFilter, rangeDimension } = useRangeDimension2D();
    const coords = useSelectionCoords(selectionFeatureCollection);
    const is3d = chart.config.dimension === "3d";
    // a 3D selection is an outline on the screen, and selects what was inside it from where it was drawn
    const selectionProperties = selectionFeatureCollection.features[0]?.properties;
    const selectionView = useMemo(
        () => (is3d ? orbitSelectionViewOf(toJS(selectionProperties)) : undefined),
        [is3d, selectionProperties],
    );
    const [cx, cy, cz] = useParamColumns();
    const filterOnScreen = useCallback(
        (points: [number, number][], view: OrbitSelectionView) => {
            if (!rangeDimension || !cx || !cy || !cz) return;
            const matrix = screenProjectionMatrix(view, modelMatrix);
            // projected and tested on the dimension's filter worker
            rangeDimension.filter("filterScreenPoly", [cx.field, cy.field, cz.field], { matrix, points });
        },
        [rangeDimension, cx, cy, cz, modelMatrix],
    );

    const getFillColor = useCallback((): [number, number, number, number] => {
        if (!editingGateId) return SELECTION_FILL_COLOR; // Neutral grey fill

        const gate = gateManager.gates.get(editingGateId);
        const [r, g, b] = gate?.color ?? DEFAULT_GATE_COLOR;

        return [r, g, b, 45]; // lower opacity of original color
    }, [editingGateId, gateManager]);

    const getLineColor = useCallback((): [number, number, number, number] => {
        if (!editingGateId) return SELECTION_LINE_COLOR; // Dark grey selection line
        const gate = gateManager.gates.get(editingGateId);
        const [r, g, b] = gate?.color ?? DEFAULT_GATE_COLOR;

        // brighten the color to make it look highlighted
        return [Math.min(255, r + 60), Math.min(255, g + 60), Math.min(255, b + 60), 255];
    }, [editingGateId, gateManager]);

    const clearSelection = useCallback(() => {
        setSelectionFeatureCollection(getEmptyFeatureCollection());
        setEditingGateId(null);
    }, [setSelectionFeatureCollection]);

    useEffect(() => {
        chart.removeFilter = clearSelection;
    }, [chart, clearSelection]);

    // "Reset All Filters" clears the dimension directly via `dataStore.removeAllFilters()`
    // without going through `chart.removeFilter`, so the drawn shape needs clearing here too.
    useEffect(() => {
        const ds = chart.dataStore;
        const listenerId = `${chart.config.id}_selection_layer`;
        ds.addListener(listenerId, (type: string, data: unknown) => {
            if (type === "filtered" && data === "all_removed") clearSelection();
        });
        return () => ds.removeListener(listenerId);
    }, [chart, clearSelection]);

    useEffect(() => {
        if (coords.length === 0) {
            removeFilter();
            chart.updateResetButtonVisibility();
            return;
        }
        chart.resetButton.style.display = "inline";

        if (is3d) {
            // a shape drawn before this chart knew its camera (or not in 3D) can't say what it selected
            if (selectionView) filterOnScreen(coords, selectionView);
            else removeFilter();
            return;
        }

        const ds = chart.dataStore;
        const regionField = ds.regions?.region_field;
        // Get the region value of the gate if it exists
        const regionValue = editingGateId ? gateManager.gates.get(editingGateId)?.region : undefined;

        // Guard async callbacks with cancellation flag so only the latest
        // effect run calls the filterPoly thereby avoiding stale polygon filter if
        // the older effect is resolved later
        let cancelled = false;

        if (regionField && regionValue) {
            const regionOptions: FilterPolyRegionOpts = { regionField, regionValue };
            if (ds.columnIndex[regionField]?.data) {
                filterPoly(coords, regionOptions);
            } else {
                // Load the region if it's not loaded
                loadColumn(ds.name, regionField)
                    .then(() => {
                        if (!cancelled) filterPoly(coords, regionOptions);
                    })
                    .catch((error) => {
                        console.error(`Could not load region column '${regionField}':`, error);
                        if (!cancelled) filterPoly(coords); // fallback: polygon-only
                    });
            }
        } else {
            // No active region (or no region metadata): plain polygon filter.
            filterPoly(coords);
        }

        return () => {
            cancelled = true;
        };
    }, [coords, filterPoly, removeFilter, chart, gateManager, editingGateId, is3d, selectionView, filterOnScreen]);
    // The collection holds at most one shape (see `onEdit`), and it is always the one being
    // edited. Selecting it on hover made translate-vs-draw depend on whether a hover pick had
    // landed before the press, which on WebGPU (async picking) it often had not.
    const hasShape = selectionFeatureCollection.features.length > 0;
    const selectedFeatureIndexes = useMemo(() => (hasShape ? [0] : []), [hasShape]);

    // we might be able to pass this to modeConfig, if it knows what to do with it?
    // const outerContainer = useOuterContainer();

    // todo: we need to do something about editing gate color
    //! When the dev tools are open, there is lagging while dragging this layer
    const editableLayer = useMemo(() => {
        return new MonkeyPatchEditableGeoJsonLayer({
            id: `selection_${getVivId(`${id}detail-react`)}`,
            data: selectionFeatureCollection,
            mode: selectionMode,
            yieldDrag: (point, viewport) => yieldDragRef.current?.(point, viewport) ?? false,
            onEmptyClick: clearSelection,
            coordinateSystem: SELECTION_COORDINATE_SYSTEM,
            getFillColor,
            getLineColor,
            getLineWidth: 2,
            lineWidthMinPixels: 2.5,
            getEditHandlePointRadius: 2,
            editHandlePointStrokeWidth: 1,
            editHandleIconSizeScale: 1,
            editHandlePointRadiusMaxPixels: 5,
            getEditHandlePointColor: (h) => {
                const intermediate = h.properties.editHandleType === "intermediate";
                return [0, 0, intermediate ? 200 : 0, 255];
            },
            selectedFeatureIndexes,
            // adding `action` here gets rid of warnings but doesn't help with performance.
            onEdit: action((edit) => {
                const selection = selectionFromEdit(edit);
                if (selection) setSelectionFeatureCollection(withSelectionView(chart, selection));
            }),
            modeConfig: {
                // dragToDraw: true,
                // hopefully there is something here we can pass outerContainer to...
                // or maybe it's enough to remake the layer when it changes?
                // it looks as though editable-layer.js _addEventHandlers() isn't called again when eventManager has been changed.
                // outerContainer
            },
        });
    }, [
        selectionFeatureCollection,
        selectionMode,
        id,
        selectedFeatureIndexes,
        setSelectionFeatureCollection,
        getFillColor,
        getLineColor,
        chart,
        clearSelection,
    ]);

    return {
        editableLayer,
        yieldDragRef,
        rangeDimension,
        selectionFeatureCollection,
        selectionMode,
        setSelectionMode,
        setSelectionFeatureCollection,
        modelMatrix,
        editingGateId,
        setEditingGateId,
        selectedTool,
        setSelectedTool,
    };
}
function useCreateMeasure() {
    const [startPixels, setStart] = useState<P>([0, 0]);
    const [endPixels, setEnd] = useState<P>([0, 0]);
    return { startPixels, setStart, endPixels, setEnd };
}
//add generic that extends SpatialConfig?
function useCreateSpatialAnnotationState(
    chart: BaseChart<any>,
    hoveredContour?: HoveredContour | null,
    separateDensityLayer = false,
) {
    // should we use zustand for this state?
    // doesn't matter too much as it's just used once by SpatialAnnotationProvider
    // consider for project-wide annotation stuff as opposed to ephemeral selections
    const rectRange = useCreateRange(chart);
    const measure = useCreateMeasure();
    const scatterProps = useScatterplotLayer(
        rectRange.modelMatrix,
        hoveredContour,
        rectRange.rangeDimension,
        separateDensityLayer,
    );
    return { rectRange, measure, scatterProps };
}

export const SpatialAnnotationProvider = observer(function SpatialAnnotationProvider({
    chart,
    children,
    hoveredContour,
    setHoveredContour,
    separateDensityLayer,
    //add generic that extends SpatialConfig?
}: {
    chart: BaseChart<any>;
    hoveredContour?: HoveredContour | null;
    setHoveredContour?: (contour: HoveredContour | null) => void;
    /** Give `scatterProps` a `densityLayer` and keep density out of the points layer. */
    separateDensityLayer?: boolean;
} & React.PropsWithChildren) {
    const annotationState = useCreateSpatialAnnotationState(chart, hoveredContour, separateDensityLayer);
    const stateWithHover = { ...annotationState, hoveredContour, setHoveredContour };
    return <SpatialAnnotationState.Provider value={stateWithHover}>{children}</SpatialAnnotationState.Provider>;
});

export function useRange() {
    const range = useContext(SpatialAnnotationState).rectRange;
    if (!range) throw new Error("no range context");
    return range;
}

export function useMeasure() {
    const measure = useContext(SpatialAnnotationState).measure;
    if (!measure) throw new Error("no measure context");
    return measure;
}

/**
 * useSpatialLayers is a hook that provides access to the spatial layers and their properties.
 *
 * This includes the selection layer, scatterplot layer (which is actually an instance
 * of a more complex `SpatialLayer` with custom density/contour rendering), etc.
 *
 * work in progress... unstable return type etc, and subject to future changes.
 */
export function useSpatialLayers() {
    const { rectRange, scatterProps } = useContext(SpatialAnnotationState);
    // this should be in the context, useScatterplotLayer currently called *only* from there...
    // if we want to use it in other places we need to refactor
    // const scatterProps = useScatterplotLayer(rectRange.modelMatrix);
    const { getTooltip } = scatterProps;
    // const layers = [rectRange.polygonLayer, scatterplotLayer]; /// should probably be in a CompositeLayer?
    return {
        getTooltip,
        scatterProps,
        selectionLayer: rectRange.editableLayer,
        selectionProps: rectRange,
    };
}
