import { createContext, useCallback, useContext, useEffect, useMemo, useState } from "react";
import type RangeDimension from "../datastore/RangeDimension";
import { type ScatterPlotConfig, useRegionScale, useScatterplotLayer } from "./scatter_state";
import { Matrix4 } from "@math.gl/core";
import { CompositeMode, type GeoJsonEditMode } from "@deck.gl-community/editable-layers";
import type { FeatureCollection, Geometry, Position } from '@turf/helpers';
import { getVivId } from "./components/avivatorish/MDVivViewer";
import { useChartID, useRangeDimension2D } from "./hooks";
import type BaseChart from "@/charts/BaseChart";
import { observer } from "mobx-react-lite";
import type { BaseConfig } from "@/charts/BaseChart";
import { action, toJS } from "mobx";
import { getEmptyFeatureCollection } from "./deck_state";
import { MonkeyPatchEditableGeoJsonLayer } from "@/lib/deckMonkeypatch";

/*****
 * Persisting some properties related to SelectionOverlay in "SpatialAnnotationProvider"... >>subject to change<<.
 * Not every type of chart will have a range dimension, and not every chart will have a selection overlay etc.
 * Needs will also get more complex, and now we have a somewhat convoluted way of doing something simple.
 * Probably going to be a zustand store in not too long.
 */

export type P = [number, number];
export type RangeState = {
    rangeDimension: RangeDimension;
    selectionFeatureCollection: FeatureCollection;
    editableLayer: MonkeyPatchEditableGeoJsonLayer;
    selectionMode: GeoJsonEditMode;
    setSelectionMode: (mode: GeoJsonEditMode) => void;
    modelMatrix: Matrix4;
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
};

// Could more usefully be thought of as SpatialContext?
const SpatialAnnotationState = createContext<SpatialAnnotationState>(undefined as any);


function useSelectionCoords(selection: FeatureCollection) {
    // where should we keep this in config for persisting?
    const feature = selection?.features[0];
    const coords = useMemo(() => {
        if (!feature) return [];
        //these casts are unsafe in a general sense, but should be ok in our editor.
        //?we could set a property in the feature to say when it's simple AABB?
        //^^ need to be careful about managing that property.
        const geometry = feature.geometry as Geometry;
        const raw = geometry.coordinates as Position[][];
        //! without toJS, this can be orders of magnitude slower than before - careful with that mobx, eugene...
        // still adds significant overhead - may need a different strategy for critical/hot paths.
        return toJS(raw[0]);
    }, [feature]);
    return coords as [number, number][];
}

/** 
 * for this to be more useful as a hook will depend on state/context...
 * and with proper spatialdata support we should review this so we have proper coordinate system.
 * hopefully we can move some things around for better HMR editing DX as well.
 */
function useScatterModelMatrix() {
    const scale = useRegionScale();
    const s = 1 / scale;
    const [modelMatrix, setModelMatrix] = useState(new Matrix4().scale(s));
    useEffect(() => {
        const m = new Matrix4().scale(s);
        setModelMatrix(m);
    }, [s]);
    return { modelMatrix, setModelMatrix };
}


function useCreateRange(chart: BaseChart<ScatterPlotConfig & BaseConfig>) {
    const id = useChartID();
    const { modelMatrix } = useScatterModelMatrix();
    // making selectionFeatureCollection part of config, so it can be persisted
    // !nb as of this writing, the scale of these features will be wrong if there is useRegionScale() / modelMatrix that compensates for image being different to 'regions'
    // so when we are persisting editable-geojson in a way that will be used elsewhere we need to address that later.
    //const [selectionFeatureCollection, setSelectionFeatureCollection] = useState<FeatureCollection>(getEmptyFeatureCollection());
    //! it appears to drastically affect performance having this in mobx config - why?
    const { selectionFeatureCollection } = chart.config;
    const setSelectionFeatureCollection = useCallback(action((newSelection: FeatureCollection) => {
        chart.config.selectionFeatureCollection = newSelection;
    }), []);
    const [selectionMode, setSelectionMode] = useState<GeoJsonEditMode>(new CompositeMode([]));
    const { filterPoly, removeFilter, rangeDimension } = useRangeDimension2D();
    const coords = useSelectionCoords(selectionFeatureCollection);

    useEffect(() => {
        console.log("pending different way of managing resetButton?");
        chart.removeFilter = () => {
            setSelectionFeatureCollection(getEmptyFeatureCollection());
        }
    }, [chart, setSelectionFeatureCollection]);
    useEffect(() => {
        if (coords.length === 0) {
            chart.resetButton.style.display = "none";
            removeFilter();
            return;
        }
        // todo - consider whether the shape is a simple AABB & apply faster filter if so.
        //rangeDimension.filterPoly(coords, [cols[0], cols[1]]); //this doesn't notify 🙄
        chart.resetButton.style.display = "inline";
        filterPoly(coords);
    }, [coords, filterPoly, removeFilter, chart]);
    const [selectedFeatureIndexes, setSelectedFeatureIndexes] = useState<number[]>([]);
    // we might be able to pass this to modeConfig, if it knows what to do with it?
    // const outerContainer = useOuterContainer();
    const editableLayer = useMemo(() => {
        return new MonkeyPatchEditableGeoJsonLayer({
            id: `selection_${getVivId(`${id}detail-react`)}`,
            data: selectionFeatureCollection as any,
            mode: selectionMode,
            getFillColor: [140, 140, 140, 50],
            getLineColor: [255, 255, 255, 200],
            getLineWidth: 1,
            getEditHandlePointRadius: 2,
            editHandlePointStrokeWidth: 1,
            editHandleIconSizeScale: 1,
            editHandlePointRadiusMaxPixels: 5,
            getEditHandlePointColor: h => {
                const intermediate = h.properties.editHandleType === "intermediate";
                return [0, 0, intermediate ? 200 : 0, 255]
            },
            lineWidthMinPixels: 1,
            selectedFeatureIndexes,
            // adding `action` here gets rid of warnings but doesn't help with performance.
            onEdit: action(({ updatedData }) => {
                // console.log("onEdit", editType, updatedData);
                const feature = updatedData.features.at(-1);
                // updatedData.features = [feature];
                setSelectionFeatureCollection({
                    ...updatedData,
                    features: feature ? [feature] : []
                });
            }),
            onHover(pickingInfo, event) {
                if ((pickingInfo as any).featureType === "points") return;
                // -- try to avoid selecting invisible features etc - refer to notes in aosta prototype
                setSelectedFeatureIndexes(pickingInfo.index !== -1 ? [pickingInfo.index] : []);
            },
            modeConfig: {
                // dragToDraw: true,
                // hopefully there is something here we can pass outerContainer to...
                // or maybe it's enough to remake the layer when it changes?
                // it looks as though editable-layer.js _addEventHandlers() isn't called again when eventManager has been changed.
                // outerContainer
            }
        })
    }, [selectionFeatureCollection, selectionMode, id, selectedFeatureIndexes,
        setSelectionFeatureCollection
    ]);
    return {
        editableLayer,
        rangeDimension,
        selectionFeatureCollection,
        selectionMode,
        setSelectionMode,
        modelMatrix
    };
}
function useCreateMeasure() {
    const [startPixels, setStart] = useState<P>([0, 0]);
    const [endPixels, setEnd] = useState<P>([0, 0]);
    return { startPixels, setStart, endPixels, setEnd };
}
//add generic that extends SpatialConfig?
function useCreateSpatialAnnotationState(chart: BaseChart<any>) {
    // should we use zustand for this state?
    // doesn't matter too much as it's just used once by SpatialAnnotationProvider
    // consider for project-wide annotation stuff as opposed to ephemeral selections
    const rectRange = useCreateRange(chart);
    const measure = useCreateMeasure();
    const scatterProps = useScatterplotLayer(rectRange.modelMatrix);
    return { rectRange, measure, scatterProps };
}

export const SpatialAnnotationProvider = observer(function SpatialAnnotationProvider({
    chart,
    children,
    //add generic that extends SpatialConfig?
}: { chart: BaseChart<any> } & React.PropsWithChildren) {
    const annotationState = useCreateSpatialAnnotationState(chart);
    return (
        <SpatialAnnotationState.Provider value={annotationState}>
            {children}
        </SpatialAnnotationState.Provider>
    );
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
        getTooltip, scatterProps, selectionLayer: rectRange.editableLayer,
        selectionProps: rectRange
    };
}
