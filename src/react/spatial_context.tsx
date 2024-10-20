import { createContext, useContext, useEffect, useMemo, useState } from "react";
import type RangeDimension from "../datastore/RangeDimension";
import type { BaseReactChart } from "./components/BaseReactChart";
import { useRegionScale, useScatterplotLayer } from "./scatter_state";
import { Matrix4 } from "@math.gl/core";
import { CompositeMode, EditableGeoJsonLayer, type GeoJsonEditMode } from "@deck.gl-community/editable-layers";
import type { FeatureCollection, Geometry, Position } from '@turf/helpers';
import { getVivId } from "./components/avivatorish/MDVivViewer";
import { useChartID, useRangeDimension2D } from "./hooks";
/*****
 * Persisting some properties related to SelectionOverlay in "SpatialAnnotationProvider"... >>subject to change<<.
 * Not every type of chart will have a range dimension, and not every chart will have a selection overlay etc.
 * Needs will also get more complex, and now we have a somewhat convoluted way of doing something simple.
 * Probably going to be a zustand store in not too long.
 */

type P = [number, number];
type RangeState = {
    rangeDimension: RangeDimension;
    selectionFeatureCollection: FeatureCollection;
    editableLayer: EditableGeoJsonLayer;
    selectionMode: GeoJsonEditMode;
    setSelectionMode: (mode: GeoJsonEditMode) => void;
    modelMatrix: Matrix4;
};
type MeasureState = {
    startPixels: P;
    setStart: (p: P) => void;
    endPixels: P;
    setEnd: (p: P) => void;
};
type SpatialAnnotationState = {
    rectRange: RangeState;
    measure: MeasureState;
};

// Could more usefully be thought of as SpatialContext?
const SpatialAnnotationState = createContext<SpatialAnnotationState>(undefined);

const getEmptyFeatureCollection = () => ({
    type: "FeatureCollection",
    features: []
} as FeatureCollection);

function useSelectionCoords(selection: FeatureCollection) {
    const feature = selection.features[0];
    const coords = useMemo(() => {
        if (!feature) return [];
        //these casts are unsafe in a general sense, but should be ok in our editor.
        //?we could set a property in the feature to say when it's simple AABB?
        //^^ need to be careful about managing that property.
        const geometry = feature.geometry as Geometry;
        const raw = geometry.coordinates as Position[][];
        return raw[0];
    }, [feature]);
    return coords as [number, number][];
}

/** for this to be more useful as a hook will depend on state/context... */
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


function useCreateRange(chart: BaseReactChart<any>) {
    const id = useChartID();
    const { modelMatrix } = useScatterModelMatrix();
    // consider making selectionFeatureCollection part of config, so it can be persisted
    // !nb as of this writing, the scale of these features will be wrong if there is useRegionScale() / modelMatrix that compensates for image being different to 'regions'
    // so when we are persisting editable-geojson in a way that will be used elsewhere we need to address that.
    const [selectionFeatureCollection, setSelectionFeatureCollection] = useState<FeatureCollection>(getEmptyFeatureCollection());
    const [selectionMode, setSelectionMode] = useState<GeoJsonEditMode>(new CompositeMode([]));
    const { filterPoly, removeFilter, rangeDimension } = useRangeDimension2D();
    const coords = useSelectionCoords(selectionFeatureCollection);
    
    useEffect(() => {
        console.log("pending different way of managing resetButton?");
        chart.removeFilter = () => {
            setSelectionFeatureCollection(getEmptyFeatureCollection());
        }
    });
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
    const editableLayer = useMemo(() => {
        return new EditableGeoJsonLayer({
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
            onEdit: ({ updatedData }) => {
                // console.log("onEdit", editType, updatedData);
                const feature = updatedData.features.pop();
                updatedData.features = [feature];
                setSelectionFeatureCollection(updatedData);
            },
            onHover(pickingInfo, event) {
                if ((pickingInfo as any).featureType === "points") return;
                // -- try to avoid selecting invisible features etc - refer to notes in aosta prototype
                setSelectedFeatureIndexes(pickingInfo.index !== -1 ? [pickingInfo.index] : []);
            },
        })
    }, [selectionFeatureCollection, selectionMode, id, selectedFeatureIndexes]);
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
function useCreateSpatialAnnotationState(chart: BaseReactChart<any>) {
    // should we use zustand for this state?
    // doesn't matter too much as it's just used once by SpatialAnnotationProvider
    // consider for project-wide annotation stuff as opposed to ephemeral selections
    const rectRange = useCreateRange(chart);
    const measure = useCreateMeasure();
    return { rectRange, measure };
}

export function SpatialAnnotationProvider({
    chart,
    children,
}: { chart: BaseReactChart<any> } & React.PropsWithChildren) {
    const annotationState = useCreateSpatialAnnotationState(chart);
    return (
        <SpatialAnnotationState.Provider value={annotationState}>
            {children}
        </SpatialAnnotationState.Provider>
    );
}

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

/** work in progress... very much unstable return type etc, but starting to make use 
 * and hopefully refactor into something coherent soon.
 */
export function useSpatialLayers() {
    const { rectRange } = useContext(SpatialAnnotationState);
    const scatterProps = useScatterplotLayer(rectRange.modelMatrix);
    const { getTooltip } = scatterProps;
    // const layers = [rectRange.polygonLayer, scatterplotLayer]; /// should probably be in a CompositeLayer?
    return { 
        getTooltip, scatterProps, selectionLayer: rectRange.editableLayer,
        selectionProps: rectRange
    };
}
