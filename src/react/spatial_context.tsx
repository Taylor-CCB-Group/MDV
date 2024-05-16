import { createContext, useContext, useEffect, useMemo, useState } from "react";
import type RangeDimension from "../datastore/RangeDimension";
import { BaseReactChart } from "./components/BaseReactChart";
import { PolygonLayer } from '@deck.gl/layers';

/*****
 * Persisting some properties related to SelectionOverlay in "SpatialAnnotationProvider"... >>subject to change<<.
 * Not every type of chart will have a range dimension, and not every chart will have a selection overlay etc.
 * Needs will also get more complex, and now we have a somewhat convoluted way of doing something simple.
 * Probably going to be a zustand store in not too long.
 */

type P = [number, number];
type RefP = React.MutableRefObject<P>;
type RangeState = {
    polygonLayer: PolygonLayer;
    rangeDimension: RangeDimension;
    start: P; setStart: (p: P) => void; startRef: RefP;
    end: P; setEnd: (p: P) => void; endRef: RefP;
};
type MeasureState = {
    startPixels: P; setStart: (p: P) => void;
    endPixels: P; setEnd: (p: P) => void;
}
type PolygonRegion = {
    coords: P[];
}
type SpatialAnnotationState = {
    rectRange: RangeState;
    measure: MeasureState;
}

// Could more usefully be thought of as SpatialContext?
const SpatialAnnotationState = createContext<SpatialAnnotationState>(undefined);

function useCreateRange(chart: BaseReactChart<any>) {
    const ds = chart.dataStore;
    // tried simpler `rangeDimesion = useMemo(...)`, but it can lead to non-destroyed rangeDimensions with HMR.
    const [rangeDimension, setRangeDimension] = useState<RangeDimension>(undefined);
    const [start, setStartX] = useState<P>([0, 0]);
    const [end, setEndX] = useState<P>([0, 0]);
    const polygonLayer = useMemo(() => {
        const data = [start, [end[0], start[1]], end, [start[0], end[1]]];
        const layer = new PolygonLayer({
            id: 'PolygonLayer', //todo: may want to be viv-like?
            data,

            // getPolygon: d => d.contour,
            // getElevation: d => d.population / d.area / 10,
            getFillColor: [140, 140, 140],
            getLineColor: [255, 255, 255],
            getLineWidth: 20,
            lineWidthMinPixels: 1,
            fillOpacity: 0.2,
            pickable: true
        });
        return layer;
    }, [start, end])
    // still not sure I want these refs
    const startRef = useMemo(() => ({ current: start }), [start]);
    const endRef = useMemo(() => ({ current: end }), [end]);
    const setStart = (p: P) => {
        startRef.current[0] = p[0];
        startRef.current[1] = p[1];
        console.log('setting start', p);
        setStartX(p);
    };
    const setEnd = (p: P) => {
        endRef.current = p;
        setEndX(p);
    };
    useEffect(() => {
        if (!ds) return;
        const rd = ds.getDimension('range_dimension') as RangeDimension;
        chart.removeFilter = () => {
            //todo this is probably bad, especially in the general case - what if there's more than one filter?
            rd.removeFilter();
            setStart([0, 0]);
            setEnd([0, 0]);
        }
        setRangeDimension(rd);

        return () => {
            chart.removeFilter = () => { };
            rd.destroy();
        }
    }, [ds]);
    return { polygonLayer, rangeDimension, start, setStart, startRef, end, setEnd, endRef };
}
function useCreateMeasure() {
    const [startPixels, setStart] = useState<P>([0, 0]);
    const [endPixels, setEnd] = useState<P>([0, 0]);
    return { startPixels, setStart, endPixels, setEnd };
}
function useCreateSpatialAnnotationState(chart: BaseReactChart<any>) {
    const rectRange = useCreateRange(chart);
    const measure = useCreateMeasure();
    return { rectRange, measure };
}


export function SpatialAnnotationProvider({ chart, children }: { chart: BaseReactChart<any> } & React.PropsWithChildren) {
    const annotationState = useCreateSpatialAnnotationState(chart);
    return (
        <SpatialAnnotationState.Provider value={annotationState}>
            {children}
        </SpatialAnnotationState.Provider>
    );
}

export function useRange() {
    const range = useContext(SpatialAnnotationState).rectRange;
    if (!range) throw new Error('no range context');
    return range;
}

export function useMeasure() {
    const measure = useContext(SpatialAnnotationState).measure;
    if (!measure) throw new Error('no measure context');
    return measure;
}