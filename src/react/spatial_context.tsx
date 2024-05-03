import { createContext, useContext, useEffect, useMemo, useState } from "react";
import type RangeDimension from "../datastore/RangeDimension";
import { BaseReactChart } from "./components/BaseReactChart";


/*****
 * Persisting some properties related to SelectionOverlay in "RangeProvider"... >>subject to change<<.
 * Not every type of chart will have a range dimension, and not every chart will have a selection overlay etc.
 * Needs will also get more complex, and now we have a somewhat convoluted way of doing something simple.
 * Probably going to be a zustand store in not too long.
 */

type P = [number, number];
type RefP = React.MutableRefObject<P>;
type RangeState = {
    rangeDimension: RangeDimension;
    start: P; setStart: (p: P) => void; startRef: RefP;
    end: P; setEnd: (p: P) => void; endRef: RefP;
};


// Could more usefully be thought of as SpatialContext?
const RangeContext = createContext<RangeState>(undefined);

function useCreateRange(chart: BaseReactChart<any>) {
    const ds = chart.dataStore;
    // tried simpler `rangeDimesion = useMemo(...)`, but it can lead to non-destroyed rangeDimensions with HMR.
    const [rangeDimension, setRangeDimension] = useState<RangeDimension>(undefined);
    const [start, setStartX] = useState<P>([0, 0]);
    const [end, setEndX] = useState<P>([0, 0]);
    // still not sure I want these refs
    const startRef = useMemo(() => ({ current: start }), [start]);
    const endRef = useMemo(() => ({ current: end }), [end]);
    const setStart = (p: P) => {
        startRef.current = p;
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
    return { rangeDimension, start, setStart, startRef, end, setEnd, endRef };
}


export function RangeProvider({ chart, children }: { chart: BaseReactChart<any> } & React.PropsWithChildren) {
    const rangeState = useCreateRange(chart);
    return (
        <RangeContext.Provider value={rangeState}>
            {children}
        </RangeContext.Provider>
    );
}

export function useRange() {
    const range = useContext(RangeContext);
    if (!range) throw new Error('no range context');
    return range;
}
