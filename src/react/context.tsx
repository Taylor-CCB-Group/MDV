import { createContext, useContext, useEffect, useMemo, useState } from "react";
import { BaseReactChart } from "./components/BaseReactChart";
// import { useLocalStore } from "mobx-react-lite";
import DataStore from "../datastore/DataStore";
import { OME_TIFF } from "./components/avivatorish/state";
import { MaterialWrapper } from "./material";
import { useOmeTiffLoader } from "./hooks"; //importing this here stops HMR working for other hooks... use smaller modules...
import RangeDimension from "../datastore/RangeDimension";

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



const ChartContext = createContext<BaseReactChart<any>>(undefined);
// Could more usefully be thought of as ScatterplotContext?
const RangeContext = createContext<RangeState>(undefined);
const DataStoreContext = createContext<DataStore>(undefined);

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
        const rd = ds.getDimension('range_dimension');
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


export function ChartProvider({ chart, children, materialui }: { chart: BaseReactChart<any>, materialui: boolean } & React.PropsWithChildren) {
    //DataStoreContext.Provider would be applied at a wider scope if we had a global root & portals.
    const rangeState = useCreateRange(chart);
    return (
    <ChartContext.Provider value={chart}>
        <DataStoreContext.Provider value={chart.dataStore}>
            <RangeContext.Provider value={rangeState}>
                {materialui && (<MaterialWrapper>
                    {children}
                </MaterialWrapper>)}
                {!materialui && children}
            </RangeContext.Provider>
        </DataStoreContext.Provider>
    </ChartContext.Provider>)
}

export function useRange() {
    const range = useContext(RangeContext);
    if (!range) throw new Error('no range context');
    return range;
}


export function useChart() {
    const chart = useContext(ChartContext);
    if (!chart) throw new Error('no chart context');
    //todo: typing...
    return chart;
}
export function useDataStore() {
    const dataStore = useContext(DataStoreContext);
    if (!dataStore) throw new Error('no data store context');
    return dataStore;
}
