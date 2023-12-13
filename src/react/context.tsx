import { createContext, useContext } from "react";
import { BaseReactChart } from "./components/BaseReactChart";
// import { useLocalStore } from "mobx-react-lite";
import DataStore from "../datastore/DataStore";
import { OME_TIFF } from "./viv_state";
import { MaterialWrapper } from "./material";
import { useOmeTiffLoader } from "./hooks"; //importing this here stops HMR working for other hooks... use smaller modules...

const ChartContext = createContext<BaseReactChart<any>>(undefined);
const DataStoreContext = createContext<DataStore>(undefined);
const OmeTiffContext = createContext<OME_TIFF | undefined>(undefined);

export function ChartProvider({ chart, children }: { chart: BaseReactChart<any>, children: any }) {
    //const mobxChart = useLocalStore(() => chart); //thought we could do without makeObserable... not helping
    //in fact, I think it leads to StackOverflows...

    //DataStoreContext.Provider would be applied at a wider scope if we had a global root & portals.

    return (<ChartContext.Provider value={chart}>
        <DataStoreContext.Provider value={chart.dataStore}>
            <MaterialWrapper>
                {children}
            </MaterialWrapper>
        </DataStoreContext.Provider>
    </ChartContext.Provider>)
}


export function OmeTiffProvider({ children }) {
    //OmeTiffContext.Provider is not always going to be there for all charts... and we might want to
    //have charts that use multiple OME_TIFFs, different formats, etc.
    const ome = useOmeTiffLoader();
    return (
        <OmeTiffContext.Provider value={ome}>
            {children}
        </OmeTiffContext.Provider>
    )
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

export function useOmeTiff() {
    const ome = useContext(OmeTiffContext);
    return ome;
}
