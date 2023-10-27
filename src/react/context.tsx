import { createContext, useContext } from "react";
import { BaseReactChart } from "./components/BaseReactChart";
// import { useLocalStore } from "mobx-react-lite";
import DataStore from "../datastore/DataStore";

const ChartContext = createContext<BaseReactChart<any>>(undefined);
const DataStoreContext = createContext<DataStore>(undefined);

export function ChartProvider({ chart, children }: { chart: BaseReactChart<any>, children: any }) {
    //const mobxChart = useLocalStore(() => chart); //thought we could do without makeObserable... not helping
    //in fact, I think it leads to StackOverflows...

    //DataStoreContext.Provider would be applied at a wider scope if we had a global root & portals.
    return <ChartContext.Provider value={chart}>
        <DataStoreContext.Provider value={chart.dataStore}>
        {children}
        </DataStoreContext.Provider>
    </ChartContext.Provider>
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
