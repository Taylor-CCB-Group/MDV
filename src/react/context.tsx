import { createContext, useContext } from "react";
import type { BaseReactChart } from "./components/BaseReactChart";
import type DataStore from "../datastore/DataStore";
import type { VivConfig } from "./components/avivatorish/state";

const ChartContext = createContext<BaseReactChart<any>>(undefined);
export const DataStoreContext = createContext<DataStore>(undefined);

export function ChartProvider({
    chart,
    children,
}: { chart: BaseReactChart<any> } & React.PropsWithChildren) {
    //DataStoreContext.Provider would be applied at a wider scope if we had a global root & portals.
    return (
        <ChartContext.Provider value={chart}>
            <DataStoreContext.Provider value={chart.dataStore}>
                {children}
            </DataStoreContext.Provider>
        </ChartContext.Provider>
    );
}

export function useChart() {
    const chart = useContext(ChartContext);
    if (!chart) throw new Error("no chart context");
    //todo: typing...
    return chart;
}
export function useDataStore(foreignDataStore?: DataStore) {
    const dataStore = useContext(DataStoreContext);
    if (foreignDataStore) return foreignDataStore;
    if (!dataStore) throw new Error("no data store context");
    return dataStore;
}
export function useVivConfig(): VivConfig {
    const chart = useContext(ChartContext);
    if (!chart) return { viewerStore: {}, channelsStore: {}, imageSettingsStore: {} };
    return chart.config.viv as VivConfig;
}
