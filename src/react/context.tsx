import { createContext, useContext } from "react";
//import type { BaseReactChart } from "./components/BaseReactChart";

import type DataStore from "../datastore/DataStore";
import type { VivConfig } from "./components/avivatorish/state";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";

export const ChartContext = createContext<BaseChart<any>>(null as any);
export const DataStoreContext = createContext<DataStore>(null as any);

export function ChartProvider<T extends BaseConfig>({
    chart,
    children,
}: { chart: BaseChart<T> } & React.PropsWithChildren) {
    //DataStoreContext.Provider would be applied at a wider scope if we had a global root & portals.
    return (
        <ChartContext.Provider value={chart}>
            <DataStoreContext.Provider value={chart.dataStore}>
                {children}
            </DataStoreContext.Provider>
        </ChartContext.Provider>
    );
}
// infer?
function typedChart<T extends BaseConfig>(chart: BaseChart<T>) {
    const t = chart.config.type;
    return chart;
}
export function useChart() {
    //todo: typing...
    const chart = typedChart(useContext(ChartContext));
    if (!chart) throw new Error("no chart context");
    return chart;
}
export function useDataStore(foreignDataStore?: DataStore) {
    const dataStore = useContext(DataStoreContext);
    if (foreignDataStore) return foreignDataStore;
    if (!dataStore) throw new Error("no data store context");
    return dataStore;
}
/**
 * If called from within a {@link ChartContext}, this will return the charts viv config.
 * Otherwise, it will return something that should function as an empty config object.
 */
export function useVivConfig(): VivConfig {
    const chart = useContext(ChartContext);
    if (!chart) return { viewerStore: {}, channelsStore: {}, imageSettingsStore: {} };
    return chart.config.viv as VivConfig;
}
