import { createContext, type PropsWithChildren, useContext } from "react";
import type ChartManager from "@/charts/ChartManager";

const ChartManagerContext = createContext<ChartManager | null>(null);

type ChartManagerProviderProps = PropsWithChildren<{
    chartManager?: ChartManager | null;
}>;

export function ChartManagerProvider({ chartManager, children }: ChartManagerProviderProps) {
    return <ChartManagerContext.Provider value={chartManager ?? null}>{children}</ChartManagerContext.Provider>;
}

export function useChartManagerContext() {
    return useContext(ChartManagerContext);
}
