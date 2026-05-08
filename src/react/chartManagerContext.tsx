import type ChartManager from "@/charts/ChartManager";
import { createContext, type PropsWithChildren, useContext } from "react";

const ChartManagerContext = createContext<ChartManager | null>(null);

type ChartManagerProviderProps = PropsWithChildren<{
    chartManager?: ChartManager | null;
}>;

export function ChartManagerProvider({
    chartManager,
    children,
}: ChartManagerProviderProps) {
    return (
        <ChartManagerContext.Provider value={chartManager ?? null}>
            {children}
        </ChartManagerContext.Provider>
    );
}

export function useChartManagerContext() {
    return useContext(ChartManagerContext);
}
