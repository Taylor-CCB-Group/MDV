import { createContext, useContext } from "react";
import { BaseReactChart } from "./components/BaseReactChart";

export const ChartContext = createContext<BaseReactChart<any>>(undefined);

export function ChartProvider({ chart, children }: { chart: BaseReactChart<any>, children: any }) {
    return <ChartContext.Provider value={chart}>
        {children}
    </ChartContext.Provider>
}

export function useChart() {
    const chart = useContext(ChartContext);
    if (!chart) throw new Error('no chart context');
    //todo: typing...
    return chart;
}