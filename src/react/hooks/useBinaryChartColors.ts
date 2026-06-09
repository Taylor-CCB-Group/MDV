import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import { useRef, useMemo } from "react";
import { useChartColors } from "./useResolvedChartColor";

export function useBinaryChartColors<T extends BaseConfig>(
    chart: BaseChart<T>
): [Uint8Array, number]{
    const colorFunction = useChartColors(chart);
    const stamp  = useRef<number>(1);

    const colorBuffer= useMemo(()=>new Uint8Array(chart.dataStore.size*3), []
    )

    return useMemo(() => {
        if (colorFunction){

            const size = chart.dataStore.size;
        
            for (let i = 0; i < size; i++) {
                const color = colorFunction(i);
                const idx = i * 3;
                colorBuffer[idx] = color[0];
                colorBuffer[idx + 1] = color[1];
                colorBuffer[idx + 2] = color[2];
            }
            stamp.current++;
        
        }
        return [colorBuffer, stamp.current];
    }, [colorFunction]);
}
