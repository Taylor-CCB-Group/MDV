import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import { useFieldSpec } from "@/react/hooks";
import { useEffect, useMemo } from "react";



export function useChartColors<T extends BaseConfig>(chart: BaseChart<T>) {
    const colorField = String(chart.config.color_by);
    const colorColumn = useFieldSpec(colorField);
    const logColorScale = chart.config.log_color_scale;
    const trimColorScale = chart.config.trim_color_scale;
    const fallbackOnZero = chart.config.fallbackOnZero;
    const hideMissing = chart.config.hideMissing;
  

    useEffect(() => {
        if (!colorField || !colorColumn || colorColumn.field !== colorField) {
            return;
        }
        if (!chart.config.color_legend) {
            chart.config.color_legend = { display: true };
        }
        const timeoutId = window.setTimeout(() => {
            chart.setColorLegend();
        }, 0);
        return () => {
            window.clearTimeout(timeoutId);
        };
    }, [chart, colorColumn, colorField,logColorScale, trimColorScale, fallbackOnZero, hideMissing]);

    const colorFunction = useMemo(() => {
        if (!colorField || !colorColumn || colorColumn.field !== colorField) {
            return null;
        }
        const conf = {
            asArray: true,
            overideValues: {
                colorLogScale: logColorScale,
                fallbackOnZero,
                hideMissing,
            },
        };
        chart._addTrimmedColor(colorField, conf);
        return chart.dataStore.getColorFunction(colorField, conf);
    }, [
        chart,
        colorColumn,
        colorField,
        fallbackOnZero,
        hideMissing,
        logColorScale,
        trimColorScale,
    ]);

    return colorFunction;
}
