import BaseChart, { type BaseConfig } from "@/charts/BaseChart";
import { type FieldSpec, type FieldSpecs, isMultiColumn } from "@/lib/columnTypeHelpers";
import { g, isArray } from "@/lib/utils";


function updateMultiParam<T extends BaseConfig>(chart: BaseChart<T>, i: number, value: FieldSpecs) {
    // const config = chart.getConfig(); //we don't want serialised config here
    // and we are currently avoiding implementing a version of getConfig() that returns the active config
    const { config } = chart;
    const chartType = BaseChart.types[config.type];
    const { params } = chartType;
    if (!params) throw new Error("No params for chart type");
    if (!isMultiColumn(params[i].type)) throw new Error("Not a multi-column parameter");
    const nParams = params.length;
    if (nParams === 1) {
        //chart.config.param = value;
        chart.setParams(value);
        return;
    }
    const currentParams = chart.activeQueries.activeParams(); //config.param;
    // we know that there is only one multi-column parameter, either at start or end of array
    const nOthers = nParams - 1;
    if (i === 0) {
        const newParams = value.concat(currentParams.slice(-nOthers));
        // chart.config.param = newParams;
        chart.setParams(newParams);
        return;
    }
    const newParams = currentParams.slice(0, nOthers).concat(value);
    // chart.config.param = newParams;
    chart.setParams(newParams);
}
function updateSingleParam<T extends BaseConfig>(chart: BaseChart<T>, i: number, value: FieldSpec) {
    // calling getConfig() here will mean that we have the config in the *serialised* form, not the active-state form.
    // perhaps if we have a `getActiveConfig()` then that might be reasonably safe.
    // Essentially, we need some way that ColumnQueryMapper can be used to get the active config.
    // we shouldn't interface with that class directly - `BaseChart` should have abstract that.
    const { config } = chart;
    const chartType = BaseChart.types[config.type];
    const { params } = chartType;
    if (!params) throw new Error("No params for chart type");
    const hasMulti = params.some(p => isMultiColumn(p.type));
    // there should be an obvious canonical way to get the active state of a chart's config
    // and use of chart.activeQueries - which is something of an internal monkeypatch - is something
    const currentParams = chart.activeQueries.activeParams(); //config.param;
    const n = currentParams.length;
    if (hasMulti) {
        const multiIndex = params.findIndex(p => isMultiColumn(p.type));
        const nBefore = multiIndex < i ? multiIndex : multiIndex - 1;
        const nAfter = n - nBefore - 1;
        const newParams = currentParams.slice(0, nBefore).concat(value).concat(currentParams.slice(n - nAfter));
        chart.setParams(newParams);
        // config.param = newParams;
    } else {
        config.param[i] = value;
        chart.setParams(config.param);
    }
}

/**
 * Generate a GUI spec folder for the parameters of a chart, based on the `params`
 * in the corresponding `BaseChart.types` entry.
 * 
 * There should be a good way to set the reaction to mutate the config object
 * in a way that will be understood by the chart (as well as triggering a redraw, and
 * also making sure that this stays in sync with other things that may cause it to change).
 */
export default function getParamsGuiSpec<T extends BaseConfig>(chart: BaseChart<T>) {
    const { config } = chart;
    const chartType = BaseChart.types[config.type];
    const { params } = chartType;
    if (!params) return null;
    const multiCount = params.reduce((acc, p) => acc + (isMultiColumn(p.type) ? 1 : 0), 0);
    if (multiCount > 1) throw new Error("More than one multi-column parameter, abandon hope");
    const currentParams = chart.activeQueries.activeParams();
    const p = params.map((param, i) => {
        const isMultiType = isMultiColumn(param.type);
        const current_value = currentParams[i];
        return g({
            type: isMultiType ? "multicolumn" : "column",
            label: param.name,
            current_value,
            columnType: param.type,
            func: v => {
                // this check doesn't give an adequate type-guard
                // if (isArray(v) !== isMultiType) throw new Error("Type mismatch");
                // also the underlying functions will mutate the config object
                // - without column data loading
                // - assuming that the chart will understand things in this format
                if (isMultiType) {
                    // if (!isArray(v)) throw new Error("Expected array");
                    //! when applying the link object value, it is not an array (although it evaluates to one)
                    // it probably should be, and then we should indeed throw an error here
                    // (there may well be a related @ts-expects-error
                    //  - need to sort out `setSelectedColumn` and associated things to be less confusing).
                    if (!isArray(v)) v = [v];
                    updateMultiParam(chart, i, v);
                } else {
                    if (isArray(v)) throw new Error("Expected single value");
                    updateSingleParam(chart, i, v);
                }
                //! maybe if we had a `setParams` method, with the appropriate decorator,
                // chart.drawChart();
            }
        });
    });
    return g({
        type: "folder",
        label: "Parameters",
        current_value: p,
    })
}