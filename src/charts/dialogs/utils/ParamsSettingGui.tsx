import BaseChart, { type BaseConfig } from "@/charts/BaseChart";
import { type FieldSpec, type FieldSpecs, isMultiColumn } from "@/lib/columnTypeHelpers";
import { g, isArray } from "@/lib/utils";

/**
 * Helper method to go figure out which entries in a config.param array correspond
 * to a given parameter at index i in the chart type's params array.
 * There is no perfect way to do this - but we should be able to rely on knowing
 * that there is at most one multi-column parameter.
 */
function getCurrentParam<T extends BaseConfig>(chart: BaseChart<T>, i: number) {
    const { config } = chart;
    const chartType = BaseChart.types[config.type];
    const { params } = chartType;
    // we only call this in the context of a chart type that has params
    if (!params) throw new Error("No params for chart type");
    const isMultiType = isMultiColumn(params[i].type);

    // const currentParams = config.param;
    const currentParams = chart.activeQueries.userValues['setParams'] || config.param;
    const n = currentParams.length;
    const hasMulti = params.some(p => isMultiColumn(p.type));
    if (isMultiType) {
        //! as far as we know, we only ever see params ordered like [_multi] or [single, _multi]
        const nOthers = params.length - 1;
        if (i === 0) {
            // we want the head of the array... the length minus however many other params there are
            const end = currentParams.length - nOthers;
            return currentParams.slice(0, end);
        }
        // If param type is _multi but it's not the first element of params array
        // we want the tail of the array, starting from where the nOthers single params end
        return currentParams.slice(nOthers);
    }
    if (hasMulti) {
        // there is a multi-column parameter, and this is not it.
        // we already handled the case where we are the multi-column parameter
        // that means we should be either first or last in the array.
        if (i === 0) {
            return currentParams[0];
        }
        return currentParams[n - 1];
    }
    // there is no multi-column parameter; happy days
    return currentParams[i];
}

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
        const current_value = getCurrentParam(chart, i);
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