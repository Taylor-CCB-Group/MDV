import BaseChart from "@/charts/BaseChart";
import type { FieldName } from "@/charts/charts";
import { flattenFields, type FieldSpec } from "@/lib/columnTypeHelpers";
import type { IReactionDisposer } from "mobx";

/**
 * Apply the {@link loadColumnData} decorator to a method on a chart class manually at runtime, rather than using the `@loadColumnData` syntax
 * which as of this writing is only working with vanilla ts code - not js or react.
 */
function decorateColumnMethod<T extends BaseChart<any>>(method: string, chart: T) {
    // could try to make `method: keyof T`, but that would require a more type juggling...
    // passing { target: chart } as context is not full information - what difference does it make?
    //@ts-expect-error
    chart[method] = loadColumnData(chart[method], { target: chart, name: method });
}

/**
 * Apply start @{link loadColumnData} decorator to all methods on a chart class that require column data to be loaded
 * and which don't specify this with the `@loadColumnData` syntax.
 */
export function decorateChartColumnMethods<T extends BaseChart<any>>(chart: T) {
    const chartType = BaseChart.types[chart.config.type];
    
    if (chart.colorByColumn) {
        decorateColumnMethod("colorByColumn", chart);
    }
    if (chart.setToolTipColumn) {
        decorateColumnMethod("setToolTipColumn", chart);
    }
    if (chart.setBackgroundFilter) { //doesn't appear in the codebase
        decorateColumnMethod(
            "setBackgroundFilter",
            chart,
        );
    }
    if (chart.changeContourParameter) {
        decorateColumnMethod(
            "changeContourParameter",
            chart,
        );
    }

    if (!chartType.methodsUsingColumns) return;
    for (const m of chartType.methodsUsingColumns) {
        decorateColumnMethod(m, chart);
    }
}


// also consider configEntriesUsingColumns...
// perhaps we can have computed getters and custom setters for the columns in the config

/**
 * Decorator for a method that requires column data to be loaded.
 * 
 * Methods decorated with this will use the first argument as a column id or array of column ids.
 * 
 * Anything specified as a virtual column that may change at runtime should cause the original method to be called again
 * as a `reaction`.
 * 
 * Considering the possibility that this might also take arguments about how this will map to the config object / params....
 * Maybe it's enough in the short term to remember which method will be called with associated column specifications.
 */
export function loadColumnData<This extends BaseChart<any>, Args extends any[], Return>(
    target: (this: This, inputCol: FieldName | FieldName[], ...args: Args) => void, // we could probably type this to have first argument specified...
    // context: ClassMethodDecoratorContext<This, (this: This, inputCol: FieldSpec, ...args: Args) => void> | { target: This },
) {
    let disposer: IReactionDisposer | null = null;
    function replacementMethod(this: This, inputCol: FieldSpec, ...args: Args) {
        // when I deserialise a chart, it should have working reactions for any live columns...
        console.log("decorated loadColumnData method called", this.config.type, inputCol, args);
        const dataSource = this.dataSource.name;
        // if we have a special value indicating live data, we can do something, perhaps with chart.mobxAutorun(), here...
        // we want to make sure that
        // - the result of the live column change will be applied ✅
        // - the value stored in the config will be the special value, rather than the result of the live column change
        // - reaction disposes when the live column is no longer needed ✅ (to be tested)
        // there could be anything happening inside the method... but perhaps if we know that the config is always a
        // proxied mobx object, we can somehow generically detect changes and re-patch???
        // perhaps safer to have a copy of the config that will be used in chart.getConfig(), which has the 'special value',
        // but then when the decorated method does things to its config, it won't be reflected in the 'special value' config that we save
        // this is really a 'config vs state' thing.
        if (inputCol == null) {
            // columns not needed
            target.call(this, inputCol, ...args);
        } else {
            // this should only happen if there are any live columns...
            // and we should dispose if replacementMethod is called again
            if (disposer) {
                disposer();
                this.reactionDisposers = this.reactionDisposers.filter(v => v !== disposer);
            }
            const multiCol = Array.isArray(inputCol);
            const colsOriginalArr = multiCol ? inputCol : [inputCol];
            // args[0] = colsOriginal;
            // how do we specify which params this will operate on?
            this.activeQueries[target.name] = [...colsOriginalArr, ...args];
            disposer = this.mobxAutorun(() => {
                const cols = flattenFields(colsOriginalArr);
                const a = multiCol ? cols : cols[0];
                // const newArgs = [a, ...args.slice(1)];
                const cm = window.mdv.chartManager;
                cm._getColumnsThen(dataSource, cols, () => {
                    target.call(this, a, ...args);
                });
            });
        }
    }

    return replacementMethod;
}