import BaseChart, { type BaseConfig } from "@/charts/BaseChart";
import { flattenFields, type FieldSpec } from "@/lib/columnTypeHelpers";
import { isArray } from "@/lib/utils";
import { RowsAsColsQuery } from "@/links/link_utils";
import { action, type IReactionDisposer } from "mobx";

/**
 * Apply the {@link loadColumnData} decorator to a method on a chart class manually at runtime, rather than using the `@loadColumnData` syntax
 * which as of this writing is only working with vanilla ts code - not js or react.
 */
function decorateColumnMethod<T extends BaseChart<any>>(method: string, chart: T) {
    // could try to make `method: keyof T`, but that would require a more type juggling...
    // passing { target: chart } as context is not full information - what difference does it make?
    //@ts-ignore metaprogramming - really not worth the effort to get this to work with typescript
    if (chart[method].name === "replacementMethod") {
        console.warn("Method already decorated", method);
        return;
    }
    //@ts-ignore metaprogramming - really not worth the effort to get this to work with typescript
    chart[method] = loadColumnData(chart[method]); //, { target: chart, name: method });
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
    console.warn("Decorating chart column methods - todo, test with @loadColumnData", chartType.methodsUsingColumns);
    for (const m of chartType.methodsUsingColumns) {
        decorateColumnMethod(m, chart);
    }
}


// also consider configEntriesUsingColumns...
// perhaps we can have computed getters and custom setters for the columns in the config

type DecoratorTarget<This extends BaseChart<any>, Args extends any[], Return> = (this: This, inputCol: FieldSpec, ...args: Args) => Return;

/**
 * Decorator for a method that requires column data to be loaded.
 * 
 * Methods decorated with this will use the first argument as a column id or array of column ids.
 * 
 * 
 * Anything specified as a virtual column that may change at runtime should cause the original method to be called again
 * as a `reaction`.
 * 
 * Considering the possibility that this might also take arguments about how this will map to the config object / params...
 * We keep track of active queries with a `ColumnQueryMapper` object. An optional argument to this decorator
 * seems neater than having to set `activeQueries.methodToConfigMap` from elsewhere, but
 * not worth the change required here as of now.
 */
export function loadColumnData<This extends BaseChart<any>, Args extends any[], Return>(
    target: DecoratorTarget<This, Args, Return>, // we could probably type this to have first argument specified...
    // context: ClassMethodDecoratorContext<This, (this: This, inputCol: FieldSpec, ...args: Args) => void> | { target: This },
) {
    function replacementMethod(this: This, inputCol: FieldSpec, ...args: Args) {
        // when I deserialise a chart, it should have working reactions for any live columns...
        console.log("decorated loadColumnData method called", this.config.type, inputCol, args);
        const dataSource = this.dataSource.name;
        // if we have a special value indicating live data, we can do something, perhaps with chart.mobxAutorun(), here...
        // we want to make sure that
        // - the result of the live column change will be applied ✅
        // - the value stored in the config will be the special value, rather than the result of the live column change
        //   it'll be mapped back during serialisation as long as there's some `methodToConfigMap` appropriately configured ✅
        // - old reaction disposes when user changes the value ✅
        // - multiple instances of the same class don't interfere with each other ✅
        //   (no longer keeping a reference to disposer in the outer scope)
        
        // there could be anything happening inside the method... but perhaps if we know that the config is always a
        // proxied mobx object, we can somehow generically detect changes and re-patch???
        // perhaps safer to have a copy of the config that will be used in chart.getConfig(), which has the 'special value',
        // but then when the decorated method does things to its config, it won't be reflected in the 'special value' config that we save
        // this is really a 'config vs state' thing.
        if (inputCol == null) {
            // columns not needed
            // should still dispose any reaction that was set up
            // ? (when) is this called? consider this untested
            this.activeQueries.setActiveQuery(target.name, [null, ...args], () => {});
            target.call(this, inputCol, ...args);
        } else {
            // this should only happen if there are any live columns...
            const multiCol = isArray(inputCol);
            const colsOriginalArr = multiCol ? inputCol : [inputCol];
            // changing this mechanism to use ColumnQueryMapper rather than managing disposer here
            const userValue = [...colsOriginalArr, ...args];
            this.activeQueries.setActiveQuery(target.name, userValue, () => {
                const cols = flattenFields(colsOriginalArr);
                const a = multiCol ? cols : cols[0];
                const cm = window.mdv.chartManager;
                cm._getColumnsThen(dataSource, cols, action(() => {
                    try {
                        //@ts-ignore array vs single column
                        target.call(this, a, ...args);
                    } catch (e) {
                        console.error("Error in loadColumnData method", e);
                    }
                }));
            });
        }
    }

    return replacementMethod;
}

type DecoratedMethodName = string;
type UserArgs = any[] | null;
/**
 * This class is used to map column queries to `@loadColumnData` methods on a chart object.
 * 
 * It will manage the reactions that will update the chart when the column query changes,
 * ensuring that they are disposed of appropriately, as well as tracking the active state
 * for serialisation.
 * 
 * As such, charts internally using `string` column names continue to internally use this format,
 * with values representing the concrete manifestation of any live queries at a given moment.
 * This `ColumnQueryMapper` will keep track of how these values relate to whatever configuration
 * a user may have set.
 */
export class ColumnQueryMapper<T extends BaseConfig> {
    reactionDisposers: Map<DecoratedMethodName, IReactionDisposer> = new Map();
    // using a `Record` so it's similar to previous `activeQueries` property
    userValues: Record<DecoratedMethodName, UserArgs> = {};
    /**
     * A map from method names to the path in the chart config object where the column query should be stored.
     * Note - a single element array is used to indicate that the value should be stored as an array.
     * It does not relate to the depth of the object in the config.
     * 
     * So `'tooltip.column'` means `config.tooltip.column` taking a non-array value,
     * while `['param']` means `config.param` taking an array value.
     */
    methodToConfigMap: Record<DecoratedMethodName, string | [string]>;
    chart: BaseChart<T>;
    constructor(chart: BaseChart<T>, methodToConfigMap: Record<DecoratedMethodName, string | [string]>) {
        this.methodToConfigMap = methodToConfigMap;
        this.chart = chart;
    }
    /**
     * Called from the implementation of `@loadColumnData` decorator to ensure proper housekeeping of reactions
     * and keeping track of object form of queries for serialisation.
     * @param methodName the name of the method being called
     * @param userValue the value passed to the method before any transformation
     */
    setActiveQuery(methodName: DecoratedMethodName, userValue: UserArgs, callback: () => void) {
        if (this.reactionDisposers.has(methodName)) {
            const disposer = this.reactionDisposers.get(methodName);
            disposer?.();
            this.reactionDisposers.delete(methodName);
            this.chart.reactionDisposers = this.chart.reactionDisposers.filter(v => v !== disposer);
        }
        this.userValues[methodName] = userValue;
        //using chart.mobxAutorun() to ensure that the reaction is disposed when the chart is disposed
        const disposer = this.chart.mobxAutorun(callback);
        this.reactionDisposers.set(methodName, disposer);
    }
    /**
     * This method is called at the start of the chart serialisation process to ensure that the chart configuration
     * has the correct values for any live queries, 
     * prior to possible further processing by overrides of `BaseChart.getConfig()`.
     * 
     * Returns something like a `SerialisedConfig<T>`, except that we don't really have a type for that yet.
     */
    serialiseChartConfig() {
        const serialized = JSON.parse(JSON.stringify(this.chart.config));
        const activeQueries = this.userValues;
        // serialized.queries = {};

        for (const k in activeQueries) {
            const methodName = k;
            const configKeyA = this.methodToConfigMap[methodName];
            // if we have something like `['param']`, we know that the value should be stored as an array
            const serializedValueA = activeQueries[k]?.map(q => {
                if (q instanceof RowsAsColsQuery) {
                    return q.toJSON();
                }
                return q;
            });
            if (!serializedValueA) continue;
            const array = isArray(configKeyA);
            const serializedValue = array ? serializedValueA : serializedValueA[0];
            const configKey = array ? configKeyA[0] : configKeyA;
            if (!configKey) {
                if (!serialized.queries) serialized.queries = {};
                serialized.queries[methodName] = serializedValue;
                continue;
            }
            const configPath = configKey.split(".");
            const configObj = serialized;
            // recursively search for the property in the config object and set the value
            function setConfigValue(obj: any, path: string[], value: any) {
                if (path.length === 1) {
                    obj[path[0]] = value;
                    return;
                }
                const nextObj = obj[path[0]];
                if (nextObj == null) {
                    obj[path[0]] = {};
                }
                setConfigValue(obj[path[0]], path.slice(1), value);
            }
            setConfigValue(configObj, configPath, serializedValue);
        }

        return serialized;
    }
}
