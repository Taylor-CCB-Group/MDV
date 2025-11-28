import BaseChart, { type BaseConfig } from "@/charts/BaseChart";
// import { deserialiseParam, type SerialisedParams } from "@/charts/chartConfigUtils";
import { type FieldSpecs, flattenFields, type FieldSpec } from "@/lib/columnTypeHelpers";
import { isArray } from "@/lib/utils";
import { RowsAsColsQuery } from "@/links/link_utils";
import { action, type IReactionDisposer } from "mobx";

/**
 * Apply the {@link loadColumnData} decorator to a method on a chart class manually at runtime,
 * rather than using the `@loadColumnData` syntax
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
    // only used by vanilla DensityScatterPlot
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
     * prevent nested methods from adding more `userValues` & confusing things? 
     * e.g. if `setParams([someQuery])` calls `setFields(params.slice(1))` internally,
     * which one should be the active query?
     * We don't want to use the inner one, which will've been called with string column names.
     * The answer seems to be that the outermost method should always be the active query...
     * and if inner previously set a query, it should be cleared.
     * but perhaps a better answer is to avoid this kind of nesting... which is what we're also doing anyway.
     *! So this is probably not needed - and not properly tested.
     * added because of `methodName: 'replacementMethod'` issue, now resolved differently.
     */
    inSetQuery = false;
    /**
     * Called from the implementation of `@loadColumnData` decorator to ensure proper housekeeping of reactions
     * and keeping track of object form of queries for serialisation.
     * @param methodName the name of the method being called
     * @param userValue the value passed to the method before any transformation
     */
    setActiveQuery(methodName: DecoratedMethodName, userValue: UserArgs, callback: () => void) {
        const alreadyInSetQuery = this.inSetQuery;
        this.inSetQuery = true;
        if (this.reactionDisposers.has(methodName)) {
            const disposer = this.reactionDisposers.get(methodName);
            disposer?.();
            this.reactionDisposers.delete(methodName);
            this.chart.reactionDisposers = this.chart.reactionDisposers.filter(v => v !== disposer);
        }
        if (!alreadyInSetQuery) {
            // Only store if userValue contains a RowsAsColsQuery (active link)
            // it would be better to avoid chart-specific logic here 
            // - and the selection dialog, as a react chart, should be able to avoid this monkey patching
            // but the way it uses params would probably need more of a refactor to avoid this
            if (this.chart.config.type === 'selection_dialog') {
                if (Array.isArray(userValue) && userValue.some(v => v instanceof RowsAsColsQuery)) {
                    this.userValues[methodName] = userValue;
                } else {
                    delete this.userValues[methodName];
                }
            } else {
                this.userValues[methodName] = userValue;
            }
        } else {
            this.userValues[methodName] = null;
        }
        //using chart.mobxAutorun() to ensure that the reaction is disposed when the chart is disposed
        const disposer = this.chart.mobxAutorun(callback);
        this.reactionDisposers.set(methodName, disposer);
        if (!alreadyInSetQuery) this.inSetQuery = false;
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
    //- not implementing for now, 
    // concentrating on more localised logic for params to reduce testing scope
    // activeChartConfig() {
    //     const config = this.chart.getConfig();
    //     //but deserialised... we don't have a clean method for that, but we can handle param specifically
    //     const param: SerialisedParams = config.param;
    //     //!! this is fishy - we don't want new objects, we want to use the activeQueries.userValues
    //     const processed = param.map(p => deserialiseParam(this.chart.dataStore, p));
    //     config.param = processed;
    //     return config;
    // }
    /**
     * Returns the active parameters for the chart.
     * 
     * This is the same as the `config.param` property, but it will be the active state of the parameters,
     * rather than the serialised state (or strings for current fields).
     * 
     * This is useful for getting the active state of the parameters, which is needed for things like
     * the parameter setting dialog.
     * 
     * We would prefer to have a method on the chart object that returns the entire active config,
     * but keeping the implementation focused on params/setParams means that we can be more confident
     * that the logic is correct.
     */
    activeParams() {
        // if there is some value in activeQueries.userValues['setParams'] that might work...
        if ('setParams' in this.userValues) {
            return this.userValues['setParams'] as FieldSpecs;
        }
        return this.chart.config.param;
    }
}
