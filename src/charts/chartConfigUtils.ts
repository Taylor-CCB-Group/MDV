import { getRandomString } from "@/utilities/Utilities";
import type BaseChart from "./BaseChart";
import { action, type IReactionDisposer, makeAutoObservable } from "mobx";
import type { FieldSpec } from "@/lib/columnTypeHelpers";
import { RowsAsColsQuery, type RowsAsColsQuerySerialized } from "@/links/link_utils";
import type { DataSource, FieldName } from "./charts";
import type DataStore from "@/datastore/DataStore";
import { isArray } from "@/lib/utils";
import type { BaseConfig } from "./BaseChart";

/**
 * This is a utility module for handling the serialisation and deserialisation of chart configurations.
 * 
 * In particular, with the introduction of dynamic column queries, we need to be able to represent these.
 * There are essentially three lifecycle stages for something that specifies a column in a chart configuration:
 * - A `string` representing a concrete column name (which may have special formatting to indicate it's a linked column)
 * - A live object representing a query that will return observable column (e.g. RowsAsColsQuery)
 * - A serialised object representing the query that can be stored in a JSON object.
 * 
 * A lot of existing charts expect strings, and we have mechanisms for instrumenting these to accept a query object - 
 * but this requires manual adaptation of the code in a way that can be error-prone.
 * 
 * It is hoped that in most cases newer charts based on React will be able to handle the query objects directly,
 * as long as it is clear at what point the configuration object is serialised and deserialised, they should be able
 * to use simple hooks, which will transparently return appropriate column objects (as in, with actual column data)
 * without being too concerned about the serialisation format etc.
 * 
 * ! as well as RowsAsColsQuery, I may implement some synthetic mock column types that can be used for testing.
 * ^^ perhaps we could try to design the error-reporting mechanism to describe something that we could reproduce without the original data???
 */


// const ParamSpec = {
//     "linkedDsName": "genes",
//         "maxItems": 10,
//             "type": "RowsAsColsQuery"
// }
// consider having synthetic columns that can be used for testing...
// possibly including things like multi-text columns that automatically mutate over time...
// - maybe serialised if it's loaded from a file, but internally do we allow it to be a live RowAsColsQuery?
type SerialisedColumnParam = (FieldName | RowsAsColsQuerySerialized);
type SerialisedParams = SerialisedColumnParam[];
export function deserialiseParam(ds: DataStore, param: SerialisedColumnParam) {
    const result = typeof param === "string" ? param : RowsAsColsQuery.fromSerialized(ds, param);
    if (!result) {
        // this happens with unexpected array...
        throw new Error(`Failed to deserialise param: ${JSON.stringify(param)}`);
    }
    return result;
}

export function getConcreteFieldNames(fieldSpec: FieldSpec) {
    if (isArray(fieldSpec)) {
        throw new Error("Not implemented");
    }
    //! nb, previously unused, now changing the signature to return string array
    //(I've had it with checking for array or string)
    return typeof fieldSpec === "string" ? [fieldSpec] : fieldSpec.fields;
}

/**
 * For charts that aren't able to accept observable props for columns,
 * this function will replace any query objects with the corresponding field name(s).
 * 
 * It takes as input a config object that may contain serialised queries, and returns a new object
 * with the queries replaced by the field names as of the moment at which the method is run.
 * 
 * @param config object in serialised form
 */
// export function getConcreteConfig(config: any) {
//     // return concrete;
// }


export function serialiseQueries(chart: BaseChart<any>) {
    const activeQueries = chart.activeQueries.userValues;
    const serialized: Record<string, any> = {};
    for (const k in activeQueries) {
        serialized[k] = activeQueries[k]?.map(q => {
            if (q instanceof RowsAsColsQuery) {
                return q.toJSON();
            }
            return q;
        });
    }
    return serialized;
}

export function serialiseConfig(config: any) {
    // we should find any mapped column queries and replace relevant values with a representation of that
    //the idea is that anywhere we previously had a column name, we can have a query object
    // Is this enough?
    // (as long as we don't have some vanilla chart with column query nonsense in it...)
    // n.b. now using `toJSON()` in favour of `toString()` where possible - but there isn't a default `toJSON()` 
    // for arbitrary objects, so still using stringification here.
    const serialized = JSON.parse(JSON.stringify(config));
    
    //! pending more generic approach to serialising queries...
    //in principle, react-based charts shouldn't have any trouble with the above...
    // if (config.color_by) {
    //     serialized.color_by = getConcreteFieldName(config.color_by);
    // }
    console.log('processed config:', serialized);
    return serialized;
}
/**
 * 
 * ! this may be precisely the kind of place where zod would be appropriate.
 * 
 * @returns resulting `config` object if there were no errors, otherwise an object with
 *  the potentially compromised result, and a list of err
 */
export function deserialiseConfig(ds: DataStore, serialConfig: any) {
    // we need to know which DataSource this is associated with to be able to deserialize
    const exceptions: { error: Error | unknown, key: string, value: any }[] = [];
    const config = JSON.parse(serialConfig, (key, value) => {
        // not sure we can have a general purpose Serializable interface with a static factory method
        // but there should be something simple along those lines we could consider...
        // ... for now, we are looking specifically for RowsAsColumnsQuerySerialized
        if (value.type === "RowsAsColsQuery") {
            try {
                return deserialiseParam(ds, value)
            } catch (error) {
                exceptions.push({
                    error,
                    key,
                    value
                });
                return value;
            }
        }
    });
    return exceptions.length ? { config, exceptions } : config;
}

export function serialiseChart<T extends BaseChart<any>>(chart: T) {
    const { config } = chart;
    // thinking about config vs state, in the context of dynamic virtual columns...
    // if we just have a record of paramSpecs, then we can use those in this representation
    // while the config object itself will have the actual evaluated runtime state of the values
    // then if we also keep a record of other `configEntriesUsingColumns`, that should, hypothetically,
    // be enough information?
    // - wouldn't like to count on that.
    // If we make sure that any special config values are of some type other than string, 
    // then we should be easily able to traverse the config object for them... however, it works the other way around:
    // the config object is liable to have computed strings & we want to get the 'special value' from that.

    // get the BaseChart.types entry for this chart, 
    // use it to determine any `configEntriesUsingColumns`...
    const serialized = JSON.parse(JSON.stringify(config));
    
    // we should find any mapped column queries and replace relevant values with a representation of that
    //the idea is that anywhere we previously had a column name, we can have a query object
    const { colorByColumn, ...serialisedQueries } = serialiseQueries(chart);
    serialized.queries = serialisedQueries;
    if (colorByColumn) {
        serialized.color_by = colorByColumn[0];
    }
    console.log('processed config:', serialized);
    return serialized;
}

/**
 * This will be called by the chart constructor to set up the config object, as well as properties on the chart object
 * for observing changes to the config object.
 * 
 * Implementation notes in progress:
 * - should we have functions that only operate on config objects?
 *   (in future we might not even need Chart objects, just config objects and components that can render them)
 * - is the config.set property safe?
 * - we should definitely have clearer expression of the notion of SerialisedConfig vs RuntimeConfig types.
 * 
 * @param originalConfig
 * @param chart
 */
export function initialiseChartConfig<C extends BaseConfig, T extends BaseChart<C>>(originalConfig: C, chart: T) {
    let config: C = JSON.parse(JSON.stringify(originalConfig));
    if (!config.id) {
        // what about when we duplicate a chart?
        config.id = getRandomString();
    }

    //we might introspect this to figure out which config entries are known to use columns...
    //const chartTypeInfo = BaseChart.types[chart.config.type];

    //or rather than relying on that property we can traverse the config object for any special values
    //for now, only operating on the param property of the config object
    //also we move around where actual loading of column data currently done by ChartManager happens
    //todo process entire config object, not just param
    //@ts-expect-error todo distinguish type of serialised vs runtime config
    const param: SerialisedParams = config.param;
    const processed = param.map(p => deserialiseParam(chart.dataStore, p));
    config.param = processed;
    //pending more generic approach to serialising queries...
    if (originalConfig.color_by) {
        //@ts-expect-error color_by
        const colorBy = isArray(originalConfig.color_by) ? deserialiseParam(chart.dataStore, originalConfig.color_by[0]) : config.color_by = deserialiseParam(chart.dataStore, originalConfig.color_by);
        config.color_by = undefined;
        setTimeout(() => {
            if (!chart.colorByColumn) {
                console.error('chart does not have colorByColumn method, but had color_by in config');
                return;
            }
            //@ts-expect-error color_by the method itself takes a string - but our decorated version takes what we're giving it...
            chart.colorByColumn?.(colorBy);
        })
    }
    console.log(config.type, 'processed config:', config);
    //temporary way of prototyping query
    //may return to something like this for non-react charts as it requires less boilerplate & chart specific code
    setTimeout(action(() => {
        const c = config as any;
        const queries = c.queries;
        for (const method in queries) {
            // if we want to handle multiple arguments that had been provided to some function, this would be the place
            const q = queries[method].map((v: any) => deserialiseParam(chart.dataStore, v));
            (chart as any)[method](q);
        }
    }));

    // makeAutoObservable(config);
    Object.defineProperty(chart, "config", {
        get: () => config,
        set: (v) => {
            // this will be invoked after initialiseConfig returns in BaseChart constructure,
            // so we no long makeAutoObservable(config) above, but this will apply it there
            config = v; // re-assigning config ref in this closure is not really relevant now;
            //we are using makeAutoObservable and not referring to config directly...
            // - this should be a good place to call a method that will update the chart based on the new config
            // current BaseChart.setParams() should probably be replaced.
            //! - maybe better to have separate methods given that we start need logic like update legend based on color/params
            // this should only happen with config objects that haven't been observed yet...
            makeAutoObservable(config);
        },
    });

    setTimeout(() => {
        // defer this until after the constructor has finished
        chart.mobxAutorun(() => {
            // setTitle also sets config.title - which is what we're observing here...
            // so we got warnings about setTitle mutating config.title 'outside an action' (at least it didn't go into infinite loop).
            // perhaps title.textContent could be a computed value and we may not need this autorun at all...
            // this.title.textContent = config.title;
            // we can now safely call this.setTitle() without warnings as it avoids unnecessary config.title changes.
            chart.setTitle(config.title);
        });
    }, 0);
    // note: a previous version of this used makeObservable for keeping track of onDataFiltered...
    // that worked, with extra extraneous number that changed to be observed by the hook...
    // What I have now done is change DataStore to be observable, and added a method for getting filtered indices
    // in a way that can be shared by different charts (react or otherwise).


    return config;
}

