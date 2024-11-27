import { getRandomString } from "@/utilities/Utilities";
import type BaseChart from "./BaseChart";
import { action, makeAutoObservable } from "mobx";
import type { FieldSpec } from "@/lib/columnTypeHelpers";
import { RowsAsColsQuery, type RowsAsColsQuerySerialized } from "@/links/link_utils";
import type { DataSource, FieldName } from "./charts";
import type DataStore from "@/datastore/DataStore";
import { isArray } from "@/lib/utils";
import type { BaseConfig } from "./BaseChart";

// const ParamSpec = {
//     "linkedDsName": "genes",
//         "maxItems": 10,
//             "type": "RowsAsColsQuery"
// }
// consider having synthetic columns that can be used for testing...
// possibly including things like multi-text columns that automatically mutate over time...
type SerialisedColumnParam = (FieldName | RowsAsColsQuerySerialized);
type SerialisedParams = SerialisedColumnParam | SerialisedColumnParam[];
export function deserialiseParam(ds: DataStore, param: SerialisedColumnParam) {
    const result = typeof param === "string" ? param : RowsAsColsQuery.fromSerialized(ds, param);
    if (!result) {
        // this happens with unexpected array...
        throw new Error(`Failed to deserialise param: ${JSON.stringify(param)}`);
    }
    return result;
}

export function getConcreteFieldName(fieldSpec: FieldSpec) {
    if (isArray(fieldSpec)) {
        throw new Error("Not implemented");
    }
    return typeof fieldSpec === "string" ? fieldSpec : fieldSpec.fields[0];
}

export class ColumnQueryMapper<T> {
    // constructor(chart: BaseChart<T>) {

    // }
}

export function serialiseQueries(chart: BaseChart<any>) {
    const { activeQueries } = chart;
    const serialized: Record<string, any> = {};
    for (const k in activeQueries) {
        serialized[k] = activeQueries[k].map(q => {
            if (q instanceof RowsAsColsQuery) {
                return q.serialized;
            }
            return q;
        });
    }
    return serialized;
}

export function serialiseConfig<T extends BaseChart<any>>(chart: T) {
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
    for (const k in chart.activeQueries) {
        // this is a different approach where we have a 'queries' property which has a method name as a key
        // so the contents of `param` are less relevant...
        // indeed, we could call this something like `columns` and use it in place of `param`...
        // would mostly apply to non-react charts though, or I'd be more keen on this approach...
        // serialized.queries = serialiseQueries(chart);
    }
    const serialisedQueries = serialiseQueries(chart);
    const colorByColumn = serialisedQueries['colorByColumn'];
    if (colorByColumn) {
        // serialized.colorByColumn = getConcreteFieldName(chart.activeQueries['colorByColumn'][0]);
        serialized.color_by = colorByColumn[0];
    }
    console.log('processed config:', serialized);
    return serialized;
}

export function initialiseConfig<C extends BaseConfig, T extends BaseChart<C>>(originalConfig: C, chart: T) {
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
    const processed = isArray(param) ? param.map(p => deserialiseParam(chart.dataStore, p)) : deserialiseParam(chart.dataStore, param);
    config.param = processed;
    //pending more generic approach to serialising queries...
    if (originalConfig.color_by) {
        //@ts-expect-error
        const colorBy = isArray(originalConfig.color_by) ? deserialiseParam(chart.dataStore, originalConfig.color_by[0]) : config.color_by = deserialiseParam(chart.dataStore, originalConfig.color_by);
        config.color_by = undefined;
        setTimeout(() => {
            if (!chart.colorByColumn) {
                console.error('chart does not have colorByColumn method, but had color_by in config');
                return;
            }
            //@ts-expect-error the method itself takes a string - but our decorated version takes what we're giving it...
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

