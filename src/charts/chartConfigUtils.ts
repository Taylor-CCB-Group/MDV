import { getRandomString } from "@/utilities/Utilities";
import type BaseChart from "./BaseChart";
import { action, makeAutoObservable } from "mobx";
import type { FieldSpec, FieldSpecs } from "@/lib/columnTypeHelpers";
import { RowsAsColsQuery, type RowsAsColsQuerySerialized } from "@/links/link_utils";
import type { FieldName } from "./charts";
import type DataStore from "@/datastore/DataStore";
import { isArray } from "@/lib/utils";
import type { BaseConfig } from "./BaseChart";
import type { TooltipConfig } from "@/react/scatter_state";

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
export type SerialisedParams = SerialisedColumnParam[];
export function deserialiseParam(ds: DataStore, param: SerialisedColumnParam) {
    //! should we consider making this whole thing async so that we can know that whatever is needed is loaded?
    // this may not be the best place to do that, but it's a thought...
    const result = typeof param === "string" ? param : RowsAsColsQuery.fromSerialized(ds, param);
    if (!result) {
        // this happens with unexpected array...
        throw new Error(`Failed to deserialise param: ${JSON.stringify(param)}`);
    }
    return result;
}

export function getConcreteFieldNames(fieldSpec: FieldSpec | FieldSpecs) {
    if (isArray(fieldSpec)) {
        // throw new Error("Not implemented");
        return fieldSpec.flatMap((f) => typeof f === "string" ? f : f.fields);
    }
    //! nb, previously unused, now changing the signature to return string array
    //(I've had it with checking for array or string)
    return typeof fieldSpec === "string" ? [fieldSpec] : fieldSpec.fields;
}


/**
 * We use this only as an intermediate step in AddChartDialogReact, as of this writing.
 */
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
 * not currently used
 * 
 * @returns resulting `config` object if there were no errors, otherwise an object with
 *  the potentially compromised result, and a list of err
 */
/**
 * Recursively deserializes a value, handling strings, objects that need deserialization, arrays, and nested objects.
 */
function deserialiseValueRecursive(ds: DataStore, value: any): any {
    // Handle null/undefined
    if (value == null) {
        return value;
    }
    
    // Handle primitives (string, number, boolean)
    if (typeof value !== 'object') {
        return value;
    }
    
    // Handle arrays - recursively deserialize each element
    if (isArray(value)) {
        return value.map(item => deserialiseValueRecursive(ds, item));
    }
    
    // Handle objects that should be deserialized (e.g., RowsAsColsQuery)
    if (value.type === "RowsAsColsQuery") {
        return deserialiseParam(ds, value);
    }
    
    // Handle plain objects - recursively deserialize all properties
    const result: any = {};
    for (const key in value) {
        if (Object.prototype.hasOwnProperty.call(value, key)) {
            result[key] = deserialiseValueRecursive(ds, value[key]);
        }
    }
    return result;
}

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

    // Recursively deserialize the entire config object
    config = deserialiseValueRecursive(chart.dataStore, config) as C;
    
    // Handle special properties that need additional processing beyond deserialization
    if (originalConfig.color_by) {
        const colorBy = (config as any).color_by;
        (config as any).color_by = undefined;
        chart.deferredInit(() => {
            if (!chart.colorByColumn) {
                console.error('chart does not have colorByColumn method, but had color_by in config');
                return;
            }
            // The method itself takes a string - but our decorated version takes what we're giving it...
            (chart.colorByColumn as any)?.(colorBy);
        })
    }
    
    const configWithTooltip = config as unknown as TooltipConfig;
    const tooltipColumn = configWithTooltip.tooltip?.column;
    if (tooltipColumn) {
        // Convert deserialized tooltip column to concrete field names
        if (isArray(tooltipColumn)) {
            (config as any).tooltip.column = getConcreteFieldNames(tooltipColumn);
        } else {
            const concreteFieldNames = getConcreteFieldNames(tooltipColumn);
            (config as any).tooltip.column = concreteFieldNames[0];
            chart.deferredInit(() => {
                if (chart.setToolTipColumn) chart.setToolTipColumn(tooltipColumn);
            });
        }
    }
    
    console.log(config.type, 'processed config:', config);
    
    // Handle queries - any methodsUsingColumns that don't have an associated ColumnQueryMapper methodToConfigMap
    // **we need better analysis of conditions under which this is still relevant.**
    chart.deferredInit(action(() => {
        const c = config as any;
        const queries = c.queries;
        if (!queries) return;
        
        for (const method in queries) {
            const sq = queries[method];
            if (!sq) {
                continue;
            }
            try {
                //TODO- this may not be correct wrt to adapting to/from array...
                //BUT IT SEEMS LIKE WE MIGHT ACTUALLY FINALLY BE ABLE TO DROP THIS SECTION ENTIRELY???
                console.warn(`Applying serialized query '${sq}' via a method we hope to remove...`);
                (chart as any)[method](sq);
            } catch (e) {
                console.error('failed to run query', method, config.type, sq, e);
            }
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

    //? had some weird issue with chart.deferredInit() here...
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
    });
    // note: a previous version of this used makeObservable for keeping track of onDataFiltered...
    // that worked, with extra extraneous number that changed to be observed by the hook...
    // What I have now done is change DataStore to be observable, and added a method for getting filtered indices
    // in a way that can be shared by different charts (react or otherwise).


    return config;
}

