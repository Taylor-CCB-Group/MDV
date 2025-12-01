import { computed, makeObservable, observable, runInAction } from "mobx";
import type ChartManager from "../charts/ChartManager";
import type { ColumnName, DataColumn, DataType, FieldName } from "../charts/charts";
import type DataStore from "@/datastore/DataStore";
import type { loadColumnData } from "../datastore/decorateColumnMethod";
import { isColumnLoaded } from "@/lib/columnTypeHelpers";
// zod?

export type ChartLink = {
    type: "chart_columnval_link";
    id: string;
    source_chart: string;
    column: ColumnName;
    target_charts: string[];
    set_color: boolean;
    set_title: boolean;
    set_tooltip: boolean;
    set_legend: boolean;
    /** if true, the parm at given parm_index will be set...*/
    set_param: boolean;
    param_index?: number;
};

export type LinkTarget = {
    target_chart: string;
    //maybe we might want to allow for a different column in the source chart to be used to set the target chart's properties...
    //not for now though...
    // source_column: string,
    set_color: boolean;
    set_title: boolean;
    set_tooltip: boolean;
    set_legend: boolean;
    set_param: boolean;
    param_index?: number;
};

export type HighlightColumnLink = {
    type: "highlight_column_link";
    id: string;
    source_ds: string;
    source_column: ColumnName;
    targets: LinkTarget[];
};

export function addHighlightColumnLink(
    link: HighlightColumnLink,
    cm: ChartManager,
) {
    /// this implementation filled in by copilot & not yet read / used...
    const ds = cm.getDataSource(link.source_ds);
    if (!ds) {
        console.error(`DataStore ${link.source_ds} not found`);
        return;
    }
    const srcCol = ds.columnIndex[link.source_column];
    if (!srcCol) {
        console.error(
            `Column ${link.source_column} not found in DataStore ${link.source_ds}`,
        );
        return;
    }
    if (!isColumnLoaded(srcCol)) {
        console.error(`${srcCol} not loaded`);
        return;
    }

    async function updateValue(newValue: string) {
        for (const target of link.targets) {
            const targetChart = cm.getChart(target.target_chart);
            if (!targetChart) {
                console.error(`Chart ${target.target_chart} not found`);
                continue;
            }
            if (target.set_color) targetChart.colorByColumn(newValue);
            if (target.set_title) targetChart.setTitle(newValue);
            if (target.set_tooltip) {
                if (targetChart.setToolTipColumn)
                    targetChart.setToolTipColumn(newValue);
                else if (targetChart.config.tooltip) {
                    runInAction(
                        () => (targetChart.config.tooltip.column = newValue),
                    );
                }
            }
            if (target.set_param) {
                runInAction(() => {
                    if (target.param_index === undefined) return;
                    targetChart.config.param[target.param_index] = newValue;
                });
            }
        }
    }

    ds.addListener(link.id, async (type: string, data: any) => {
        if (type === "data_highlighted") {
            const newValue = srcCol.values[srcCol.data[data.indexes[0]]];
            updateValue(newValue);
        }
    });
}

/**
 * Add a link between two charts. When the source chart is clicked, the target charts will be updated based on the value of the source chart's column.
 * The column should be a categorical/text column in the source_chart, the value of the highlighted item to be used to set properties on target_charts.
 */
export function addChartLink(link: ChartLink, cm: ChartManager) {
    const chart = cm.charts[link.source_chart];
    if (!chart) {
        console.error(`Chart ${link.source_chart} not found`);
        return;
    }
    const ds = chart.dataSource.dataStore;
    const srcCol = ds.columnIndex[link.column];
    if (!srcCol) {
        console.error(
            `"chart_columnval_link" link column "${link.column}" not found in dataStore "${ds.name}"`,
        );
        return;
    }
    if (!isColumnLoaded(srcCol)) {
        console.error(`${srcCol} not loaded`);
        return;
    }

    async function updateValue(newValue: string) {
        await cm._getColumnsAsync(ds.name, [link.column]); //may not be necessary as `colorByColumn` is 'decorated' to load column
        // ^^ also, what about the column for the destination chart/ds?
        for (const cid of link.target_charts) {
            const target = cm.getChart(cid);
            if (link.set_color) target.colorByColumn(newValue);
            if (link.set_title) target.setTitle(newValue);
            // if (link.set_legend) chart.setLegend(newValue); //todo
            if (link.set_tooltip) {
                if (target.setToolTipColumn)
                    target.setToolTipColumn(link.column);
                else if (target.config.tooltip) {
                    runInAction(
                        () => (target.config.tooltip.column = newValue),
                    );
                }
            }
            if (link.set_param) {
                runInAction(() => {
                    if (link.param_index === undefined) return;
                    target.config.param[link.param_index] = newValue;
                    //not sure what the most appropriate way to update non-react charts is... needs some work
                    if (target.config.type === "table_chart") {
                        // good luck setting grid columns etc...
                        // const cols = target.grid.getColumns();
                        // const oldCol = cols[link.param_index];
                        // cols.splice[link.param_index] = { name: newValue, field: newValue, id: newValue, sortable: true, ...oldCol };
                        // target.grid.setColumns(cols);
                        // target.grid.invalidate();
                        // // target.dataModel.setColumns(target.config.param);
                        // target.dataModel.createColumn(newValue);
                        // target.dataModel.removeColumn(oldCol.field, true, true);
                    }
                });
            }
        }
    }

    chart.chart.addListener(link.id, async (type: string, data: any) => {
        if (type === "cell_clicked") {
            // this is specific to HeatMap cells...
            updateValue(data.row);
        }
        if (type === "data_highlighted" && srcCol) {
            const newValue = srcCol.values[srcCol.data[data.indexes[0]]];
            updateValue(newValue);
        }
    });
}

interface IRowAsColumn {
    /** Index of the row in the target data-store, from which other properties can be derived */
    index: number;
    /** Human-readable form of column name */
    name: ColumnName;
    /**
     * Unique identifier for the column, in the form of `"subgroup|name (subgroup)|index"`,
     * used internally for querying the associated column in the parent {@link DataStore}.
     */
    fieldName: FieldName;
    /**
     * The actual {@link DataColumn} instance, which will be added to the parent {@link DataStore} when accessed.
     * It should be possible for consumers to access this property directly in order to have an object that can be used
     * for accessing the column (meta)data.
     */
    column: DataColumn<DataType>;
}

export type RowsAsColslink = {
    name_column: string;
    name: string;
    subgroups: {
        [sgName: string]: {
            name: string;
            label: string;
            type: string;
        }
    },
    /** 
     * added at runtime, observable array corresponding to currently selected/highlighted items.
     * In order to support multiple subgroups, should this be Map<string, IRowAsColumn[]> (or Record, not important)
     */
    observableFields: IRowAsColumn[];
    //^^ what if it was iterator instead of Array? and/or `column` can be a computed value?
    // observableColumns: DataColumn<DataType>[]; //computed?
    /** 
     * Added at runtime, given a value appearing in `name_column`, the index of a row corresponding to that value.
     * This can be used in the formation of a {@link FieldName}
     * ! note, we should be able to use 'unique' values in the name_column, and would like to do away with this
     * ^^ this is essentially a bug
     * It is now optional to reflect that it will not be available until `initPromise` resolves.
     */
    valueToRowIndex?: Map<string, number>;
    initPromise: Promise<void>;
}
/**
 * Represents an active query for a set of columns based on the currently highlighted or filtered rows.
 * ? in cases where only a single (virtual) column is needed, we can still use this interface 
 * but only access the first element of the `columns` array?
 * probably better to have a distinct interface for that.
 */
export interface MultiColumnQuery {
    columns: DataColumn<DataType>[]; //could have a generic type for this?
    fields: FieldName[];
    /**
     * When the query is first created, it may need to perform some async initialization:
     * in the case of a `RowsAsColsQuery`, it needs to know the values of the `name_column` 
     * in the linked dataSource.
     * 
     * Implementation note: as far as this interface is concerned, maybe this
     * could be optional - we could have a conformant `MultiColumnQuery`
     * that doesn't need any async initialization. For example, if it was for
     * generating mock data.
     */
    initialize(): Promise<void>;
}

// unused
// export interface ColumnQuery {
//     column: DataColumn<DataType>;
//     field: FieldName;
// }

type RowsAsColsPrefs = {
    //todo: properties for setting the link to pause / choose between highlighted or filtered data.
    //if paused, that implies we might remember the concrete values of the columns.
    //could maybe have a firstIndex as well as maxItems, for pagination
    maxItems: number;
}
export type RowsAsColsQuerySerialized = {
    linkedDsName: string;
    type: "RowsAsColsQuery";
} & RowsAsColsPrefs;

class DeserializedQueryError extends Error {
    constructor(
        // maybe something more generic or tailored in future
        public readonly serialized: RowsAsColsQuerySerialized,
        public readonly dataStore: DataStore,
        message: string
    ) {
        super(message);
        this.name = 'DeserializedQueryError';
    }

    // Helper to suggest fixes in UI (not really useful but something to maybe build on later)
    getSuggestions() {
        return {
            missingLink: this.serialized.linkedDsName,
            availableAlternatives: getRowsAsColumnsLinks(this.dataStore),
            // Could even include Levenshtein distance suggestions
        };
    }
}

/**
 * Represents an active query for a set of columns based on the currently highlighted or filtered rows
 * in a linked {@link DataSource}. Can be used in place of a {@link string} representing a column {@link FieldName}
 * 
 * Lazily evaluates with side-effects that will add column metadata into the parent {@link DataStore}, up
 * to the number specified by `maxItems`. Loading column data is handled separately as of now - for example with
 * by a {@link loadColumnData} decorator or `methods_using_columns`.
 * 
 * @param link - the link configuration object, with observable properties
 * @param maxItems - the maximum number of columns to return
 * 
 * ! we should really have sg as well as linkedDsName
 */
export class RowsAsColsQuery implements MultiColumnQuery {
    // it may be logical to make this class be the thing that encapsulates a listener,
    // don't want to introduce extra overhead... but it's probably small.
    // perhaps if rather than just link.observableFields, we have two separate arrays (!or iterators!)
    // and the computed properties can choose between them based on config flags.
    @observable accessor maxItems: number;
    constructor(public link: RowsAsColslink, public linkedDsName: string, maxItems = 1) {
        this.maxItems = maxItems;
    }
    @computed
    get columns() {
        //how would we pause this?
        return this.link.observableFields.slice(0, this.maxItems).map(f => f.column);
    }
    @computed
    get fields() {
        return this.columns.map(c => c.field);
    }
    async initialize() {
        return this.link.initPromise;
    }
    /**
     * JSON serialisation that can be used to recreate the query object.
     */
    toJSON() {
        // I was previously naively using toString() for serialization - now I know better.
        return { linkedDsName: this.linkedDsName, maxItems: this.maxItems, type: "RowsAsColsQuery" };
    }
    /**
     * @internal This should only be called by the `deserialiseParam` factory method.
     */
    static fromSerialized(ds: DataStore, serialized: RowsAsColsQuerySerialized) {
        const link = getRowsAsColumnsLinks(ds).find(l => l?.linkedDs.name === serialized.linkedDsName)?.link;
        if (!link) {
            throw new DeserializedQueryError(serialized, ds, `Link not found for ${serialized.linkedDsName}`);
        }
        return new RowsAsColsQuery(link, serialized.linkedDsName, serialized.maxItems);
    }
}
export function getFieldName(sg: string, value: string, index: number) {
    return `${sg}|${value} (${sg})|${index}`;
}
/** 
 * a given {@param link} will currently have a single listener associated with it,
 * instantiated by this function.
 * ! design needs review to account for multiple subgroups
 */
async function initRacListener(link: RowsAsColslink, ds: DataStore, tds: DataStore) {
    if (link.initPromise !== undefined) return link.initPromise;
    link.initPromise = new Promise<void>((resolve, reject) => {
        initRacListenerImpl(link, ds, tds).then(resolve, reject);
    });
    return link.initPromise;
}
async function initRacListenerImpl(link: RowsAsColslink, ds: DataStore, tds: DataStore) {
    // didn't want to re-indent this block, so it's only valid to call once as in the above function.
    if (link.initPromise !== undefined) throw new Error("initPromise already set");
    const cm = window.mdv.chartManager;
    const nameCol = tds.columnIndex[link.name_column];
    if (!nameCol) {
        console.error(`Column ${link.name_column} not found in DataStore ${ds.name}`);
        return;
    }
    // make sure this happens before any async stuff so the check above is valid and subsequent calls return early.
    link.observableFields = []; //this will be populated by setFieldsFromFilter below
    makeObservable(link, { observableFields: true });
    await cm.loadColumnSetAsync([link.name_column], tds.name);
    if (!isColumnLoaded(nameCol)) {
        throw new Error(`Column ${link.name_column} not loaded`);
    }

    // we should also add a data structure mapping each name_column value (string) to the index of the corresponding row
    // which seems to assume a 1:1 mapping between name_column values and rows? Or not? 
    // While we're here, check for anything that may be inconsistent with our assumptions.
    // Maybe we can have a lot of rows with the same name_column value, returning equivalent rows_as_columns data?
    const valueToRowIndex = new Map<string, number>();
    //! this is O(n) when we shouldn't need it...
    nameCol.data.forEach((valueIndex, rowIndex) => {
        const value = nameCol.values[valueIndex];
        if (valueToRowIndex.has(value)) {
            console.warn(`Multiple rows with the same value '${value}' in column '${link.name_column}'`);
        }
        valueToRowIndex.set(value, rowIndex);
    });
    link.valueToRowIndex = valueToRowIndex;
    
    
    
    //! we have an assumption here that there is only one subgroup
    //(nb, I think this is the thing that there'd be a radio-button for in the UI)
    //(whereas if there were multiple linked dataSources, they'd appear as different buttons for opening the dialog)
    const sgKeys = Object.keys(link.subgroups);
    const sg = Object.keys(link.subgroups)[0];
    if (sgKeys.length > 1) {
        //! this should be fixed!
        console.warn("Multiple subgroups not supported", link);
        console.warn("Using first subgroup", sg);
        //return;
    }
    if (!sg) {
        console.error("No subgroups found for link", link);
        // return;
    }
    /** ephemeral object representing a virtual column entry in a rows_as_columns link.
     * accessing the `column` property will add the column to the parent DataStore.
     */
    class RAColumn implements IRowAsColumn {
        constructor(public index: number) {}
        @computed
        get name() {
            // if I don't have `as string` here, it's inferred as string | number
            return tds.getRowText(this.index, link.name_column) as string;
        }
        @computed
        get fieldName(): FieldName {
            return getFieldName(sg, this.name, this.index);
        }
        @computed
        get column(): DataColumn<DataType> {
            //https://mobx.js.org/computeds.html
            //"They should not have side effects or update other observables."
            //this will have side effects - I think it's ok, subsequent calls will be idempotent
            //(not sure idempotent is right, thanks copilot, but it will find the columnIndex
            //and return the same DataColumn instance each time)
            //"Avoid creating and returning new observables." - ok
            //"They should not depend on non-observable values."
            return ds.addColumnFromField(this.fieldName);
        }
    }
    const getField = (index: number) => {
        return new RAColumn(index);
    };
    async function setFieldsFromFilter() {
        //! this Array.from could be suboptimal for large numbers of indices
        //at risk of overthinking, but we could have a lazy filteredIndices iterator...
        const filteredIndices = await tds.getFilteredIndices();
        const vals = Array.from(filteredIndices).map(getField);
        runInAction(() => {
            link.observableFields = vals;
        });
    }
    //todo check link.name is a good id for the listener
    tds.addListener(link.name, async (type, data) => {
        if (type === "data_highlighted") {
            const vals = data.indexes.map(getField);
            runInAction(() => {
                link.observableFields = vals;
            });
        } else if (type === "filtered") {
            // 'data' is a Dimension in this case - so we want to zip filteredIndices with the values
            
            setFieldsFromFilter();
        }
    });
    console.log("settings initial fields for link...", link);
    runInAction(() => {
        //initial population so sync callers will have something to work with
        link.observableFields = [getField(0)];
    });
    await setFieldsFromFilter();
    console.log("link initialized", link);
}
export function getRowsAsColumnsLinks(dataStore: DataStore) {
    const dataSources = window.mdv.chartManager.dataSources;
    if (dataStore.links) {
        const result = Object.keys(dataStore.links).map((linkedDsName) => {
            if (!dataStore.links) throw "unreachable";
            const links = dataStore.links[linkedDsName];
            if (links.rows_as_columns) {
                // first pass... there can be only one or zero.
                // how often will users actually want to link multiple dataSources in this way?
                // perhaps not often - but let's handle it so we don't have to change it later or have bugs.
                // There are also multiple subgroups within a single linked dataSource - also important to support.
                // UI should be simpler for the common case with a single linked dataSource/sg.
                // Are there any crazy edge cases we should consider - like indirect links? links to self?
                // !! this should be ok now, but we should test with multiple links.
                const linkedDs = dataSources.find(
                    (ds) => ds.name === linkedDsName,
                );
                if (!linkedDs) {
                    throw new Error();
                }
                // todo make sure the link is reasonably typed
                const link = links.rows_as_columns as RowsAsColslink;
                // we could be lazier about this, especially if we have a lot of links.
                initRacListener(link, dataStore, linkedDs.dataStore);
                return { linkedDs, link };
            }
        }).filter(l => l !== undefined); //todo type-safety of filter
        return result;
    }
    return [];
}
