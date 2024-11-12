import { computed, makeAutoObservable, makeObservable, runInAction } from "mobx";
import type ChartManager from "../charts/ChartManager";
import type { ColumnName, DataColumn, DataType, FieldName } from "../charts/charts";
import type DataStore from "@/datastore/DataStore";
// zod?

type ChartLink = {
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

type LinkTarget = {
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

type HighlightColumnLink = {
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
                    targetChart.config.param[target.param_index] = newValue;
                });
            }
        }
    }

    ds.addListener(link.id, async (type, data) => {
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
        console.log(
            `"chart_columnval_link" link column "${link.column}" not found in dataStore "${ds.name}"`,
        );
        // return;
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

    chart.chart.addListener(link.id, async (type, data) => {
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
    index: number;
    name: ColumnName;
    fieldName: FieldName;
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
    /** added at runtime, observable array corresponding to currently selected/highlighted items */
    observableFields: IRowAsColumn[];
    //^^ what if it was Enumerable instead of Array? and/or `column` can be a computed value?
    // observableColumns: DataColumn<DataType>[]; //computed?
}

export class RowsAsColsQuery {
    constructor(public link: RowsAsColslink, public maxItems = 1) {}
    @computed
    get columns() {
        return this.link.observableFields.slice(0, this.maxItems).map(f => f.column);
    }
    @computed
    get fields() {
        return this.columns.map(c => c.field);
    }
}


async function initRacListener(link: RowsAsColslink, ds: DataStore, tds: DataStore) {
    if (link.observableFields !== undefined) return;
    const cm = window.mdv.chartManager;
    const nameCol = tds.columnIndex[link.name_column];
    if (!nameCol) {
        console.error(`Column ${link.name_column} not found in DataStore ${ds.name}`);
        return;
    }
    link.observableFields = []; //maybe initialize this with current values
    makeObservable(link, { observableFields: true });
    
    await cm.loadColumnSetAsync([link.name_column], tds.name);
    const sg = Object.keys(link.subgroups)[0];
    class RAColumn implements IRowAsColumn {
        index: number;
        @computed
        get name() {
            return tds.getRowText(this.index, link.name_column) as string;
        }
        @computed
        get fieldName(): FieldName {
            return `${sg}|${this.name} (${sg})|${this.index}`;
        }
        @computed
        get column(): DataColumn<DataType> {
            return ds.addColumnFromField(this.fieldName);
        }
        constructor(index: number) {
            this.index = index;
        }
    }
    const getField = (index: number) => {
        // if I don't have `as string` here, it's inferred as string | number, incompatible with the type of 'value'
        // const name = tds.getRowText(index, link.name_column) as string;
        // const fieldName = `${sg}|${name} (${sg})|${index}`;
        // const column = ds.addColumnFromField(fieldName);
        // return ({ index, value: name, fieldName, column });
        return new RAColumn(index);
    };
    //todo check link.name is a good id for the listener
    tds.addListener(link.name, async (type, data) => {
        if (type === "data_highlighted") {
            const vals = data.indexes.map(getField);
            runInAction(() => {
                link.observableFields = vals;
            });
        } else if (type === "filtered") {
            // 'data' is a Dimension in this case - so we want to zip filteredIndices with the values
            const filteredIndices = await tds.getFilteredIndices();
            //! this Array.from could be suboptimal for large numbers of indices
            // we should do fieldName formatting here... should we also call ds.addColumnFromField?
            // or do we want to limit that?
            // slightly depends on our ability to display appropriately grouped options in the dropdown...
            const vals = Array.from(filteredIndices).map(getField);
            runInAction(() => {
                link.observableFields = vals;
            });
        }
    });
}

export function getRowsAsColumnsLinks(dataStore: DataStore) {
    const dataSources = window.mdv.chartManager.dataSources;
    if (dataStore.links) {
        return Object.keys(dataStore.links).map((linkedDsName) => {
            const links = dataStore.links[linkedDsName];
            if (links.rows_as_columns) {
                // first pass... there can be only one or zero.
                // how often will users actually want to link multiple dataSources in this way?
                // perhaps not often - but let's handle it so we don't have to change it later or have bugs.
                // UI should be simpler for the common case with a single linked dataSource.
                // Are there any crazy edge cases we should consider - like indirect links? links to self?
                // !! this should be right now, but we should test with multiple links.
                const linkedDs = dataSources.find(
                    (ds) => ds.name === linkedDsName,
                );
                if (!linkedDs) {
                    throw new Error();
                }
                // todo make sure the link is reasonably typed
                const link = links.rows_as_columns as RowsAsColslink;
                initRacListener(link, dataStore, linkedDs.dataStore);
                return { linkedDs, link };
            }
        });
    }
    return [];
}
