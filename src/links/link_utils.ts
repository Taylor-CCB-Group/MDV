import { runInAction } from "mobx";
import ChartManager from "../charts/ChartManager";
import { ColumnName } from "../charts/charts";

type ChartLink = {
    type: "chart_columnval_link",
    id: string,
    source_chart: string,
    column: ColumnName,
    target_charts: string[],
    set_color: boolean,
    set_title: boolean,
    set_tooltip: boolean,
    set_legend: boolean,
    /** if true, the parm at given parm_index will be set...*/
    set_param: boolean,
    param_index?: number,
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
        console.log(`"chart_columnval_link" link column "${link.column}" not found in dataStore "${ds.name}"`);
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
                if (target.setToolTipColumn) target.setToolTipColumn(link.column);
                else if (target.config.tooltip) {
                    runInAction(() => target.config.tooltip.column = newValue);
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