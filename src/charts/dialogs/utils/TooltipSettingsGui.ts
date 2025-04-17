import BaseChart from "@/charts/BaseChart";
import { DataColumn, DataType } from "@/charts/charts";
import { loadColumn } from "@/dataloaders/DataLoaderUtil";
import { g } from "@/lib/utils";
import type { ScatterPlotConfig } from "@/react/scatter_state";

// this will be for more general tooltip settings rather than just scatterplot
export default function getTooltipSettings(config: ScatterPlotConfig, chart: BaseChart<any>) {
    // really don't want to be messing around with dataStore and such here...
    // first step is to separate out the tooltip settings from the chart config
    // and then we can sort out the implementation details
    const { tooltip } = config;
    const { dataStore } = chart;
    const cols = dataStore.getColumnList() as DataColumn<DataType>[];
    return g({
        type: "folder",
        label: "Tooltip",
        current_value: [
            g({
                type: "check",
                label: "Show Tooltip",
                current_value: tooltip.show,
                func: async (x: boolean) => {
                    tooltip.show = x;
                    if (!tooltip.column) {
                        const columnName = cols[0].field;
                        console.log(
                            "No tooltip column set, using first column:",
                            columnName,
                        );
                        await loadColumn(dataStore.name, cols[0].field);
                        tooltip.column = cols[0].field;
                    }
                },
            }),
            g({
                type: "dropdown",
                label: "Tooltip value",
                current_value: tooltip.column || cols[0].field,
                //@ts-expect-error - should be using "column" for the GiuSpecType rather than "dropdown"
                values: [cols, "name", "field"],
                func: async (c) => {
                    await loadColumn(dataStore.name, c);
                    tooltip.column = c;
                },
            }),
        ]
    });
}

function getTooltipSettings2(config: ScatterPlotConfig, callback?: () => void) {
    return g({
        type: "folder",
        label: "Tooltip",
        current_value: [
            g({
                type: "check",
                label: "Show tooltip",
                current_value: config.tooltip.show,
                func: (x) => {
                    config.tooltip.show = x;
                    callback?.();
                },
            }),
            g({
                type: "column",
                label: "Tooltip value",
                //@ts-expect-error - need a way of dealing with optional column...
                current_value: config.tooltip.column,
                func: async (x) => {
                    //@ts-expect-error pending tooltip column being FieldSpec(s)
                    config.tooltip.column = x;
                    callback?.();
                }
            }),
        ]
    });
}
