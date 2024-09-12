import BaseChart from "../../charts/BaseChart";
import { type BaseConfig, BaseReactChart } from "./BaseReactChart";
import { action, makeObservable, observable } from "mobx";
import type { ColumnName, DataColumn, DataType } from "../../charts/charts";
import { loadColumn } from "@/dataloaders/DataLoaderUtil";
import { observer } from "mobx-react-lite";
import { scatterDefaults, type ScatterPlotConfig } from "../scatter_state";
import DeckScatterComponent from "./DeckScatterComponent";

const MainChart = observer(() => {
    return <DeckScatterComponent />;
});


class DeckScatterReact extends BaseReactChart<ScatterPlotConfig> {
    /** set to true when this is the source of a viewState change etc to prevent circular update */
    ignoreStateUpdate = false;
    constructor(
        dataStore,
        div,
        originalConfig: ScatterPlotConfig & BaseConfig,
    ) {
        // config.tooltip = config.tooltip || { show: false, column: config.param[0] }; //todo fix this
        // config.opacity = config.opacity || scatterDefaults.opacity; // todo establish a better method for defaults
        const config = { ...scatterDefaults, ...originalConfig };
        super(dataStore, div, config, MainChart);
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
        });
    }
    colorBy?: (i: number) => [r: number, g: number, b: number];
    colorByColumn(col?: ColumnName) {
        if (!col) return this.colorByDefault();
        this.config.color_by = col;
        this.colorBy = this.getColorFunction(col, true);
    }
    colorByDefault() {
        this.config.color_by = null;
        this.colorBy = null;
    }
    getColorOptions() {
        return {
            colorby: "all",
        };
    }

    getSettings() {
        const c = this.config;
        const { tooltip } = c;
        const cols = this.dataStore.getColumnList() as DataColumn<DataType>[];
        //const catCols = cols.filter((c) => c.datatype.match(/text/i));
        const settings = super.getSettings();


        return settings.concat([
            {
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
                        await loadColumn(this.dataStore.name, cols[0].field);
                        tooltip.column = cols[0].field;
                    }
                },
            },
            {
                type: "dropdown",
                label: "Tooltip value",
                current_value: c.tooltip.column || cols[0].field,
                values: [cols, "name", "field"],
                func: async (c) => {
                    await loadColumn(this.dataStore.name, c);
                    tooltip.column = c;
                },
            },
            {
                type: "radiobuttons",
                label: "course radius",
                current_value: c.course_radius || 1,
                choices: [
                    [0.1, 0.1],
                    [1, 1],
                    [10, 10],
                    [100, 100],
                ],
                func: (x) => {
                    c.course_radius = x;
                },
            },
            {
                type: "slider",
                label: "radius",
                current_value: c.radius || 5,
                min: 0,
                max: 20,
                continuous: true,
                func: (x) => {
                    c.radius = x;
                },
            },
            {
                type: "slider",
                label: "opacity",
                current_value: Math.sqrt(c.opacity || scatterDefaults.opacity),
                min: 0,
                max: 1,
                continuous: true,
                func: (x) => {
                    c.opacity = x * x;
                },
            },
            {
                type: "check",
                label: "zoom on filter",
                current_value: c.zoom_on_filter || false,
                func: (x) => {
                    c.zoom_on_filter = x;
                },
            },
        ]);
    }
}

BaseChart.types["DeckScatter"] = {
    name: "2D Scatter Plot (new)",
    class: DeckScatterReact,
    params: [
        {
            type: "number",
            name: "x axis",
        },
        {
            type: "number",
            name: "y axis",
        },
    ]
};
BaseChart.types["DeckScatter3D"] = {
    name: "3D Scatter Plot (new)",
    class: DeckScatterReact,
    params: [
        {
            type: "number",
            name: "x axis",
        },
        {
            type: "number",
            name: "y axis",
        },
        {
            type: "number",
            name: "z axis",
        },
    ],
    init: (config: ScatterPlotConfig) => {
        config.dimension = "3d";
    }
};

export type VivMdvReactType = typeof DeckScatterReact;
// we rely on the side-effect of this import to register the chart type
///- will register otherwise, but this tricks HMR into working...
export default 42;
