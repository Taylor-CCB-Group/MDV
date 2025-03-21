import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import { action, makeObservable, observable } from "mobx";
import type { ColumnName, DataColumn, DataType } from "../../charts/charts";
import { loadColumn } from "@/dataloaders/DataLoaderUtil";
import { observer } from "mobx-react-lite";
import { scatterDefaults, type ScatterPlotConfig } from "../scatter_state";
import DeckScatterComponent from "./DeckScatterComponent";
import type { OrthographicViewState, OrbitViewState } from "deck.gl";
import { g } from "@/lib/utils";
import type DataStore from "@/datastore/DataStore";

const MainChart = observer(() => {
    return <DeckScatterComponent />;
});


export type DeckScatterConfig = ScatterPlotConfig & {
    viewState: OrthographicViewState | OrbitViewState;
} & BaseConfig;
const defaultViewState = {
        viewState: {
        target: [0, 0, 0],
        zoom: 0,
        minZoom: -50,
    }
}

class DeckScatterReact extends BaseReactChart<DeckScatterConfig> {
    /** set to true when this is the source of a viewState change etc to prevent circular update */
    ignoreStateUpdate = false;
    pendingRecenter = false;
    constructor(
        dataStore: DataStore,
        div: HTMLDivElement,
        originalConfig: DeckScatterConfig & BaseConfig,
    ) {
        // config.tooltip = config.tooltip || { show: false, column: config.param[0] }; //todo fix this
        const config = { ...scatterDefaults, ...defaultViewState, ...originalConfig };
        super(dataStore, div, config, MainChart);
        if (!originalConfig.viewState) this.pendingRecenter = true;
        //@ts-expect-error - pending colorBy type fix
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
            pendingRecenter: observable,
        });
    }
    colorBy?: (i: number) => [r: number, g: number, b: number];
    colorByColumn(col?: ColumnName) {
        if (!col) return this.colorByDefault();
        this.config.color_by = col;
        this.colorBy = this.getColorFunction(col, true);
    }
    colorByDefault() {
        //@ts-expect-error - pending colorBy type fix
        this.config.color_by = null;
        //@ts-expect-error - pending colorBy type fix
        this.colorBy = null;
    }
    getColorOptions() {
        return {
            colorby: "all",
        };
    }

    getConfig() {
        const c = super.getConfig() as DeckScatterConfig;
        // if we serialize the viewState, we get copies of things like `transition` properties
        // that don't de-serialize properly, so we explicity extract the properties we want
        const viewState = {
            target: c.viewState.target,
            zoom: c.viewState.zoom,
        }
        if (c.dimension === "3d") {
            const v = c.viewState as OrbitViewState;
            //@ts-expect-error - need to think about viewState types...
            viewState["rotationOrbit"] = v.rotationOrbit;
            //@ts-expect-error - need to think about viewState types...
            viewState["rotationX"] = v.rotationX;
        }
        return { ...c, viewState };
    }

    getSettings() {
        const c = this.config;
        const { tooltip } = c;
        const cols = this.dataStore.getColumnList() as DataColumn<DataType>[];
        //const catCols = cols.filter((c) => c.datatype.match(/text/i));
        const settings = super.getSettings();


        return settings.concat([
            //todo standard tooltip setting (with multi-choice)
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
                        await loadColumn(this.dataStore.name, cols[0].field);
                        tooltip.column = cols[0].field;
                    }
                },
            }),
            g({
                type: "dropdown",
                label: "Tooltip value",
                current_value: c.tooltip.column || cols[0].field,
                //@ts-expect-error - should be using "column" for the GiuSpecType rather than "dropdown"
                values: [cols, "name", "field"],
                func: async (c) => {
                    await loadColumn(this.dataStore.name, c);
                    tooltip.column = c;
                },
            }),
            g({
                type: "radiobuttons",
                label: "course radius",
                //@ts-expect-error - maybe we should allow numbers in radiobuttons? what actually happens internally?
                current_value: c.course_radius || 1,
                choices: [
                    ["0.1", "0.1"],
                    ["1", "1"],
                    ["10", "10"],
                    ["100", "100"],
                ],
                func: (x) => {
                    c.course_radius = Number.parseFloat(x);
                },
            }),
            g({
                type: "slider",
                label: "radius",
                current_value: c.radius || 5,
                min: 0,
                max: 20,
                continuous: true,
                func: (x) => {
                    c.radius = x;
                },
            }),
            g({
                type: "slider",
                label: "opacity",
                current_value: Math.sqrt(c.opacity || scatterDefaults.opacity),
                min: 0,
                max: 1,
                continuous: true,
                func: (x) => {
                    c.opacity = x * x;
                },
            }),
            g({
                type: "radiobuttons",
                label: "Action on Filter",
                choices: [
                    ["Hide Points", "hide"],
                    ["Gray Out Points", "grey"],
                ],
                current_value: c.on_filter,
                func: (x) => {
                    //@ts-ignore x is a string, but we have a narrow "hide" | "grey" type
                    c.on_filter = x;
                },
            }),
            g({
                type: "check",
                label: "zoom on filter",
                current_value: c.zoom_on_filter || false,
                func: (x) => {
                    c.zoom_on_filter = x;
                },
            }),
            g({
                type: "button",
                label: "Center Plot",
                //@ts-expect-error - no nay never no more
                current_value: null,
                func: () => {
                    // what should we do to trigger a re-center?
                    // we also want some logic so that this will happen when a new chart is loaded
                    // we don't really want it to be a serialized property... so we have observable state
                    // that isn't in the config, but is in the chart.
                    this.pendingRecenter = true;
                },
            })
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
