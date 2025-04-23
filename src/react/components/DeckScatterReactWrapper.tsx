import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import { action, makeObservable, observable } from "mobx";
import type { ColumnName, DataColumn, DataType } from "../../charts/charts";
import { loadColumn } from "@/dataloaders/DataLoaderUtil";
import { observer } from "mobx-react-lite";
import { scatterAxisDefaults, scatterDefaults, type ScatterPlotConfig2D, type ScatterPlotConfig3D, type ScatterPlotConfig } from "../scatter_state";
import DeckScatterComponent from "./DeckScatterComponent";
import type { OrthographicViewState, OrbitViewState } from "deck.gl";
import { g } from "@/lib/utils";
import type DataStore from "@/datastore/DataStore";
import getAxisGuiSpec from "@/charts/dialogs/utils/AxisSettingsGui";
import getTooltipSettings from "@/charts/dialogs/utils/TooltipSettingsGui";

const MainChart = observer(() => {
    return <DeckScatterComponent />;
});


export type DeckScatterConfig = ScatterPlotConfig2D | ScatterPlotConfig3D;
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
        originalConfig: DeckScatterConfig,
    ) {
        // config.tooltip = config.tooltip || { show: false, column: config.param[0] }; //todo fix this
        // there is probably a less confusing way of writing this...
        //! originalConfig.dimension may be undefined, which lead to a bug with axis settings & broken charts
        // so if it is explicitly "3d", we have no axis settings, otherwise it will be "2d" | undefined
        const defaults = originalConfig.dimension !== "3d" ? {axis: scatterAxisDefaults} : {};
        const config = { ...scatterDefaults, ...defaults, ...defaultViewState, ...originalConfig };
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
        // we should be able to make use of discriminated type here
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
        const settings = super.getSettings();

        if (c.dimension === "2d") {
            const axisSettings = getAxisGuiSpec(c);
            settings.push(axisSettings);
        }
        return settings.concat([
            getTooltipSettings(c),
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
            }),            
        ]);
    }
}

BaseChart.types["DeckScatter"] = {
    name: "2D Scatter Plot (new)",
    class: DeckScatterReact,
    allow_user_add: false,
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
    allow_user_add: false,
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
BaseChart.types["wgl_scatter_plot"] = {
    name: "2D Scatter Plot",
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
BaseChart.types["wgl_3d_scatter_plot"] = {
    name: "3D Scatter Plot",
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
