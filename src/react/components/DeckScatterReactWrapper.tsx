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
        config: ScatterPlotConfig & BaseConfig,
    ) {
        config.tooltip = config.tooltip || { show: false, column: config.param[0] }; //todo fix this
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
        const catCols = cols.filter((c) => c.datatype.match(/text/i));
        const settings = super.getSettings();

        let cats = this.dataStore.getColumnValues(c.param[2]) || [];
        cats = cats.map((x) => {
            return { t: x };
        });
        cats.push({ t: "None" });
        const catsValues = observable.array([cats, "t", "t"]);

        // What I would like is ability to
        // - change selected image at runtime.
        // - choose multiple categories on which to filter.
        const filters = c.category_filters.map((f) => {
            // what we really want is to have a type that is fairly specific to category filters...
            // able to handle an array of them with controls for adding/removing...
            const values = this.dataStore.columnIndex[f.column].values.slice();
            values.unshift("all");
            return {
                type: "multidropdown",
                label: `'${f.column}' filter`,
                current_value: f.category,
                values: [values],
                func: (v) => {
                    f.category = v;
                    c.category_filters = c.category_filters.slice();
                },
            };
        });
        //   ^^ kinda want a more react-y SettingsDialog for that...
        // todo this code should be shared with VivMDVReact
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
                type: "dropdown",
                label: "Shape",
                current_value: c.point_shape,
                values: [["circle", "square", "gaussian"]], //ugh
                func: (x) => {
                    c.point_shape = x;
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
            {
                type: "folder",
                label: "Contour Settings",
                current_value: [
                    {
                        type: "folder",
                        label: "Category selection",
                        current_value: [
                            //maybe 2-spaces format is better...
                            {
                                type: "dropdown",
                                label: "Contour parameter",
                                // current_value: c.contourParameter || this.dataStore.getColumnName(c.param[2]),
                                current_value: c.contourParameter || c.param[2],
                                values: [catCols, "name", "field"],
                                func: (x) => {
                                    if (x === c.contourParameter) return;
                                    // could we change 'cats' and have the dropdowns update?
                                    // was thinking this might mean a more general refactoring of the settings...
                                    c.contourParameter = c.param[2] = x; //this isn't causing useParamColumns to update...
                                    // but maybe it's not necessary if 'cats' is observable... fiddly to get right...
                                    const newCats = (
                                        this.dataStore.getColumnValues(x) || []
                                    ).map((t) => ({ t }));
                                    newCats.push({ t: "None" });
                                    catsValues[0] = newCats;
                                    //ru-roh, we're not calling the 'func's... mostly we just care about reacting to the change...
                                    //but setting things on config doesn't work anyway, because the dialog is based on this settings object...
                                    // c.category1 = c.category2 = null;  //maybe we can allow state to be invalid?
                                    //the dropdowns can set values to null if they're invalid rather than throw error?
                                    //is that a good idea?
                                },
                            },
                            {
                                type: "multidropdown",
                                label: "Contour Category 1",
                                current_value: c.category1 || "None",
                                // values: [cats, "t", "t"],
                                values: catsValues,
                                func: (x) => {
                                    if (x === "None") x = null;
                                    c.category1 = x;
                                },
                            },
                            {
                                type: "multidropdown",
                                label: "Contour Category 2",
                                current_value: c.category2 || "None",
                                // values: [cats, "t", "t"],
                                values: catsValues,
                                func: (x) => {
                                    if (x === "None") x = null;
                                    c.category2 = x;
                                },
                            },
                        ],
                        func: (x) => { },
                    },
                    {
                        type: "slider",
                        max: 25,
                        min: 1,

                        // doc: this.__doc__, //why?
                        current_value: c.contour_bandwidth,
                        label: "KDE Bandwidth",
                        continuous: true,
                        func: (x) => {
                            c.contour_bandwidth = x;
                        },
                    },
                    {
                        label: "Fill Contours",
                        type: "check",
                        current_value: c.contour_fill,
                        func: (x) => {
                            c.contour_fill = x;
                        },
                    },
                    {
                        type: "slider",
                        max: 1,
                        min: 0,
                        current_value: c.contour_intensity,
                        continuous: true,
                        label: "Fill Intensity",
                        func: (x) => {
                            c.contour_intensity = x;
                        },
                    },
                    {
                        type: "slider",
                        max: 1,
                        min: 0,

                        doc: this.__doc__,
                        current_value: c.contour_opacity,
                        continuous: false, //why so slow?
                        label: "Contour opacity",
                        func: (x) => {
                            c.contour_opacity = x ** 3;
                        },
                    },
                ],
                func: (x) => { },
            },
            {
                type: "folder",
                label: "Category Filters",
                current_value: filters,
                func: (x) => { },
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

export type VivMdvReactType = typeof DeckScatterReact;
// we rely on the side-effect of this import to register the chart type
///- will register otherwise, but this tricks HMR into working...
export default 42;
