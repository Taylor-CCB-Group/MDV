import BaseChart from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import { action, makeObservable, observable } from "mobx";
import {
    type ROI,
    type VivConfig,
    type VivContextType,
    VivProvider,
    useViewerStore,
    useViewerStoreApi,
    applyDefaultChannelState,
    createVivStores,
} from "./avivatorish/state";
import "../../charts/VivScatterPlot"; //because we use the BaseChart.types object, make sure it's loaded.
import { useEffect } from "react";
import type { ColumnName } from "../../charts/charts";
import { useImage } from "./avivatorish/hooks";
import { VivScatter } from "./VivScatterComponent";
import { useImgUrl } from "../hooks";
import ColorChannelDialogReactWrapper from "./ColorChannelDialogReactWrapper";
import { observer } from "mobx-react-lite";
import { useChart } from "../context";
import type DataStore from "@/datastore/DataStore";
import { g, isArray, toArray } from "@/lib/utils";
import { scatterDefaults, type ScatterPlotConfig } from "../scatter_state";
import getTooltipSettings from "@/charts/dialogs/utils/TooltipSettingsGui";
import { allNumeric } from "@/lib/columnTypeHelpers";

function VivScatterChartRoot() {
    // to make this look like Avivator...
    // we use a VivProvider, with hooks that are mofified versions of Avivator's
    // VivProvider makes the vivStores available to the chart so that the chart can update & use the viewerStore et al.
    // in a way that is similar to Avivator - with the caveat that use of state in callbacks needs to be done
    // with saving a reference to e.g. `useViewerStoreApi()` and calling the `setState` method on that
    // rather than `useViewerStore.setState()`.
    const { vivStores } = useChart() as VivMdvReact;
    return (
        <VivProvider vivStores={vivStores}>
            <MainChart />
        </VivProvider>
    );
}

/** comparable to main `<Avivator />` component */
const MainChart = observer(() => {
    const imgUrl = useImgUrl();
    const isViewerLoading = useViewerStore((store) => store.isViewerLoading);
    const viewerStore = useViewerStoreApi();

    useEffect(() => {
        if (!imgUrl) throw "no image url";
        const source = { urlOrFile: imgUrl, description: "test" };
        viewerStore.setState({ source });
    }, [imgUrl, viewerStore.setState]);

    const source = useViewerStore((store) => store.source);
    // if (!source) throw "no image source"; //this is allowed to be undefined
    useImage(source);
    return !isViewerLoading && <VivScatter />;
});

export type TooltipConfig = {
    tooltip: {
        show: boolean;
        column?: ColumnName;
    };
};
export type CategoryFilter = {
    column: ColumnName;
    category: string | string[];
    // consider properties like 'invert' or 'exclude', or 'color'...
};
export type VivRoiConfig = {
    // making this 'type' very specific will let us infer the rest of the type, i.e.
    // `if (config.type === 'VivMdvRegionReact')` will narrow the type to VivScatterConfig
    // ... except that we also need to check the other condition, because 'string' could also be that.
    // so it's not completely ideal.
    type: "VivMdvRegionReact" | "viv_scatter_plot";
    region: string;
    background_filter: CategoryFilter;
    roi: ROI;
    viv: VivConfig;
    showJson: boolean;
    // json should come from associated config.region
    // json?: string, //for extra e.g. cell segmentation data - but we might want more than just a string...
    //image_properties: ChannelsState,
} & ScatterPlotConfig;

export type VivMdvReactConfig = ScatterPlotConfig &
    // 'simpler' version where you can just pass an image URL, don't need all the regions config...
    // not currently used.
    // { type: 'VivMdvReact', imageURL: string,
    //     overviewOn: boolean,
    //     image_properties: ChannelsState
    // }
    VivRoiConfig;
export type VivMDVReact = VivMdvReact;

function adaptConfig(originalConfig: VivMdvReactConfig, dataStore: DataStore) {
    const config = { ...scatterDefaults, ...originalConfig };
    // in future we might have something like an array of layers with potentially ways of describing parameters...
    // also need to address 3D case - intention to migrate to a different config schema before then.
    //@ts-expect-error contourParameter type
    if (!config.contourParameter) config.contourParameter = config.param[2];
    if (!config.contourParameter) throw "unexpected: no contourParameter";
    // **we need to mitigate the risk of this being non-categorical**
    const column = dataStore.columnIndex[config.contourParameter];
    if (!column) throw `unexpected: no column for contourParameter '${config.contourParameter}'`;
    // column type helpers could be more user friendly/comprehensive...
    if (allNumeric([column])) {
        console.warn(`contourParameter '${config.contourParameter}' is not categorical, using region field '${dataStore.getColumnName(dataStore.regions?.region_field)}' instead`);
        const defaultRegionField = dataStore.regions?.region_field;
        if (!defaultRegionField) throw "no default region field";
        config.contourParameter = defaultRegionField;
        config.param[2] = defaultRegionField;
    }
    // if (!config.)
    // === some dead code ===
    // if (config.type === 'VivMdvRegionReact') {
    //     // we don't use viv.image_properties, we use viv.channelsStore et al.
    //     // if (!config.viv.image_properties) config.viv.image_properties = DEFAUlT_CHANNEL_STATE;
    //     // else config.viv.image_properties = {...DEFAUlT_CHANNEL_STATE, ...config.viv.image_properties};
    //     // if (config.viv.image_properties) config.viv.image_properties = undefined;
    // } else if (config.type === 'VivMdvReact') {
    //     //unused
    //     if (config.overviewOn === undefined) config.overviewOn = false;
    //     // if (config.image_properties === undefined) config.image_properties = DEFAUlT_CHANNEL_STATE;
    //     // else config.viv.image_properties = {...DEFAUlT_CHANNEL_STATE, ...config.viv.image_properties};
    //     // if (config.viv.image_properties) config.viv.image_properties = undefined;
    // }

    // consider adapting `channels` from old format to new format...
    config.viv = applyDefaultChannelState(config.viv);
    return config;
}

class VivMdvReact extends BaseReactChart<VivMdvReactConfig> {
    colorDialog?: ColorChannelDialogReactWrapper;
    declare dataStore: DataStore;

    vivStores: VivContextType;
    get viewerStore() {
        return this.vivStores?.viewerStore;
    }

    /** set to true when this is the source of a viewState change etc to prevent circular update */
    ignoreStateUpdate = false;
    constructor(dataStore: DataStore, div: HTMLDivElement, originalConfig: VivMdvReactConfig) {
        // is this where I should be initialising vivStores? (can't refer to 'this' before super)
        // this.vivStores = createVivStores(this);
        const config = adaptConfig(originalConfig, dataStore);
        super(dataStore, div, config, VivScatterChartRoot);
        //@ts-expect-error color_by legacy options
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
        });
        this.addMenuIcon("fas fa-palette", "Alter Channels").addEventListener(
            "click",
            () => {
                if (!this.colorDialog) {
                    this.colorDialog = new ColorChannelDialogReactWrapper(this);
                    this.dialogs.push(this.colorDialog);
                }
            },
        );
        this.vivStores = createVivStores();
    }
    colorBy?: (i: number) => [r: number, g: number, b: number];
    colorByColumn(col?: ColumnName) {
        if (!col) return this.colorByDefault();
        this.config.color_by = col;
        this.colorBy = this.getColorFunction(col, true);
    }
    colorByDefault() {
        this.config.color_by = undefined;
        this.colorBy = undefined;
    }
    getColorOptions() {
        return {
            colorby: "all",
        };
    }

    getSettings() {
        const c = this.config;
        const { tooltip } = c;
        const cols = this.dataStore.getColumnList();// as DataColumn<DataType>[];
        const catCols = cols.filter((c) => c.datatype.match(/text/i));
        const settings = super.getSettings();

        //@ts-expect-error category column values...
        const ocats = this.dataStore.getColumnValues(c.param[2])?.slice() || [];
        const cats = ocats.map((x) => {
            return { t: x };
        });
        // could've sworn mobx observable had been working here at some point
        // (changing contourParameter should immediately update "Contour Category" dropdowns)... it isn't now.
        // and the type is dodgy - need to get on top of that with mobx in general.
        const catsValues = observable.array([cats, "t", "t"]) as unknown as [{t: string}[], "t", "t"];
        // const catsValues = [cats, "t", "t"];

        // What I would like is ability to
        // - change selected image at runtime.
        // - choose multiple categories on which to filter.
        const filters = c.category_filters.map((f) => {
            // what we really want is to have a type that is fairly specific to category filters...
            // able to handle an array of them with controls for adding/removing...
            const values = this.dataStore.columnIndex[f.column]?.values?.slice();
            if (!values) throw `failed assertion that we should have a categorical '${f.column}' here`;
            values.unshift("all");
            return g({
                type: "multidropdown",
                label: `'${f.column}' filter`,
                current_value: toArray(f.category),
                values: [values],
                func: (v) => {
                    f.category = v;
                    c.category_filters = c.category_filters.slice();
                },
            });
        });
        //   ^^ kinda want a more react-y SettingsDialog for that...
        // todo make sure associated json etc switches when region changes
        const ds = this.dataStore;
        const imageRegionKeys = Object.keys(ds.regions?.all_regions).filter(
            (r) => ds.regions?.all_regions[r].viv_image,
        );
        const images = imageRegionKeys.map((r) => ({ name: r, value: r }));


        return settings.concat([
            g({
                type: "dropdown",
                label: `Image (${ds.getColumnName(ds.regions?.region_field)})`,
                current_value: c.region,
                values: [images, "name", "value"],
                func(v) {
                    console.log("setting image region:", v);
                    //nb, 'this' is not the chart...
                    if (c.title === c.region) {
                        c.title = v;
                    }
                    // c.viv.url = v;
                    // ideally, c.region could be enough to also know what json etc is relevant...
                    c.region = v;
                    // background_filter.category should be inferred from the region...
                    // (this shouldn't be the responsibility of this function)
                    c.background_filter.category = v;
                },
            }),
            getTooltipSettings(c),
            g({
                type: "dropdown",
                label: "Shape",
                current_value: c.point_shape,
                values: [["circle", "square", "gaussian"]], //ugh
                func: (x) => {
                    //@ts-ignore we have a very restricted type for point_shape - could think about supporting that
                    c.point_shape = x;
                },
            }),
            g({
                type: "radiobuttons",
                label: "course radius",
                current_value: `${c.course_radius || 1}`,
                choices: [
                    [0.1, 0.1],
                    [1, 1],
                    [10, 10],
                    [100, 100],
                ].map(a => [`${a[0]}`, `${a[1]}]`]),
                func: (x) => {
                    //@ts-check !todo - maybe we should allow radiobuttons to use different datatypes
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
                type: "check",
                label: "zoom on filter",
                current_value: c.zoom_on_filter || false,
                func: (x) => {
                    c.zoom_on_filter = x;
                },
            }),
            g({
                type: "check",
                label: "show json layer",
                current_value: c.showJson || false,
                func: (x) => {
                    c.showJson = x;
                },
            }),
            g({
                type: "folder",
                label: "Density Visualisation",
                current_value: [
                    g({
                        type: "folder",
                        label: "Category selection",
                        current_value: [
                            //maybe 2-spaces format is better...
                            g({
                                type: "dropdown", //todo, make this "column" and fix odd behaviour with showing the value...
                                //todo: make the others be "category_selection" or something (which we don't have yet as a GuiSpec type)
                                label: "Contour parameter",
                                // current_value: c.contourParameter || this.dataStore.getColumnName(c.param[2]),
                                //@ts-expect-error contourParameter type
                                current_value: c.contourParameter || c.param[2],
                                values: [catCols, "name", "field"],
                                func: (x) => {
                                    if (x === c.contourParameter) return;
                                    // could we change 'cats' and have the dropdowns update?
                                    // was thinking this might mean a more general refactoring of the settings...
                                    if (!isArray(c.param)) throw "expected param array";
                                    c.contourParameter = c.param[2] = x; //this isn't causing useParamColumns to update...
                                    // but maybe it's not necessary if 'cats' is observable... fiddly to get right...
                                    const newCats = (
                                        this.dataStore.getColumnValues(x) || []
                                    ).map((t) => ({ t }));
                                    // newCats.push({ t: "None" });
                                    console.warn("changing contour parameter isn't properly updating dropdowns as of this writing...");
                                    catsValues[0] = newCats;
                                    //ru-roh, we're not calling the 'func's... mostly we just care about reacting to the change...
                                    //but setting things on config doesn't work anyway, because the dialog is based on this settings object...
                                    // c.category1 = c.category2 = null;  //maybe we can allow state to be invalid?
                                    //the dropdowns can set values to null if they're invalid rather than throw error?
                                    //is that a good idea?
                                },
                            }),
                            g({
                                type: "multidropdown",
                                label: "Contour Category 1",
                                current_value: toArray(c.category1 || "None"),
                                // values: [cats, "t", "t"],
                                values: catsValues,
                                func(x) {
                                    // if (x === "None") x = null;
                                    c.category1 = x;
                                },
                            }),
                            g({
                                type: "multidropdown",
                                label: "Contour Category 2",
                                current_value: toArray(c.category2 || "None"),
                                // values: [cats, "t", "t"],
                                values: catsValues,
                                func(x) {
                                    // if (x === "None") x = null;
                                    c.category2 = x;
                                },
                            }),
                        ],
                    }),
                    g({
                        type: "slider",
                        max: 25,
                        min: 1,

                        // doc: this.__doc__, //why?
                        current_value: c.contour_bandwidth,
                        label: "KDE Bandwidth",
                        continuous: true,
                        func(x) {
                            c.contour_bandwidth = x;
                        },
                    }),
                    g({
                        label: "Fill Contours",
                        type: "check",
                        current_value: c.contour_fill,
                        func(x) {
                            c.contour_fill = x;
                        },
                    }),
                    g({
                        type: "slider",
                        max: 1,
                        min: 0,
                        current_value: c.contour_intensity,
                        continuous: true,
                        label: "Fill Intensity",
                        func(x) {
                            c.contour_intensity = x;
                        },
                    }),
                    g({
                        type: "slider",
                        max: 1,
                        min: 0,
                        current_value: c.contour_opacity,
                        continuous: false, //why so slow?
                        label: "Contour opacity",
                        func(x) {
                            c.contour_opacity = x ** 3;
                        },
                    }),
                ],
            }),
            g({
                type: "folder",
                label: "Category Filters",
                current_value: filters,
            }),
            // ...filters,
            // no longer using PictureInPictureViewer - up for review as could be useful
            // {
            //     type: "check",
            //     label: "overview",
            //     current_value: c.overviewOn || false,
            //     func: x => {
            //         c.overviewOn = x;
            //     }
            // }
        ]);
    }
    getConfig() {
        const config = super.getConfig();
        // todo this is an experimental version of serialisation, WIP / subject to change.
        if (this.vivStores) {
            const { viewerStore, channelsStore, imageSettingsStore } =
                this.vivStores;
            const channels = channelsStore.getState();
            const viewer = viewerStore.getState();
            // const imageSettings = imageSettingsStore.getState();
            // omit some things like width, height...
            // for 3D we'd want rotation etc, but if we keep entire viewState it causes problems (stuck zoom)
            const viewState = {
                target: viewer.viewState.target,
                zoom: viewer.viewState.zoom,
            };
            const viv = {
                viewerStore: {
                    // could be that we want this to be where `'source'` comes from...
                    // we could be interested in use3d, useLens, useLinkedView... those would come from imageSettingsStore
                    // we might want an expanded version of 'pixelValues'...
                    viewState,
                },
                // we could call this `image_properties` and make it compatible with the 'legacy format' from VivScatterPlot
                // or we could parse it into a nicer viv.channels[] array-of-objects format, and parse back when we load.
                channelsStore: {
                    channelsVisible: channels.channelsVisible,
                    colors: channels.colors,
                    contrastLimits: channels.contrastLimits,
                    brightness: channels.brightness,
                    contrast: channels.contrast,
                    domains: channels.domains,
                    selections: channels.selections,
                },
                imageSettingsStore: {
                    // could be interested in lensEnabled, zoomLock, panLock, etc.
                },
            };
            // console.table(viv);
            config.viv = {
                ...config.viv,
                ...viv,
            };
        }
        return config;
    }
}

BaseChart.types["VivMdvRegionReact"] = {
    ...BaseChart.types["viv_scatter_plot"], //this is doing something that means my default radius isn't being used...
    init: (config, ds, ec) => {
        const base = BaseChart.types["viv_scatter_plot"];
        if (!base || !base.init) throw "no base viv_scatter_plot"; //may well want to change this behaviour soon
        base.init(config, ds, ec);
        config.radius = scatterDefaults.radius;
    },
    class: VivMdvReact,
    name: "Viv Scatter Plot (react)",
};

export type VivMdvReactType = typeof VivMdvReact;
// we rely on the side-effect of this import to register the chart type
///- will register otherwise, but this tricks HMR into working...
export default 42;
