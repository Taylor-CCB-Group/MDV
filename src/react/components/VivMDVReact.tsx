import BaseChart from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import { action, makeObservable, observable } from "mobx";
import { BaseDialog } from "../../utilities/Dialog";
import { DEFAUlT_CHANNEL_STATE, ROI, VivConfig, VivContextType, VivProvider, useViewerStore, useViewerStoreApi } from "./avivatorish/state";
import "../../charts/VivScatterPlot"; //because we use the BaseChart.types object, make sure it's loaded.
import { useEffect } from "react";
import type { ColumnName, DataColumn } from "../../charts/charts";
import { useImage } from "./avivatorish/hooks";
import { VivScatter } from "./VivScatterComponent";
import { useImgUrl } from "../hooks";
import ColorChannelDialogReactWrapper from "./ColorChannelDialogReactWrapper";

function ReactTest() {
    // to make this look more like Avivator...
    // we probably don't want OmeTiffProvider to be a thing...
    // we should use a VivProvider, with hooks that look more like Avivator's
    // so VivProvider should have whatever is necessary to adapt our config to that
    // and we'd useLoader() as opposed to useOmeTiff()
    // ... and hopefully our version of Avivator hooks will have better types ...
    return (
    <VivProvider>
        <MainChart />
    </VivProvider>
    )
}



/** comparable to main `<Avivator />` component */
const MainChart = () => {
    const imgUrl = useImgUrl();
    const isViewerLoading = useViewerStore(store => store.isViewerLoading);
    const viewerStore = useViewerStoreApi();

    useEffect(() => {
        if (!imgUrl) throw 'no image url';
        const source = { urlOrFile: imgUrl, description: 'test' };
        viewerStore.setState({ source });
    }, [imgUrl]);

    const source = useViewerStore(store => store.source);
    useImage(source);
    return (!isViewerLoading && <VivScatter />);
};



export type TooltipConfig = {
    tooltip: {
        show: boolean,
        column?: ColumnName,
    }
};
type CategoryFilter = {
    column: ColumnName,
    category: string | string[],
    // consider properties like 'invert' or 'exclude', or 'color'...
}
//viewState should be persisted... maybe a way of saving different snapshots?
//could we infer or something to avoid having to repeat this?
export type ScatterPlotConfig = {
    radius: number,
    opacity: number,
    color_by: ColumnName,
    color_legend: {
        display: boolean,
        // todo: add more options here...
    },
    category_filters: Array<CategoryFilter>,
    zoom_on_filter: boolean,
    point_shape: "circle" | "square" | "gaussian"
} & TooltipConfig;
const scatterDefaults: ScatterPlotConfig = {
    radius: 10,
    opacity: 1,
    color_by: null,
    color_legend: {
        display: false,
    },
    tooltip: {
        show: false,
    },
    category_filters: [],
    zoom_on_filter: false,
    point_shape: "circle",
};
export type VivRoiConfig = {
    // making this 'type' very specific will let us infer the rest of the type, i.e.
    // `if (config.type === 'VivMdvRegionReact')` will narrow the type to VivScatterConfig
    // ... except that we also need to check the other condition, because 'string' could also be that.
    // so it's not completely ideal.
    type: "VivMdvRegionReact" | "viv_scatter_plot",
    background_filter: CategoryFilter,
    roi: ROI,
    viv: VivConfig,
    //image_properties: ChannelsState,
} & ScatterPlotConfig;

export type VivMdvReactConfig = ScatterPlotConfig & (
    // 'simpler' version where you can just pass an image URL, don't need all the regions config...
    // not currently used.
    // { type: 'VivMdvReact', imageURL: string, 
    //     overviewOn: boolean, 
    //     image_properties: ChannelsState 
    // }
    | VivRoiConfig
) & { channel: number };
export type VivMDVReact = VivMdvReact;
class VivMdvReact extends BaseReactChart<VivMdvReactConfig> {
    colorDialog: any;

    vivStores?: VivContextType;
    get viewerStore() {
        return this.vivStores?.viewerStore;
    }

    /** set to true when this is the source of a viewState change etc to prevent circular update */
    ignoreStateUpdate = false;
    constructor(dataStore, div, config) {
        // todo better default config
        config = {...scatterDefaults, ...config};
        if (!config.channel) config.channel = 0;
        if (config.type === 'VivMdvRegionReact') {
            if (!config.viv.image_properties) config.viv.image_properties = DEFAUlT_CHANNEL_STATE;
        } else if (config.type === 'VivMdvReact') {
            if (config.overviewOn === undefined) config.overviewOn = false;
            if (config.image_properties === undefined) config.image_properties = DEFAUlT_CHANNEL_STATE;
        }
        // is this where I should be initialising vivStores? (can't refer to 'this' before super)
        // this.vivStores = createVivStores(this);
        super(dataStore, div, config, ReactTest);
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
        });
        this.addMenuIcon("fas fa-palette", "Alter Channels").addEventListener("click", (e) => {
            if (!this.colorDialog) {
                this.colorDialog = new ColorChannelDialogReactWrapper(this);
                this.dialogs.push(this.colorDialog);
            }
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
            colorby: "all"
        }
    }

    getSettings() {
        const c = this.config;
        const { tooltip } = c;
        const cols = this.dataStore.getColumnList() as DataColumn<any>[];
        const settings = super.getSettings();
        // What I would like is ability to
        // - change selected image at runtime.
        // - choose multiple categories on which to filter.
        const filters = c.category_filters.map(f => {
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
                }
            }
        });
        //   ^^ kinda want a more react-y SettingsDialog for that...
        // todo switch image, <- more coherent region logic...
        // const ds = this.dataStore;
        // const images = [];
        // for (const r in ds.regions.all_regions) {
        //     const viv = ds.regions.all_regions[r].viv_image;
        //     if (viv) {
        //         const x = viv.url || viv.file;
        //         images.push({name: x, value: x});
        //     }
        // }
        return settings.concat([
            // {
            //     type: "dropdown",
            //     label: `Image (${ds.getColumnName(ds.regions.region_field)})`,
            //     current_value: c.viv.url,
            //     values: [images, 'name', 'value'],
            //     func: (v) => {
            //         c.viv.url = v;
            //         c.background_filter.category = v;
            //     }
            // },
            {
                type: "check",
                label: "Show Tooltip",
                current_value: tooltip.show,
                func: (x: boolean) => {
                    tooltip.show = x;
                }
            },
            {
                type: "dropdown",
                label: "Tooltip value",
                current_value: c.tooltip.column || cols[0].field,
                values: [cols, "name", "field"],
                func: (c) => {
                    tooltip.column = c;
                }
            },
            {
                type: "dropdown",
                label: "Shape",
                current_value: c.point_shape,
                values: [["circle", "square", "gaussian"]], //ugh
                func: x => {
                    c.point_shape = x;
                }
            },
            {
                type: "slider",
                label: "radius",
                current_value: c.radius || 10,
                min: 0,
                max: 500,
                continuous: true,
                func: x => {
                    c.radius = x;
                }
            },
            {
                type: "slider",
                label: "opacity",
                current_value: Math.sqrt(c.opacity || scatterDefaults.opacity),
                min: 0,
                max: 1,
                continuous: true,
                func: x => {
                    c.opacity = x*x;
                }
            },
            {
                type: "check",
                label: "zoom on filter",
                current_value: c.zoom_on_filter || false,
                func: x => {
                    c.zoom_on_filter = x;
                }
            },
            {
                type: 'folder',
                label: 'Category Filters',
                current_value: filters,
                func: (x) => {},
            },
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
            const { viewerStore, channelsStore, imageSettingsStore } = this.vivStores;
            const channels = channelsStore.getState();
            const viewer = viewerStore.getState();
            const imageSettings = imageSettingsStore.getState();
            const viv = {
                viewerStore: {
                    // could be that we want this to be where `'source'` comes from...
                    // we could be interested in use3d, useLens, useLinkedView... those would come from imageSettingsStore
                    // we might want an expanded version of 'pixelValues'...
                    viewState: viewer.viewState,
                },
                // we could call this `image_properties` and make it compatible with the 'legacy format' from VivScatterPlot
                // or we could parse it into a nicer viv.channels[] array-of-objects format, and parse back when we load.
                channelsStore: {
                    channelsVisible: channels.channelsVisible,
                    colors: channels.colors,
                    contrastLimits: channels.contrastLimits,
                    domains: channels.domains,
                    selections: channels.selections,
                },
                imageSettingsStore: {
                    // could be interested in lensEnabled, zoomLock, panLock, etc.
                },
            }
            console.table(viv);
            config.viv = {
                ...config.viv,
                ...viv,
            }
        }
        return config;
    }
}

BaseChart.types["VivMdvRegionReact"] = {
    ...BaseChart.types["viv_scatter_plot"], //this is doing something that means my default radius isn't being used...
    init: (config, ds, ec) => {
        BaseChart.types["viv_scatter_plot"].init(config, ds, ec);
        config.radius = scatterDefaults.radius;
    },
    "class": VivMdvReact,
    name: "Viv Scatter Plot (react)",
}

export type VivMdvReactType = typeof VivMdvReact;
// we rely on the side-effect of this import to register the chart type
///- will register otherwise, but this tricks HMR into working...
export default 42;