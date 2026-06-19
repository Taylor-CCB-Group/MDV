import type { RenderStack } from "@spatialdata/layers";
import { action, makeObservable, observable } from "mobx";

import BaseChart from "../../charts/BaseChart";
import type DataStore from "@/datastore/DataStore";
import { allNumeric } from "@/lib/columnTypeHelpers";
import { g, toArray } from "@/lib/utils";
import { scatterDefaults } from "../scatter_state";
import { BaseReactChart } from "./BaseReactChart";
import "../../charts/VivScatterPlot";
import {
    type VivContextType,
    applyDefaultChannelState,
    createVivStores,
} from "./avivatorish/state";
import type { VivMdvReactConfig } from "./VivMDVReact";
import { getSharedScatterSettings } from "./sharedScatterSettings";
import SpatialLayerDialogReactWrapper from "./SpatialLayerDialogReactWrapper";
import SpatialDataChartRoot from "./SpatialDataMDVReactComponent";

export type SpatialDataMdvReactConfig = VivMdvReactConfig & {
    renderStack?: RenderStack;
};

function adaptSpatialDataConfig(
    originalConfig: SpatialDataMdvReactConfig,
    dataStore: DataStore,
) {
    const config = { ...scatterDefaults, ...originalConfig };
    if (!dataStore.regions) {
        throw new Error("unexpected attempt to load spatial chart with no regions in datasource");
    }
    config.param = [...dataStore.regions.position_fields];
    if (typeof config.contourParameter === "string") {
        const column = dataStore.columnIndex[config.contourParameter];
        if (!column || allNumeric([column])) {
            config.contourParameter = dataStore.regions.region_field;
        }
    }
    config.viv = applyDefaultChannelState(config.viv);
    return config;
}

class SpatialDataMdvReact extends BaseReactChart<SpatialDataMdvReactConfig> {
    declare dataStore: DataStore;
    vivStores: VivContextType;
    layerDialog?: SpatialLayerDialogReactWrapper;
    ignoreStateUpdate = false;

    get viewerStore() {
        return this.vivStores?.viewerStore;
    }

    constructor(
        dataStore: DataStore,
        div: HTMLDivElement,
        originalConfig: SpatialDataMdvReactConfig,
    ) {
        const config = adaptSpatialDataConfig(originalConfig, dataStore);
        super(dataStore, div, config, SpatialDataChartRoot);
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
        });
        this.vivStores = createVivStores();
        this.addMenuIcon("fas fa-layer-group", "Manage Layers").addEventListener(
            "click",
            () => {
                if (!this.layerDialog) {
                    this.layerDialog = new SpatialLayerDialogReactWrapper(this);
                    this.dialogs.push(this.layerDialog);
                }
            },
        );
    }

    colorBy?: (i: number) => [r: number, g: number, b: number];

    colorByColumn(col?: VivMdvReactConfig["color_by"]) {
        if (!col) return this.colorByDefault();
        this.config.color_by = col;
        //@ts-expect-error legacy color_by options are normalised at runtime by BaseChart.
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
        const config = this.config;
        const settings = super.getSettings();
        const filters = config.category_filters.map((filter) => {
            const values = this.dataStore.columnIndex[filter.column]?.values?.slice();
            if (!values) throw `failed assertion that we should have a categorical '${filter.column}' here`;
            values.unshift("all");
            return g({
                type: "multidropdown",
                label: `'${filter.column}' filter`,
                current_value: toArray(filter.category),
                values: [values],
                func: (value) => {
                    filter.category = value;
                    config.category_filters = config.category_filters.slice();
                },
            });
        });
        const dataStore = this.dataStore;
        const imageRegionKeys = Object.keys(dataStore.regions?.all_regions ?? {}).filter(
            (regionKey) => dataStore.regions?.all_regions[regionKey].spatial,
        );
        const images = imageRegionKeys.map((regionKey) => ({
            name: regionKey,
            value: regionKey,
        }));

        return settings.concat([
            g({
                type: "dropdown",
                label: `SpatialData (${dataStore.getColumnName(dataStore.regions?.region_field)})`,
                current_value: config.region,
                values: [images, "name", "value"],
                func(value) {
                    if (config.title === config.region) {
                        config.title = value;
                    }
                    config.region = value;
                    config.background_filter.category = value;
                },
            }),
            ...getSharedScatterSettings(config, {
                chart: this,
                includeDensitySettings: true,
                includePointShape: true,
            }),
            g({
                type: "folder",
                label: "Category Filters",
                current_value: filters,
            }),
        ]);
    }

    getConfig() {
        const config = super.getConfig();
        if (this.vivStores) {
            const viewer = this.vivStores.viewerStore.getState();
            config.viv = {
                ...config.viv,
                viewerStore: {
                    viewState: viewer.viewState
                        ? {
                              target: viewer.viewState.target,
                              zoom: viewer.viewState.zoom,
                          }
                        : null,
                },
            };
        }
        if (this.config.renderStack) {
            config.renderStack = this.config.renderStack;
        }
        return config;
    }
}

BaseChart.types.SpatialDataMdvRegionReact = {
    ...BaseChart.types.VivMdvRegionReact,
    init: (config, dataStore, extraConfig) => {
        BaseChart.types.VivMdvRegionReact.init?.(config, dataStore, extraConfig);
        config.type = "SpatialDataMdvRegionReact";
    },
    class: SpatialDataMdvReact,
    name: "SpatialData.js Image Viewer (experimental)",
};

export { SpatialDataMdvReact };
export type SpatialDataMdvReactType = typeof SpatialDataMdvReact;
export default SpatialDataMdvReact;
