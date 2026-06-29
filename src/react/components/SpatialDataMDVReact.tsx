import type { RenderStack } from "@spatialdata/layers";
import { action, makeObservable, observable } from "mobx";

import BaseChart from "../../charts/BaseChart";
import type DataStore from "@/datastore/DataStore";
import { allNumeric } from "@/lib/columnTypeHelpers";
import { g, toArray } from "@/lib/utils";
import {
    type CategoryFilter,
    type ScatterPlotConfig,
    scatterDefaults,
} from "../scatter_state";
import { BaseReactChart } from "./BaseReactChart";
import {
    type VivContextType,
    createVivStores,
} from "./avivatorish/state";
import { getSharedScatterSettings } from "./sharedScatterSettings";
import type { ImageLayerRegistry } from "@/react/spatialdata/image_layer_registry";
import { createHostOnlyRenderStack } from "@/react/spatialdata/render_stack_defaults";
import SpatialLayerDialogReactWrapper from "./SpatialLayerDialogReactWrapper";
import SpatialDataChartRoot from "./SpatialDataMDVReactComponent";
import {
    type SpatialDataSerializableViewState,
    toSerializableSpatialDataViewState,
} from "@/react/spatialdata/spatialdata_config";

export type SpatialDataMdvReactConfig = ScatterPlotConfig & {
    type: "SpatialDataMdvRegionReact";
    region: string;
    background_filter: CategoryFilter;
    renderStack?: RenderStack;
    viewState?: SpatialDataSerializableViewState | null;
};

function adaptSpatialDataConfig(
    originalConfig: SpatialDataMdvReactConfig,
    dataStore: DataStore,
) {
    const config = { ...scatterDefaults, ...originalConfig };
    if (!dataStore.regions) {
        throw new Error("unexpected attempt to load spatial chart with no regions in datasource");
    }
    // Host-only stack must exist before makeAutoObservable(config) so in-place
    // entry edits stay reactive. Only synthesise a stack when the config never had one.
    if (config.renderStack === undefined) {
        config.renderStack = createHostOnlyRenderStack();
    }
    config.param = [...dataStore.regions.position_fields];
    if (typeof config.contourParameter === "string") {
        const column = dataStore.columnIndex[config.contourParameter];
        if (!column || allNumeric([column])) {
            config.contourParameter = dataStore.regions.region_field;
        }
    }
    return config;
}

class SpatialDataMdvReact extends BaseReactChart<SpatialDataMdvReactConfig> {
    declare dataStore: DataStore;
    vivStores: VivContextType;
    layerDialog?: SpatialLayerDialogReactWrapper;
    ignoreStateUpdate = false;
    /**
     * Intentional version token for in-place renderStack edits. The adapter reads this
     * so cosmetic edits can refresh viewer inputs without replacing stack objects and
     * causing layer/data churn.
     */
    renderStackGeneration = 0;
    /** True only until the first default image layer seed runs for a brand-new chart. */
    seedDefaultSpatialLayers: boolean;

    get viewerStore() {
        return this.vivStores?.viewerStore;
    }

    bumpRenderStackGeneration() {
        this.renderStackGeneration++;
    }

    finishDefaultSpatialLayerSeed() {
        this.seedDefaultSpatialLayers = false;
    }

    imageLayerRegistry?: ImageLayerRegistry;
    setImageLayerRegistry(registry: ImageLayerRegistry | undefined) {
        this.imageLayerRegistry = registry;
    }

    constructor(
        dataStore: DataStore,
        div: HTMLDivElement,
        originalConfig: SpatialDataMdvReactConfig,
    ) {
        const seedDefaultSpatialLayers = originalConfig.renderStack === undefined;
        const config = adaptSpatialDataConfig(originalConfig, dataStore);
        super(dataStore, div, config, SpatialDataChartRoot);
        this.seedDefaultSpatialLayers = seedDefaultSpatialLayers;
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
            renderStackGeneration: observable,
            bumpRenderStackGeneration: action,
            seedDefaultSpatialLayers: observable,
            finishDefaultSpatialLayerSeed: action,
            imageLayerRegistry: observable.ref,
            setImageLayerRegistry: action,
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

    colorByColumn(col?: SpatialDataMdvReactConfig["color_by"]) {
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
        const images = spatialRegionOptions(dataStore);

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
            config.viewState = toSerializableSpatialDataViewState(viewer.viewState);
        }
        if (this.config.renderStack) {
            config.renderStack = this.config.renderStack;
        }
        return config;
    }
}

function spatialRegionOptions(dataStore: DataStore) {
    return Object.keys(dataStore.regions?.all_regions ?? {})
        .filter((regionKey) => dataStore.regions?.all_regions[regionKey].spatial)
        .map((regionKey) => ({ name: regionKey, value: regionKey }));
}

BaseChart.types.SpatialDataMdvRegionReact = {
    required: (dataStore) => spatialRegionOptions(dataStore).length > 0,
    extra_controls: (dataStore) => {
        const regions = dataStore.regions;
        if (!regions) return [];
        return [
            {
                type: "dropdown",
                name: "region",
                label: dataStore.getColumnName(regions.region_field) ?? regions.region_field,
                values: spatialRegionOptions(dataStore),
            },
        ];
    },
    init: (config, dataStore, extraConfig) => {
        const regions = dataStore.regions;
        if (!regions) {
            throw new Error("unexpected attempt to initialise spatialdata chart with no regions in datasource");
        }
        const regionKey =
            extraConfig.region ?? spatialRegionOptions(dataStore)[0]?.value;
        if (!regionKey) {
            throw new Error("unexpected attempt to initialise spatialdata chart with no spatial regions in datasource");
        }
        config.color_by = regions.default_color;
        const colorColumn = dataStore.columnIndex[regions.default_color];
        if (!allNumeric([colorColumn])) {
            config.contourParameter = config.color_by;
        }
        config.param = [...regions.position_fields];
        config.background_filter = {
            column: regions.region_field,
            category: regionKey,
        };
        config.color_legend = { display: false };
        config.region = regionKey;
        config.title = regionKey;
        config.type = "SpatialDataMdvRegionReact";
        // nb, we're not initialising config.renderStack here, and it will need to be done during chart initialisation
        // but this means that if e.g. ChatMDV was generating a config object, it wouldn't necessarily need to fill in all those details
    },
    class: SpatialDataMdvReact,
    name: "SpatialData.js Image Viewer (experimental)",
};

export { SpatialDataMdvReact };
export type SpatialDataMdvReactType = typeof SpatialDataMdvReact;
export default SpatialDataMdvReact;
