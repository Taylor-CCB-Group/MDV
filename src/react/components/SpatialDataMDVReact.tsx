import { findRenderStackEntry, type RenderStack } from "@spatialdata/layers";
import type { EntryStatus } from "@spatialdata/vis";
import { action, comparer, computed, makeObservable, observable, reaction } from "mobx";
import type { SpatialDataGridConfig } from "@/charts/schemas/ChartConfigSchema";
import type DataStore from "@/datastore/DataStore";
import { allNumeric } from "@/lib/columnTypeHelpers";
import { g, toArray } from "@/lib/utils";
import { splitStorageIssues } from "@/react/spatialdata/data_quality";
import type { DeckBackend } from "@/react/spatialdata/deck_backend";
import { entryStatusErrorKey } from "@/react/spatialdata/entry_status";
import {
    DATA_POINTS_GROUP_ID,
    DECK_OVERLAY_VISIBILITY_FOLLOWERS,
    deckHostLayerId,
    deckIdFromHostLayerId,
} from "@/react/spatialdata/host_overlay_ids";
import type { ImageLayerRegistry } from "@/react/spatialdata/image_layer_registry";
import type { CanvasError, LayerFailure } from "@/react/spatialdata/layer_failures";
import type { PointsLayerRegistry } from "@/react/spatialdata/points_layer_registry";
import {
    createHostOnlyRenderStack,
    type DefaultImageLayerProps,
    migrateSavedRenderStack,
    syncOverlayFollowerVisibility,
} from "@/react/spatialdata/render_stack_defaults";
import { renderStackEntryDisplayName } from "@/react/spatialdata/render_stack_display";
import { sharedSpatialDataReader } from "@/react/spatialdata/spatial_data_reader";
import {
    type SpatialDataSerializableViewState,
    toSerializableSpatialDataViewState,
} from "@/react/spatialdata/spatialdata_config";
import BaseChart, { type MenuPip } from "../../charts/BaseChart";
import { getDensitySettings } from "../contour_state";
import { type CategoryFilter, type ScatterPlotConfig, scatterDefaults } from "../scatter_state";
import { createVivStores, type VivContextType } from "./avivatorish/state";
import { BaseReactChart } from "./BaseReactChart";
import { getCrosshairSetting } from "./pointerCrosshair";
import { getSpatialDataGridSettingsFolder } from "./regionGridUtils";
import SpatialDataChartRoot from "./SpatialDataMDVReactComponent";
import SpatialLayerDialogReactWrapper from "./SpatialLayerDialogReactWrapper";
import SpatialStatusPanel from "./SpatialStatusPanel";
import { getRenderBackendSettings, getSharedScatterSettings } from "./sharedScatterSettings";

export type SpatialDataMdvReactConfig = ScatterPlotConfig & {
    type: "SpatialDataMdvRegionReact";
    /** A cell per density field, or a cell per region (`SpatialRegionGrid`). */
    grid?: SpatialDataGridConfig | null;
    region: string;
    background_filter: CategoryFilter;
    renderStack?: RenderStack;
    viewState?: SpatialDataSerializableViewState | null;
    /** Experimental; null (or absent) follows the browser's default (see `deck_backend.ts`). */
    render_backend?: DeckBackend | null;
    /**
     * Props for the image layer added when the chart has no render stack; read when the chart is
     * created and not kept in its config. Set when a Viv chart is converted, for its channel settings.
     */
    defaultImageLayer?: DefaultImageLayerProps;
};

function adaptSpatialDataConfig(originalConfig: SpatialDataMdvReactConfig, dataStore: DataStore) {
    const config = { ...scatterDefaults, ...originalConfig };
    if (!dataStore.regions) {
        throw new Error("unexpected attempt to load spatial chart with no regions in datasource");
    }
    // Host-only stack must exist before makeAutoObservable(config) so in-place
    // entry edits stay reactive. Only synthesise a stack when the config never had one.
    if (config.renderStack === undefined) {
        config.renderStack = createHostOnlyRenderStack();
    } else {
        config.renderStack = migrateSavedRenderStack(config.renderStack);
        syncOverlayFollowerVisibility(config.renderStack);
    }
    // Must exist before makeAutoObservable(config) for later changes to be observed, and
    // survive initialiseChartConfig's JSON copy, so null rather than undefined.
    config.render_backend = config.render_backend ?? null;
    config.grid = config.grid ?? null;
    config.crosshair = config.crosshair ?? null;
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
     * so structural edits can refresh viewer inputs without replacing stack objects.
     */
    renderStackGeneration = 0;
    /**
     * Lightweight version token for in-place renderStack prop edits. This deliberately
     * does not recreate the spatial canvas passthrough used by image layers.
     */
    renderStackPropsGeneration = 0;
    /** True only until the first default image layer seed runs for a brand-new chart. */
    seedDefaultSpatialLayers: boolean;
    /** Props for that default image layer, from the config's `defaultImageLayer`. */
    readonly defaultImageLayerProps?: DefaultImageLayerProps;

    get viewerStore() {
        return this.vivStores?.viewerStore;
    }

    bumpRenderStackGeneration() {
        this.renderStackGeneration++;
    }

    bumpRenderStackPropsGeneration() {
        this.renderStackPropsGeneration++;
    }

    finishDefaultSpatialLayerSeed() {
        this.seedDefaultSpatialLayers = false;
    }

    imageLayerRegistry?: ImageLayerRegistry;
    setImageLayerRegistry(registry: ImageLayerRegistry | undefined) {
        this.imageLayerRegistry = registry;
    }

    pointsLayerRegistry?: PointsLayerRegistry;
    setPointsLayerRegistry(registry: PointsLayerRegistry | undefined) {
        this.pointsLayerRegistry = registry;
    }

    /**
     * The backend the canvas is drawing with, after any fallback to WebGL; null until the
     * device is chosen. Runtime only. For UI outside the canvas tree, such as the layer
     * dialog's image tone controls, which have no effect on WebGPU.
     */
    activeBackend: DeckBackend | null = null;
    /** Why WebGPU was requested but not used, when it wasn't. */
    backendFallbackReason: string | null = null;
    setActiveBackend(backend: DeckBackend | null, fallbackReason?: string) {
        this.activeBackend = backend;
        this.backendFallbackReason = fallbackReason ?? null;
    }

    get rendererLabel() {
        if (this.activeBackend === "webgpu") return "WebGPU";
        if (this.activeBackend === null) return "not chosen yet";
        if (this.backendFallbackReason) return `WebGL (WebGPU unavailable: ${this.backendFallbackReason})`;
        return "WebGL";
    }

    /** The GPU deck's device runs on, as luma reports it; for error reports. Runtime only. */
    gpuDescription: string | null = null;
    setGpuDescription(description: string | null) {
        this.gpuDescription = description;
    }

    /**
     * Render-stack entries that failed to draw and were turned off, by entry id. Runtime
     * only, never saved: turning the layer back on clears it, so a retry can fail and be
     * reported again. See `layer_failures.ts`.
     */
    layerFailures = new Map<string, LayerFailure>();
    /** Render errors that could not be tied to a layer. */
    canvasErrors: CanvasError[] = [];

    /** Record a failure and turn the entry off, so the rest of the canvas can draw again. */
    recordLayerFailure(failedId: string, failure: Omit<LayerFailure, "count" | "dismissed">) {
        // The data points' group has no row of its own in the layer dialog: its failure is
        // the data points', so turning them back on there retries it. It also outlives the
        // group, which comes and goes with the compositing settings.
        const entryId = failedId === DATA_POINTS_GROUP_ID ? deckHostLayerId("scatter") : failedId;
        const existing = this.layerFailures.get(entryId);
        if (existing) {
            existing.count++;
            return;
        }
        this.layerFailures.set(entryId, { ...failure, count: 1, dismissed: false });
        const stack = this.config.renderStack;
        if (!stack) return;
        let hidden = false;
        // An overlay that follows this one (the grey scatter) goes off with it, as its checkbox does.
        const deckId = deckIdFromHostLayerId(entryId);
        const followers = deckId ? (DECK_OVERLAY_VISIBILITY_FOLLOWERS[deckId] ?? []).map(deckHostLayerId) : [];
        for (const id of [entryId, ...followers]) {
            const entry = findRenderStackEntry(stack, id);
            if (entry?.visible) {
                entry.visible = false;
                hidden = true;
            }
        }
        if (hidden) this.bumpRenderStackGeneration();
    }

    clearLayerFailure(entryId: string) {
        this.layerFailures.delete(entryId);
    }

    recordCanvasError(message: string, stack?: string) {
        const existing = this.canvasErrors.find((item) => item.message === message);
        if (existing) existing.count++;
        else this.canvasErrors.push({ message, ...(stack ? { stack } : {}), count: 1, dismissed: false });
    }

    /**
     * What sd.js reports about each SpatialData layer, by entry id: why it could not load,
     * and what is degraded but still drawing (see `entry_status.ts`). Runtime only.
     */
    entryStatuses = new Map<string, EntryStatus>();
    setEntryStatus(entryId: string, status: EntryStatus | undefined) {
        if (status) this.entryStatuses.set(entryId, status);
        else this.entryStatuses.delete(entryId);
    }
    clearEntryStatuses() {
        this.entryStatuses.clear();
    }
    /** Load errors taken off the menu pip (`entryStatusErrorKey`); they stay listed. */
    dismissedLoadErrors = new Set<string>();

    /** Load a layer's element again (the renderer's `retryEntry`), while the canvas is up. */
    retryEntry?: (entryId: string) => Promise<void>;
    setRetryEntry(retry: ((entryId: string) => Promise<void>) | undefined) {
        this.retryEntry = retry;
    }

    /**
     * The chart's read of its store, for the canvas and the layer dialog (`useChartSpatialData`).
     * The page's reader, so charts on one store and Settings › Data quality share one read.
     * Not observable: the hook subscribes to the reader itself.
     */
    readonly spatialDataReader = sharedSpatialDataReader();

    /**
     * A layer's status without its element's storage issues, which stay quiet outside
     * Settings › Data quality (`data_quality.ts`), and those issues.
     */
    visibleEntryStatus(entryId: string) {
        return splitStorageIssues(this.entryStatuses.get(entryId));
    }

    /** The store the canvas is showing, for linking to it. Runtime only. */
    spatialDataUrl: string | null = null;
    setSpatialDataUrl(url: string | null) {
        this.spatialDataUrl = url;
    }

    /**
     * Failed layers and layer statuses by display name, and a one-line summary of the
     * chart's problems.
     */
    get renderStatus() {
        const stack = this.config.renderStack;
        const nameOf = (entryId: string) => {
            const entry = stack && findRenderStackEntry(stack, entryId);
            return entry ? renderStackEntryDisplayName(entry, stack) : entryId;
        };
        const failures = [...this.layerFailures].map(([entryId, failure]) => ({
            entryId,
            failure,
            name: nameOf(entryId),
        }));
        // storage issues are for Settings › Data quality: they neither count nor show here
        const statuses = [...this.entryStatuses.keys()].flatMap((entryId) => {
            const { status } = this.visibleEntryStatus(entryId);
            return status ? [{ entryId, status, name: nameOf(entryId) }] : [];
        });
        const renderErrorCount = failures.filter((item) => !item.failure.dismissed).length + this.canvasErrors.length;
        const failedToLoad = statuses.filter(({ status }) => status.errors.length > 0);
        const loadErrorCount = failedToLoad.filter(({ entryId, status }) =>
            status.errors.some((error) => !this.dismissedLoadErrors.has(entryStatusErrorKey(entryId, error))),
        ).length;
        const noticeCount = statuses.filter(({ status }) => status.notices.length > 0).length;
        const summary = [
            renderErrorCount > 0 && `${renderErrorCount} rendering ${renderErrorCount === 1 ? "problem" : "problems"}`,
            renderErrorCount === 0 &&
                failures.length > 0 &&
                `${failures.length} ${failures.length === 1 ? "layer" : "layers"} turned off`,
            failedToLoad.length > 0 &&
                `${failedToLoad.length} ${failedToLoad.length === 1 ? "layer" : "layers"} could not load`,
            noticeCount > 0 && `${noticeCount} layer ${noticeCount === 1 ? "notice" : "notices"}`,
        ]
            .filter(Boolean)
            .join(" · ");
        return { failures, statuses, errorCount: renderErrorCount + loadErrorCount, summary };
    }

    get statusPip(): MenuPip | null {
        const { failures, statuses, errorCount, summary } = this.renderStatus;
        if (errorCount > 0) return { level: "error", tooltip: summary };
        if (failures.length > 0 || statuses.length > 0) return { level: "notice", tooltip: summary };
        return null;
    }

    getDebugPanel() {
        return <SpatialStatusPanel chart={this} />;
    }

    getDebugSettings() {
        return getRenderBackendSettings(this.config);
    }

    /** Clear the error count on the menu pip; failing layers stay marked in the layer dialog. */
    dismissRenderErrors() {
        for (const failure of this.layerFailures.values()) failure.dismissed = true;
        this.canvasErrors = [];
        for (const [entryId, status] of this.entryStatuses) {
            for (const error of status.errors) this.dismissedLoadErrors.add(entryStatusErrorKey(entryId, error));
        }
    }

    constructor(dataStore: DataStore, div: HTMLDivElement, originalConfig: SpatialDataMdvReactConfig) {
        const seedDefaultSpatialLayers = originalConfig.renderStack === undefined;
        const config = adaptSpatialDataConfig(originalConfig, dataStore);
        const { defaultImageLayer } = config;
        config.defaultImageLayer = undefined;
        super(dataStore, div, config, SpatialDataChartRoot);
        this.seedDefaultSpatialLayers = seedDefaultSpatialLayers;
        if (seedDefaultSpatialLayers) this.defaultImageLayerProps = defaultImageLayer;
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
            renderStackGeneration: observable,
            bumpRenderStackGeneration: action,
            renderStackPropsGeneration: observable,
            bumpRenderStackPropsGeneration: action,
            seedDefaultSpatialLayers: observable,
            finishDefaultSpatialLayerSeed: action,
            imageLayerRegistry: observable.ref,
            setImageLayerRegistry: action,
            // `.ref` only: the engine notifies its own subscribers, and making it
            // deeply observable would have MobX walk a live render-path object.
            pointsLayerRegistry: observable.ref,
            setPointsLayerRegistry: action,
            activeBackend: observable,
            backendFallbackReason: observable,
            setActiveBackend: action,
            rendererLabel: computed,
            gpuDescription: observable,
            setGpuDescription: action,
            layerFailures: observable,
            canvasErrors: observable,
            recordLayerFailure: action,
            clearLayerFailure: action,
            recordCanvasError: action,
            dismissRenderErrors: action,
            // `.shallow`: a status is replaced whole when sd.js reports a change
            entryStatuses: observable.shallow,
            setEntryStatus: action,
            clearEntryStatuses: action,
            dismissedLoadErrors: observable,
            retryEntry: observable.ref,
            setRetryEntry: action,
            spatialDataUrl: observable,
            setSpatialDataUrl: action,
            renderStatus: computed,
            statusPip: computed,
        });
        this.reactionDisposers.push(
            reaction(
                () => this.statusPip,
                (pip) => this.setMenuPip(pip),
                { fireImmediately: true, equals: comparer.structural },
            ),
        );
        this.vivStores = createVivStores();
        this.addMenuIcon("fas fa-layer-group", "Manage Layers").addEventListener("click", () => this.openLayerDialog());
    }

    openLayerDialog() {
        if (!this.layerDialog) {
            this.layerDialog = new SpatialLayerDialogReactWrapper(this);
            this.dialogs.push(this.layerDialog);
        }
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

    /**
     * Chart-wide settings only. What the MDV data points and their density look like belongs
     * to their layers, and is edited with every other layer in the layer dialog
     * (`getDataLayerSettings`, `getDensityLayerSettings`), so the two are not split across dialogs.
     */
    getSettings() {
        const config = this.config;
        const dataStore = this.dataStore;
        const images = spatialRegionOptions(dataStore);

        return this.getGeneralSettings().concat([
            g({
                type: "button",
                label: "Layers",
                current_value: "Edit layers and data display…",
                func: () => this.openLayerDialog(),
            }),
            getCrosshairSetting(config),
            getSpatialDataGridSettingsFolder(config),
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
        ]);
    }

    /** Settings for the density layer (the `deck:density` overlay), shown in the layer dialog. */
    getDensityLayerSettings() {
        return getDensitySettings(this.config).current_value;
    }

    /** Settings for the MDV data points layer (the `deck:scatter` overlay), shown in the layer dialog. */
    getDataLayerSettings() {
        const config = this.config;
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
        const colorSettings = this.getColorSettingsFolder();

        return [
            ...(colorSettings ? [colorSettings] : []),
            ...getSharedScatterSettings(config, {
                chart: this,
                includePointShape: true,
            }),
            g({
                type: "folder",
                label: "Category Filters",
                current_value: filters,
            }),
        ];
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
        // null only exists to keep the key observable; leave it out of saved views
        if (config.render_backend === null) config.render_backend = undefined;
        if (config.grid === null) config.grid = undefined;
        if (config.crosshair === null) config.crosshair = undefined;
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
        const regionKey = extraConfig.region ?? spatialRegionOptions(dataStore)[0]?.value;
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
    name: "SpatialData.js viewer",
};

export { SpatialDataMdvReact };
export type SpatialDataMdvReactType = typeof SpatialDataMdvReact;
export default SpatialDataMdvReact;
