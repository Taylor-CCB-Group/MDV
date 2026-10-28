import BaseChart from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import { action, makeObservable, observable } from "mobx";
import type { ColumnName } from "../../charts/charts";
import { observer } from "mobx-react-lite";
import {
    scatterAxisDefaults,
    scatterDefaults,
    type ScatterPlotConfig2D,
    type ScatterPlotConfig3D,
    type ScatterPlotConfig,
} from "../scatter_state";
import DeckScatterComponent from "./DeckScatterComponent";
import type { OrbitViewState } from "deck.gl";
import { g } from "@/lib/utils";
import type DataStore from "@/datastore/DataStore";
import getAxisGuiSpec from "@/charts/dialogs/utils/AxisSettingsGui";
import { getRenderBackendSettings, getSharedScatterSettings } from "./sharedScatterSettings";
import type { DeckBackend } from "@/react/spatialdata/deck_backend";

const MainChart = observer(() => {
    return <DeckScatterComponent />;
});

export type DeckScatterConfig = (ScatterPlotConfig2D | ScatterPlotConfig3D) & {
    /** null (or absent) follows the browser's default (see `deck_backend.ts`). */
    render_backend?: DeckBackend | null;
};
export type DeckScatterConfigWithRegion = DeckScatterConfig & { region?: string };
export const DECK_CONTOUR_SCATTER_TYPE = "DeckContourScatter";
const deckContourScatterTypes = new Set([DECK_CONTOUR_SCATTER_TYPE, "DeckDensity"]);

const defaultViewState = {
    viewState: {
        target: [0, 0, 0],
        zoom: 0,
        /** not expected to be user-tweakable as of this writing */
        minZoom: -50,
    },
};

function shouldMigrateLegacyContourParameter(originalConfig: DeckScatterConfig) {
    if (originalConfig.dimension === "3d") return false;
    if (!deckContourScatterTypes.has(originalConfig.type)) return false;
    return "category1" in originalConfig || "category2" in originalConfig;
}

function adaptConfig(originalConfig: DeckScatterConfig) {
    if (originalConfig.type === "wgl_3d_scatter_plot") {
        //! charts loaded from other configs may not have a dimension set
        originalConfig.dimension = "3d";
        if (originalConfig.course_radius === undefined) {
            originalConfig.course_radius = 10;
        }
    }

    // there is probably a less confusing way of writing this...
    //! originalConfig.dimension may be undefined, which lead to a bug with axis settings & broken charts
    // so if it is explicitly "3d", we have no axis settings, otherwise it will be "2d" | undefined
    const defaults = originalConfig.dimension !== "3d" ? { axis: scatterAxisDefaults } : {};
    const config = { ...scatterDefaults, ...defaults, ...defaultViewState, ...originalConfig };
    // Must exist before the config is made observable for later changes to be observed, and
    // survive initialiseChartConfig's JSON copy, so null rather than undefined.
    config.render_backend = config.render_backend ?? null;
    const legacyContourParameter = originalConfig.param?.[2];
    if (
        !config.contourParameter &&
        shouldMigrateLegacyContourParameter(originalConfig) &&
        legacyContourParameter &&
        !Array.isArray(legacyContourParameter)
    ) {
        config.contourParameter = legacyContourParameter;
    }
    return config;
}

const RendererStatus = observer(({ chart }: { chart: DeckScatterReact }) => {
    const { activeBackend, backendFallbackReason } = chart;
    const label =
        activeBackend === "webgpu"
            ? "WebGPU"
            : activeBackend === null
              ? "not chosen yet"
              : backendFallbackReason
                ? `WebGL (WebGPU unavailable: ${backendFallbackReason})`
                : "WebGL";
    return <p className="text-sm">Renderer: {label}</p>;
});

export class DeckScatterReact extends BaseReactChart<DeckScatterConfig> {
    /** set to true when this is the source of a viewState change etc to prevent circular update */
    ignoreStateUpdate = false;
    pendingRecenter = false;
    constructor(dataStore: DataStore, div: HTMLDivElement, originalConfig: DeckScatterConfig) {
        const config = adaptConfig(originalConfig);
        super(dataStore, div, config, MainChart);
        if (!originalConfig.viewState) this.pendingRecenter = true;
        //@ts-expect-error - pending colorBy type fix
        this.colorByColumn(config.color_by);
        makeObservable(this, {
            colorBy: observable,
            colorByColumn: action,
            colorByDefault: action,
            pendingRecenter: observable,
            activeBackend: observable,
            backendFallbackReason: observable,
            setActiveBackend: action,
        });
    }
    /** The backend the canvas is drawing with, after any fallback to WebGL; null until the device is chosen. */
    activeBackend: DeckBackend | null = null;
    /** Why WebGPU was requested but not used, when it wasn't. */
    backendFallbackReason: string | null = null;
    setActiveBackend(backend: DeckBackend | null, fallbackReason?: string) {
        this.activeBackend = backend;
        this.backendFallbackReason = fallbackReason ?? null;
    }
    getDebugPanel() {
        return <RendererStatus chart={this} />;
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
        };
        if (c.dimension === "3d") {
            const v = c.viewState;
            const orbitSerialized = viewState as OrbitViewState;
            orbitSerialized.rotationOrbit = v.rotationOrbit;
            orbitSerialized.rotationX = v.rotationX;
        }
        // null only exists to keep the key observable; leave it out of saved views
        const render_backend = c.render_backend ?? undefined;
        return { ...c, viewState, render_backend };
    }

    getSettings() {
        const c = this.config;
        const settings = super.getSettings();

        if (c.dimension === "2d") {
            const axisSettings = getAxisGuiSpec(c);
            settings.push(axisSettings);
        }
        return settings.concat([
            ...getSharedScatterSettings(c, {
                chart: this,
                includeDensitySettings: deckContourScatterTypes.has(c.type),
                includePointShape: true,
            }),
            ...getRenderBackendSettings(c),
            g({
                type: "button",
                label: "Center Plot",
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
        // {
        //     type: "_multi_column:number",
        //     name: "density fields"
        // }
    ],
};
BaseChart.types["DeckDensity"] = {
    name: "Scatter Plot + Density",
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
    ],
};
// The 2D scatter plot offered in the Add Chart dialog; density is an option that is off by default.
BaseChart.types[DECK_CONTOUR_SCATTER_TYPE] = {
    ...BaseChart.types["DeckDensity"],
    name: "2D Scatter Plot",
    allow_user_add: true,
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
        // todo 3d density with isosurface rendering...
    ],
    init: (config: ScatterPlotConfig) => {
        config.dimension = "3d";
    },
};
BaseChart.types["wgl_scatter_plot_dev"] = {
    name: "2D Scatter Plot (no density)",
    class: DeckScatterReact,
    legacy: true,
    params: [
        {
            type: "number",
            name: "x axis",
        },
        {
            type: "number",
            name: "y axis",
        },
    ],
};
BaseChart.types["wgl_3d_scatter_plot_dev"] = {
    // todo: remove (Deck) later, added to differentiate it with the existing 3d scatter plot
    name: "3D Scatter Plot (Deck)",
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
        //! charts loaded from other configs may not have a dimension set
        // and won't pass through the init function
        config.dimension = "3d";
    },
};

export type VivMdvReactType = typeof DeckScatterReact;
// we rely on the side-effect of this import to register the chart type
///- will register otherwise, but this tricks HMR into working...
export default 42;
