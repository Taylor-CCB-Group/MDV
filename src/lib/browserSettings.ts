import { action, observable } from "mobx";

/**
 * Settings for this browser rather than for a project: stored in localStorage, so they
 * apply to every project opened here and are never saved with a project or view. They are
 * observable, so anything reading them in an `observer` updates when one changes. The Debug
 * dialog lists every setting in `browserSettings`.
 *
 * Only choice settings exist so far; add other kinds alongside `defineChoiceSetting` when a
 * setting needs one.
 */

export type SettingChoice<T extends string> = { value: T; label: string };

export type ChoiceSettingDefinition<T extends string> = {
    /** localStorage key. */
    key: string;
    label: string;
    description?: string;
    choices: readonly SettingChoice<T>[];
    defaultValue: T;
    /** A value stored by an older form of this setting, read when there is no current value. */
    legacy?: () => T | undefined;
};

export type BrowserSetting<T extends string = string> = ChoiceSettingDefinition<T> & {
    get(): T;
    set(value: T): void;
};

function readStorage(key: string): string | null {
    try {
        return window.localStorage?.getItem(key) ?? null;
    } catch {
        return null;
    }
}

function writeStorage(key: string, value: string) {
    try {
        window.localStorage?.setItem(key, value);
    } catch (error) {
        console.warn(`Could not save browser setting '${key}'`, error);
    }
}

export function defineChoiceSetting<T extends string>(definition: ChoiceSettingDefinition<T>): BrowserSetting<T> {
    const isChoice = (value: string | null): value is T => definition.choices.some((choice) => choice.value === value);
    const stored = readStorage(definition.key);
    const box = observable.box<T>(isChoice(stored) ? stored : (definition.legacy?.() ?? definition.defaultValue));
    return {
        ...definition,
        get: () => box.get(),
        set: action((value: T) => {
            if (!isChoice(value)) return;
            box.set(value);
            writeStorage(definition.key, value);
        }),
    };
}

export type RenderBackend = "webgl" | "webgpu";

/** Which deck.gl backend SpatialData charts and scatter plots use when a chart does not choose one itself. */
export const renderBackendSetting = defineChoiceSetting<RenderBackend>({
    key: "mdv.settings.renderBackend",
    label: "Default rendering backend",
    description:
        "Used by SpatialData charts and scatter plots whose own Rendering backend setting is left " +
        "at the default. Where WebGPU is unavailable or fails, those charts use WebGL. On " +
        "WebGPU, image contrast and brightness are not yet supported. ?webgpu or ?webgpu=0 " +
        "in the URL overrides this for one page load.",
    choices: [
        { value: "webgpu", label: "WebGPU" },
        { value: "webgl", label: "WebGL" },
    ],
    defaultValue: "webgpu",
    // the dev flag this setting replaced
    legacy: () => (readStorage("MDV_WEBGPU") === "1" ? "webgpu" : undefined),
});

/** Diagnostics for WebGPU validation errors; see `src/react/spatialdata/webgpu_trace.ts`. */
export const webgpuTraceSetting = defineChoiceSetting<"off" | "on">({
    key: "mdv.settings.webgpuTrace",
    label: "WebGPU call tracing",
    description:
        "For diagnosing WebGPU errors: records recent command encoder calls, reports passes " +
        "left open when an encoder is finished, and adds both to a chart's Copy details. " +
        "Slows rendering a little. Takes effect when the page is reloaded.",
    choices: [
        { value: "off", label: "Off" },
        { value: "on", label: "On (after reload)" },
    ],
    defaultValue: "off",
});

/** Whether the Add Chart dialog offers chart types marked `legacy` in their `BaseChart.types` entry. */
export const legacyChartsSetting = defineChoiceSetting<"hide" | "show">({
    key: "mdv.settings.legacyCharts",
    label: "Legacy charts",
    description:
        "Show chart types that have been superseded (such as the classic scatter plots) in the " +
        "Add Chart dialog. Existing charts of those types are unaffected.",
    choices: [
        { value: "hide", label: "Hide" },
        { value: "show", label: "Show" },
    ],
    defaultValue: "hide",
});

/** Whether charts of a `legacy` type are loaded as their current equivalent; see `src/charts/legacyCharts.ts`. */
export const convertLegacyChartsSetting = defineChoiceSetting<"off" | "on">({
    key: "mdv.settings.convertLegacyCharts",
    label: "Convert legacy charts",
    description:
        "Load charts of a legacy type (such as the classic scatter plots) as their current " +
        "equivalent. The conversion is saved if the view is saved. Charts with no equivalent, " +
        "or that the equivalent could not show faithfully, are left as they are.",
    choices: [
        { value: "off", label: "Off" },
        { value: "on", label: "On" },
    ],
    defaultValue: "off",
});

/** Every browser setting, in the order the Debug dialog lists them. */
export const browserSettings: readonly BrowserSetting[] = [
    renderBackendSetting,
    webgpuTraceSetting,
    legacyChartsSetting,
    convertLegacyChartsSetting,
];
