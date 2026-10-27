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

/** Which deck.gl backend SpatialData charts use when a chart does not choose one itself. */
export const renderBackendSetting = defineChoiceSetting<RenderBackend>({
    key: "mdv.settings.renderBackend",
    label: "Default rendering backend",
    description:
        "Used by SpatialData charts whose own Rendering backend setting is left at the default. " +
        "Where WebGPU is unavailable or fails, those charts use WebGL. On WebGPU, image " +
        "contrast and brightness are not yet supported. ?webgpu or ?webgpu=0 in the URL " +
        "overrides this for one page load.",
    choices: [
        { value: "webgpu", label: "WebGPU" },
        { value: "webgl", label: "WebGL" },
    ],
    defaultValue: "webgpu",
    // the dev flag this setting replaced
    legacy: () => (readStorage("MDV_WEBGPU") === "1" ? "webgpu" : undefined),
});

/** Every browser setting, in the order the Debug dialog lists them. */
export const browserSettings: readonly BrowserSetting[] = [renderBackendSetting];
