import { useMemo } from "react";
import { makeAutoObservable } from "mobx";
import type { AnyGuiSpec, GuiSpec } from "../../charts/charts";

type StableGuiSpec = AnyGuiSpec & { _stableId?: string };

function stampStableIds(specs: StableGuiSpec[], parentPath = "") {
    specs.forEach((spec, index) => {
        if (!spec._stableId) {
            spec._stableId = `${parentPath}${spec.type}:${spec.label}:${index}`;
        }
        if (spec.type === "folder" && Array.isArray(spec.current_value)) {
            stampStableIds(spec.current_value as StableGuiSpec[], `${spec._stableId}/`);
        }
    });
}

function filterSettingsByLabel(specs: StableGuiSpec[], searchTerm: string): StableGuiSpec[] {
    const query = searchTerm.trim().toLowerCase();
    if (!query) return specs;

    return specs.reduce<StableGuiSpec[]>((acc, spec) => {
        const labelMatches = spec.label.toLowerCase().includes(query);

        if (spec.type !== "folder") {
            // Include the setting only if the label matches
            if (labelMatches) acc.push(spec);
            return acc;
        }

        // Folder setting: recurse into children so nested matches are preserved.
        const filteredChildren = filterSettingsByLabel(spec.current_value as StableGuiSpec[], query);
        if (!labelMatches && filteredChildren.length === 0) {
            // Neither this folder nor any descendant matches.
            return acc;
        }

        const filteredFolder: GuiSpec<"folder"> & { _stableId?: string } = {
            ...spec,
            // If the folder label itself matches, keep all children so folder-name searches
            // show the full folder instead of an empty one.
            current_value: labelMatches ? spec.current_value : filteredChildren,
        };

        acc.push(filteredFolder);
        return acc;
    }, []);
}

/**
 * Returns render-ready top-level settings:
 * - filters by `label` (recursively for nested folders)
 * - wraps top-level non-folder elements into a single `General` folder
 */
export function useFilteredGroupedSettings(rawSettings: AnyGuiSpec[], searchTerm: string) {
    const topLevelSettings = useMemo(() => {
        // Split root settings into existing folders and top level settings.
        const groupedSettings = rawSettings.reduce<{
            folderSettings: AnyGuiSpec[];
            topLevelSettings: AnyGuiSpec[];
        }>(
            (acc, setting) => {
                if (setting.type === "folder") {
                    acc.folderSettings.push(setting);
                } else {
                    acc.topLevelSettings.push(setting);
                }
                return acc;
            },
            { folderSettings: [], topLevelSettings: [] },
        );

        // Keep UI tidy by collecting top level settings under a single "General" folder.
        const wrappedTopLevelSettings: StableGuiSpec[] =
            groupedSettings.topLevelSettings.length > 0
                ? [
                    {
                        type: "folder",
                        label: "General",
                        current_value: groupedSettings.topLevelSettings,
                    },
                ]
                : [];

        const settings = [...wrappedTopLevelSettings, ...groupedSettings.folderSettings] as StableGuiSpec[];
        stampStableIds(settings);
        return settings;
    }, [rawSettings]);

    const filteredSettings = useMemo(() => {
        // Filter after creating "General" folder, so folder-name searches includes it
        return filterSettingsByLabel(topLevelSettings, searchTerm);
    }, [topLevelSettings, searchTerm]);

    return useMemo(() => {
        const settings = filteredSettings.map(
            (setting) => ({
                setting,
                // Stable across filtering/search transitions.
                id: setting._stableId || `${setting.type}:${setting.label}`,
            }),
        );

        const wrap = { settings };
        makeAutoObservable(wrap);
        return wrap.settings;
    }, [filteredSettings]);
}

