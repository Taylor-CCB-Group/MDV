import { useMemo } from "react";
import { makeAutoObservable } from "mobx";
import { v4 as uuid } from "uuid";
import type { AnyGuiSpec, GuiSpec } from "../../charts/charts";

function filterSettingsByLabel(specs: AnyGuiSpec[], searchTerm: string): AnyGuiSpec[] {
    const query = searchTerm.trim().toLowerCase();
    if (!query) return specs;

    return specs.reduce<AnyGuiSpec[]>((acc, spec) => {
        const labelMatches = spec.label.toLowerCase().includes(query);

        if (spec.type !== "folder") {
            if (labelMatches) acc.push(spec);
            return acc;
        }

        const filteredChildren = filterSettingsByLabel(spec.current_value, query);
        if (!labelMatches && filteredChildren.length === 0) {
            return acc;
        }

        const filteredFolder: GuiSpec<"folder"> = {
            ...spec,
            current_value: filteredChildren,
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
    const filteredRawSettings = useMemo(() => {
        return filterSettingsByLabel(rawSettings, searchTerm);
    }, [rawSettings, searchTerm]);

    return useMemo(() => {
        const groupedSettings = filteredRawSettings.reduce<{
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

        const wrappedTopLevelSettings: AnyGuiSpec[] =
            groupedSettings.topLevelSettings.length > 0
                ? [
                    {
                        type: "folder",
                        label: "General",
                        current_value: groupedSettings.topLevelSettings,
                    },
                ]
                : [];

        const settings = [...wrappedTopLevelSettings, ...groupedSettings.folderSettings].map(
            (setting) => ({ setting, id: uuid() }),
        );

        const wrap = { settings };
        makeAutoObservable(wrap);
        return wrap.settings;
    }, [filteredRawSettings]);
}

