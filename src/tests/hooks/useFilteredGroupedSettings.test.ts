import { describe, expect, test } from "vitest";
import { renderHook } from "@testing-library/react";
import { g } from "@/lib/utils";
import type { AnyGuiSpec, GuiSpec } from "@/charts/charts";
import { useFilteredGroupedSettings } from "../../react/hooks/useFilteredGroupedSettings";

function getFolder(
    settings: Array<{ setting: AnyGuiSpec; id: string }>,
    label: string,
): GuiSpec<"folder"> | undefined {
    const match = settings.find(
        (entry) => entry.setting.type === "folder" && entry.setting.label === label,
    )?.setting;
    if (!match || match.type !== "folder") return undefined;
    return match;
}

describe("useFilteredGroupedSettings", () => {
    test("wraps top-level non-folder settings into General folder", () => {
        const rawSettings: AnyGuiSpec[] = [
            g({ type: "slider", label: "Opacity", current_value: 0.6 }),
            g({
                type: "folder",
                label: "Axes",
                current_value: [g({ type: "check", label: "Show grid", current_value: true })],
            }),
            g({ type: "check", label: "Show legend", current_value: true }),
        ];

        const { result } = renderHook(() =>
            useFilteredGroupedSettings(rawSettings, ""),
        );

        const general = getFolder(result.current, "General");
        const axes = getFolder(result.current, "Axes");

        expect(general).toBeDefined();
        expect(general?.current_value.map((item) => item.label)).toEqual([
            "Opacity",
            "Show legend",
        ]);
        expect(axes).toBeDefined();
    });

    test("matches General folder by name when searching", () => {
        const rawSettings: AnyGuiSpec[] = [
            g({ type: "slider", label: "Opacity", current_value: 0.4 }),
            g({ type: "check", label: "Show legend", current_value: true }),
        ];

        const { result } = renderHook(() =>
            useFilteredGroupedSettings(rawSettings, "general"),
        );

        const general = getFolder(result.current, "General");
        expect(result.current).toHaveLength(1);
        expect(general).toBeDefined();
        expect(general?.current_value.map((item) => item.label)).toEqual([
            "Opacity",
            "Show legend",
        ]);
    });

    test("keeps full folder content when folder label matches search", () => {
        const rawSettings: AnyGuiSpec[] = [
            g({
                type: "folder",
                label: "Axes",
                current_value: [
                    g({ type: "slider", label: "X min", current_value: 0 }),
                    g({ type: "slider", label: "Y min", current_value: 0 }),
                ],
            }),
        ];

        const { result } = renderHook(() =>
            useFilteredGroupedSettings(rawSettings, "axes"),
        );

        const axes = getFolder(result.current, "Axes");
        expect(axes).toBeDefined();
        expect(axes?.current_value.map((item) => item.label)).toEqual([
            "X min",
            "Y min",
        ]);
    });

    test("keeps only matching descendants when folder does not match", () => {
        const rawSettings: AnyGuiSpec[] = [
            g({
                type: "folder",
                label: "Axes",
                current_value: [
                    g({ type: "slider", label: "X min", current_value: 0 }),
                    g({ type: "slider", label: "Y min", current_value: 0 }),
                ],
            }),
        ];

        const { result } = renderHook(() =>
            useFilteredGroupedSettings(rawSettings, "x min"),
        );

        const axes = getFolder(result.current, "Axes");
        expect(axes).toBeDefined();
        expect(axes?.current_value.map((item) => item.label)).toEqual(["X min"]);
    });

    test("keeps stable ids across filtering changes", () => {
        const rawSettings: AnyGuiSpec[] = [
            g({ type: "slider", label: "Opacity", current_value: 0.5 }),
            g({
                type: "folder",
                label: "Axes",
                current_value: [g({ type: "slider", label: "X min", current_value: 0 })],
            }),
        ];

        const { result, rerender } = renderHook(
            ({ term }) => useFilteredGroupedSettings(rawSettings, term),
            { initialProps: { term: "" } },
        );

        const axesBefore = result.current.find((entry) => entry.setting.label === "Axes");
        expect(axesBefore).toBeDefined();

        rerender({ term: "axes" });

        const axesAfter = result.current.find((entry) => entry.setting.label === "Axes");
        expect(axesAfter).toBeDefined();
        expect(axesAfter?.id).toBe(axesBefore?.id);
    });
});

