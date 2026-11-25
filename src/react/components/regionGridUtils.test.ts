import type { RenderStack, RenderStackEntry } from "@spatialdata/layers";
import { describe, expect, test } from "vitest";
import type { GuiSpec } from "@/charts/charts";
import type { SpatialDataGridConfig } from "@/charts/schemas/ChartConfigSchema";
import {
    getSpatialDataGridSettingsFolder,
    isRegionGridShown,
    moveRegionViews,
    spatialRegionCells,
} from "./regionGridUtils";

describe("spatialRegionCells", () => {
    test("a cell per region with a store and coordinate system", () => {
        const cells = spatialRegionCells({
            r1: { spatial: { file: "one.zarr", coordinate_system: "global" } },
            r2: { viv_image: { file: "x" } },
            r3: { spatial: { file: "three.zarr", coordinate_system: "global", url: "https://example.org/three.zarr" } },
        });
        expect(cells.map((cell) => cell.key)).toEqual(["r1", "r3"]);
        expect(cells[1]).toMatchObject({
            file: "three.zarr",
            coordinateSystem: "global",
            url: "https://example.org/three.zarr",
        });
    });
});

describe("isRegionGridShown", () => {
    test("only the SpatialData chart, with grid by region", () => {
        expect(isRegionGridShown({ type: "SpatialDataMdvRegionReact", grid: { by: "region" } })).toBe(true);
        expect(isRegionGridShown({ type: "SpatialDataMdvRegionReact", grid: { by: "densityFields" } })).toBe(false);
        expect(isRegionGridShown({ type: "DeckScatter", grid: { by: "region" } })).toBe(false);
    });
});

describe("moveRegionViews", () => {
    const a = { target: [0, 0] as [number, number], zoom: 0 };
    const b = { target: [100, 50] as [number, number], zoom: -2 };

    test("unlinked, only the cell that moved changes, and the others keep their identity", () => {
        const moved = moveRegionViews({ a, b }, "a", { target: [10, 0], zoom: 1 }, false);
        expect(moved.a).toEqual({ target: [10, 0], zoom: 1 });
        expect(moved.b).toBe(b);
    });

    test("linked, the others change zoom by as much and pan as far on screen", () => {
        // a pans 10 world units at zoom 1: 20 pixels
        const moved = moveRegionViews({ a, b }, "a", { target: [10, 0], zoom: 1 }, true);
        expect(moved.b.zoom).toBe(-1);
        // 20 pixels at zoom -1 is 40 world units
        expect(moved.b.target).toEqual([140, 50]);
    });

    test("linked, a zoom about a point zooms the others about the same point on screen", () => {
        // a zooms in by 1 about the point 8 world units right of its centre (8 pixels at zoom 0)
        const next = { target: [4, 0] as [number, number], zoom: 1 };
        const moved = moveRegionViews({ a, b }, "a", next, true);
        // b's point 8 pixels right of centre, at zoom -2, is 32 world units right of its centre
        const pointX = b.target[0] + 8 / 2 ** b.zoom;
        // after the zoom it is still 8 pixels right of b's centre
        expect((pointX - moved.b.target[0]) * 2 ** moved.b.zoom).toBeCloseTo(8);
    });

    test("cells without a view are left to be fitted", () => {
        const moved = moveRegionViews({ a }, "a", { target: [1, 1], zoom: 0 }, true);
        expect(Object.keys(moved)).toEqual(["a"]);
    });
});

describe("getSpatialDataGridSettingsFolder", () => {
    const layer = (id: string, elementType: "shapes" | "labels" | "image", elementKey: string): RenderStackEntry => ({
        kind: "spatial",
        id,
        visible: true,
        props: {},
        source: { elementType, elementKey },
    });
    const renderStack: RenderStack = {
        schemaVersion: 2,
        entries: [
            layer("he", "image", "he"),
            layer("circles", "shapes", "cell_circles"),
            layer("boundaries", "shapes", "cell_boundaries"),
        ],
    };
    const layerChoices = (folder: GuiSpec<"folder">) => {
        const [show] = folder.current_value;
        if (show?.type !== "dropdown") return [];
        const options: unknown[] = show.values?.[0] ?? [];
        return options.flatMap((option) =>
            option &&
            typeof option === "object" &&
            "value" in option &&
            typeof option.value === "string" &&
            option.value.startsWith("layer_fill:")
                ? [option.value]
                : [],
        );
    };
    const annotated = ({ elementKey }: { elementKey: string }) => elementKey === "cell_circles";

    test("offers a colour grid for each shapes or labels layer a table annotates", () => {
        const folder = getSpatialDataGridSettingsFolder({ renderStack }, annotated);
        expect(layerChoices(folder)).toEqual(["layer_fill:columns:circles", "layer_fill:categories:circles"]);
    });

    test("keeps the layer the chart is a grid of, and every layer when it cannot tell", () => {
        const grid: SpatialDataGridConfig = { by: "layer_fill", layer: "boundaries" };
        expect(layerChoices(getSpatialDataGridSettingsFolder({ renderStack, grid }, annotated))).toHaveLength(4);
        expect(layerChoices(getSpatialDataGridSettingsFolder({ renderStack }))).toHaveLength(4);
    });
});
