import { describe, expect, test } from "vitest";
import {
    COLOR_GRID_CHOICES,
    colorByField,
    colorGridCells,
    colorGridChoice,
    colorGridChoices,
    colorGridConfig,
    colorGridSplit,
    colorGridSplitFor,
    gridChoiceValue,
    partitionRowsByCategory,
    supportsColorGridMode,
} from "./colorGridUtils";

describe("partitionRowsByCategory", () => {
    test("splits rows by category, in order, with their positions", () => {
        // row -> category code
        const codes = [0, 1, 0, 2, 1, 0];
        const values = ["a", "b", "c"];
        const rows = new Uint32Array([5, 1, 2, 4, 0]);
        const parts = partitionRowsByCategory(rows, codes, values);
        expect([...parts.keys()]).toEqual(["a", "b"]);
        expect([...(parts.get("a")?.rows ?? [])]).toEqual([5, 2, 0]);
        expect([...(parts.get("a")?.positions ?? [])]).toEqual([0, 2, 4]);
        expect([...(parts.get("b")?.rows ?? [])]).toEqual([1, 4]);
        expect([...(parts.get("b")?.positions ?? [])]).toEqual([1, 3]);
        expect(parts.has("c")).toBe(false);
    });

    test("leaves out rows whose code has no value", () => {
        const parts = partitionRowsByCategory([0, 1], [0, 7], ["a"]);
        expect([...(parts.get("a")?.rows ?? [])]).toEqual([0]);
        expect(parts.size).toBe(1);
    });
});

describe("colorGridSplit", () => {
    test("any 2D deck scatter plot, and the SpatialData chart, with a grid by color_fields", () => {
        const config = { type: "DeckScatter", grid: { by: "color_fields" } } as const;
        expect(colorGridSplit(config)).toBe("columns");
        expect(colorGridSplit({ ...config, grid: { by: "color_fields", split: "categories" } })).toBe("categories");
        expect(colorGridSplit({ ...config, type: "DeckContourScatter" })).toBe("columns");
        expect(colorGridSplit({ ...config, type: "SpatialDataMdvRegionReact" })).toBe("columns");
        expect(colorGridSplit({ ...config, dimension: "3d" })).toBeNull();
        expect(colorGridSplit({ ...config, grid: { by: "densityFields" } })).toBeNull();
        expect(colorGridSplit({ ...config, grid: null })).toBeNull();
        expect(supportsColorGridMode("DeckScatter3D")).toBe(false);
        expect(supportsColorGridMode(undefined)).toBe(false);
    });
});

describe("grid settings choices", () => {
    test("name each kind of grid, and a colour grid's split", () => {
        expect(gridChoiceValue({ grid: null })).toBe("none");
        expect(gridChoiceValue({ grid: { by: "densityFields" } })).toBe("densityFields");
        expect(gridChoiceValue({ grid: { by: "color_fields" } })).toBe("color_fields:columns");
        expect(gridChoiceValue({ grid: { by: "color_fields", split: "categories" } })).toBe("color_fields:categories");
        for (const choice of COLOR_GRID_CHOICES) {
            const split = colorGridChoice(choice.value);
            expect(split).not.toBeNull();
            if (split) expect(gridChoiceValue({ grid: colorGridConfig(split) })).toBe(choice.value);
        }
        expect(colorGridChoice("region")).toBeNull();
    });

    test("offer the split the colour column's type allows, and the current one", () => {
        const values = (choices: { value: string }[]) => choices.map((choice) => choice.value);
        expect(values(colorGridChoices("columns", null))).toEqual(["color_fields:columns"]);
        expect(values(colorGridChoices("categories", null))).toEqual(["color_fields:categories"]);
        expect(values(colorGridChoices(null, null))).toEqual([]);
        expect(values(colorGridChoices("columns", "categories"))).toEqual([
            "color_fields:columns",
            "color_fields:categories",
        ]);
    });
});

describe("colorGridSplitFor", () => {
    test("by category for one category per row, by column for numbers, else none", () => {
        expect(colorGridSplitFor({ datatype: "text" })).toBe("categories");
        expect(colorGridSplitFor({ datatype: "text16" })).toBe("categories");
        expect(colorGridSplitFor({ datatype: "double" })).toBe("columns");
        expect(colorGridSplitFor({ datatype: "integer" })).toBe("columns");
        expect(colorGridSplitFor({ datatype: "int32" })).toBe("columns");
        expect(colorGridSplitFor({ datatype: "multitext" })).toBeNull();
        expect(colorGridSplitFor({ datatype: "unique" })).toBeNull();
        expect(colorGridSplitFor(undefined)).toBeNull();
    });
});

describe("colorByField", () => {
    test("a field, or the field of a column", () => {
        expect(colorByField("cluster")).toBe("cluster");
        expect(colorByField({ column: { field: "cluster" } })).toBe("cluster");
        expect(colorByField("")).toBeUndefined();
        expect(colorByField(null)).toBeUndefined();
        expect(colorByField({ fields: ["a"] })).toBeUndefined();
    });
});

describe("colorGridCells", () => {
    test("a cell per category, in its colour", () => {
        const cells = colorGridCells(["a", "b"], [[255, 0, 0], "#00ff00"]);
        expect(cells).toEqual([
            { key: "a", label: "a", color: [255, 0, 0] },
            { key: "b", label: "b", color: [226, 232, 240] },
        ]);
    });
});
