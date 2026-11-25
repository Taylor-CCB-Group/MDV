import { describe, expect, test } from "vitest";
import { colorLegendField, migrateColorBy } from "./colorFields";

describe("migrateColorBy", () => {
    test("a saved color_by becomes the one colour column, and is dropped", () => {
        const config = migrateColorBy({ type: "DeckScatter", color_by: "cluster" });
        expect(config.color_fields).toEqual(["cluster"]);
        expect("color_by" in config).toBe(false);
    });

    test("the legacy object form names its column's field", () => {
        expect(migrateColorBy({ color_by: { column: { field: "cluster" } } }).color_fields).toEqual(["cluster"]);
    });

    test("an active link in color_by stays a link, as saved or live", () => {
        const saved = { type: "RowsAsColsQuery", linkedDsName: "genes", subgroupName: "scores", maxItems: 10 };
        expect(migrateColorBy({ color_by: saved }).color_fields).toEqual([saved]);
        const live = { fields: ["CD3E"], columns: [] };
        expect(migrateColorBy({ color_by: live }).color_fields).toEqual([live]);
    });

    test("no colour, no columns; color_fields already set wins", () => {
        expect(migrateColorBy({}).color_fields).toEqual([]);
        const config = migrateColorBy({ color_by: "old", color_fields: ["a", "b"] });
        expect(config.color_fields).toEqual(["a", "b"]);
        expect("color_by" in config).toBe(false);
    });
});

describe("colorLegendField", () => {
    test("the first colour column, unless each cell has a column of its own", () => {
        const config = { type: "DeckScatter", color_fields: ["a", "b"] };
        expect(colorLegendField(config)).toBe("a");
        expect(colorLegendField({ ...config, grid: { by: "color_fields", split: "categories" } })).toBe("a");
        expect(colorLegendField({ ...config, grid: { by: "color_fields", split: "columns" } })).toBeUndefined();
        expect(colorLegendField({ ...config, color_fields: [] })).toBeUndefined();
    });
});
