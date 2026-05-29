import { describe, expect, test, vi } from "vitest";
import DataStore from "@/datastore/DataStore";
import { buildColorLegendSpec } from "./buildColorLegendSpec";

function createStore(
    columns: Record<string, Record<string, unknown>>,
): DataStore {
    return Object.assign(Object.create(DataStore.prototype), {
        columnIndex: columns,
        getColumnColors: vi.fn((column: string) => {
            const col = columns[column];
            if (col?.datatype === "double") {
                return ["#000000", "#ffffff"];
            }
            return ["#ff0000", "#00ff00"];
        }),
        getMinMaxForColumn: vi.fn(() => [0, 10]),
    });
}

describe("buildColorLegendSpec", () => {
    test("returns null for missing column", () => {
        const store = createStore({});
        expect(buildColorLegendSpec(store, "missing")).toBeNull();
    });

    test("builds categorical spec for text column", () => {
        const store = createStore({
            cell_type: {
                datatype: "text",
                name: "Cell type",
                values: ["A", "B"],
            },
        });
        const spec = buildColorLegendSpec(store, "cell_type");
        expect(spec).toEqual({
            kind: "categorical",
            label: "Cell type",
            column: "cell_type",
            items: [
                { color: "#ff0000", name: "A", value: "A" },
                { color: "#00ff00", name: "B", value: "B" },
            ],
        });
    });

    test("builds continuous spec for numeric column", () => {
        const store = createStore({
            expression: {
                datatype: "double",
                name: "Expression",
            },
        });
        const spec = buildColorLegendSpec(store, "expression");
        expect(spec).toMatchObject({
            kind: "continuous",
            label: "Expression",
            colors: ["#000000", "#ffffff"],
            range: [0, 10],
        });
    });

    test("applies override min and max", () => {
        const store = createStore({
            expression: {
                datatype: "double",
                name: "Expression",
            },
        });
        const spec = buildColorLegendSpec(store, "expression", {
            overideValues: { min: 2, max: 8 },
        });
        expect(spec).toMatchObject({
            kind: "continuous",
            range: [2, 8],
        });
    });

    test("uses config name override", () => {
        const store = createStore({
            cell_type: {
                datatype: "text",
                name: "Cell type",
                values: ["A"],
            },
        });
        const spec = buildColorLegendSpec(store, "cell_type", { name: "Custom" });
        expect(spec).toMatchObject({ label: "Custom" });
    });
});
