import type DataStore from "@/datastore/DataStore";
import { fillColorSchemeFromDataStore } from "@/react/spatialdata/fill_color_scheme";
import { describe, expect, test } from "vitest";

/**
 * The contract here is parity: a column handed to the viewer must come back in the
 * colours MDV would have drawn it in itself. If it does not, the same `cell_type`
 * reads one way on the spatial canvas and another on every other chart over the
 * same data — which looks like a data problem, not a rendering one.
 */

type StubColumn = {
    datatype: string;
    values?: string[];
    minMax?: [number, number];
};

function stubDataStore(columns: Record<string, StubColumn>, colorsByColumn: Record<string, number[][]>): DataStore {
    return {
        columnIndex: columns,
        getColumnColors: (columnName: string) => colorsByColumn[columnName],
    } as unknown as DataStore;
}

describe("fillColorSchemeFromDataStore", () => {
    test("names each category's colour, so the viewer cannot reorder them", () => {
        const store = stubDataStore(
            { Leiden: { datatype: "text", values: ["stroma", "tumour"] } },
            {
                Leiden: [
                    [10, 20, 30],
                    [40, 50, 60],
                ],
            },
        );

        expect(fillColorSchemeFromDataStore(store, "Leiden")).toEqual({
            columnName: "Leiden",
            mode: "categorical",
            categoricalPalette: {
                byValue: { stroma: [10, 20, 30], tumour: [40, 50, 60] },
            },
        });
    });

    test("keys an integer-coded category by its string form", () => {
        // The viewer normalises every cell to a string before matching, so a
        // cluster spelled `3` has to be `"3"` in the palette or it goes unnamed.
        const store = stubDataStore(
            { cluster: { datatype: "text", values: ["3", "10"] } },
            {
                cluster: [
                    [1, 1, 1],
                    [2, 2, 2],
                ],
            },
        );

        expect(fillColorSchemeFromDataStore(store, "cluster")).toMatchObject({
            categoricalPalette: { byValue: { "3": [1, 1, 1], "10": [2, 2, 2] } },
        });
    });

    test("pins a numeric column to the column's own range, not the view's", () => {
        const store = stubDataStore(
            { area: { datatype: "double", minMax: [0, 500] } },
            {
                area: [
                    [0, 0, 255],
                    [255, 255, 255],
                    [255, 0, 0],
                ],
            },
        );

        expect(fillColorSchemeFromDataStore(store, "area")).toEqual({
            columnName: "area",
            mode: "continuous",
            numericRamp: [
                [0, 0, 255],
                [255, 255, 255],
                [255, 0, 0],
            ],
            numericDomain: [0, 500],
        });
    });

    test("never sets numericScale, because the log remap is already in the stops", () => {
        // `getColumnColors` returns bins that have already been remapped when the
        // column is on a log scale. Asking the viewer to apply symlog as well would
        // apply it twice, which is a subtler bug than not applying it at all.
        const store = stubDataStore(
            { counts: { datatype: "integer", minMax: [0, 1000] } },
            {
                counts: [
                    [0, 0, 0],
                    [9, 9, 9],
                ],
            },
        );

        expect(fillColorSchemeFromDataStore(store, "counts")).not.toHaveProperty("numericScale");
    });

    test("declines rather than guessing when MDV cannot answer yet", () => {
        const noMinMax = stubDataStore({ area: { datatype: "double" } }, { area: [[0, 0, 0]] });
        const noValues = stubDataStore({ label: { datatype: "text", values: [] } }, {});
        const unknown = stubDataStore({}, {});

        // A numeric column's minMax is only computed once its data lands. Until
        // then the viewer measuring the loaded extent is a better answer than a
        // domain MDV made up.
        expect(fillColorSchemeFromDataStore(noMinMax, "area")).toBeUndefined();
        expect(fillColorSchemeFromDataStore(noValues, "label")).toBeUndefined();
        expect(fillColorSchemeFromDataStore(unknown, "missing")).toBeUndefined();
    });

    test("passes over a datatype it has no ramp or palette for", () => {
        const store = stubDataStore({ barcode: { datatype: "unique" } }, {});
        expect(fillColorSchemeFromDataStore(store, "barcode")).toBeUndefined();
    });
});
