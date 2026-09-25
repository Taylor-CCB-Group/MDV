import DataStore from "@/datastore/DataStore";
import { describe, expect, test } from "vitest";

describe("DataStore.getColorFunction", () => {
    test("keeps legacy single-slot multitext coloring", () => {
        const column = {
            field: "__tags__",
            name: "tags",
            datatype: "multitext" as const,
            delimiter: ";",
            stringLength: 1,
            values: ["a", "b", "a; b"],
            data: new Uint16Array([0, 1, 2]),
        };

        const store = {
            columnIndex: { __tags__: column },
            getColumnColors: () => ["#ff0000", "#00ff00", "#0000ff"],
        };

        const colorFn = DataStore.prototype.getColorFunction.call(
            store,
            "__tags__",
            {},
        );

        expect(colorFn(0)).toBe("#ff0000");
        expect(colorFn(1)).toBe("#00ff00");
        expect(colorFn(2)).toBe("#0000ff");
    });

    test("reads packed multitext colors using row stride", () => {
        const column = {
            field: "__gates__",
            name: "gates",
            datatype: "multitext" as const,
            delimiter: ",",
            stringLength: 3,
            values: ["N/A", "a", "b"],
            data: new Uint16Array([
                0, 65535, 65535,
                2, 65535, 65535,
                1, 2, 65535,
            ]),
        };

        const store = {
            columnIndex: { __gates__: column },
            getColumnColors: () => ["#999999", "#ff0000", "#00ff00"],
        };

        const colorFn = DataStore.prototype.getColorFunction.call(
            store,
            "__gates__",
            {},
        );

        expect(colorFn(0)).toBe("#999999");
        expect(colorFn(1)).toBe("#00ff00");
        expect(colorFn(2)).toBe("#ff0000");
    });
});
