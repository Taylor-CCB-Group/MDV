import { describe, expect, test } from "vitest";

import { splitTooltipFields, tooltipItemsFromDataStore } from "@/react/spatialdata/table_association";

describe("splitTooltipFields", () => {
    test("fields in obs go to the viewer, the rest to MDV, in order", () => {
        expect(splitTooltipFields(["cluster", "CD3E", "area", "mdv_cell_id"], new Set(["cluster", "area"]))).toEqual({
            viewer: ["cluster", "area"],
            mdv: ["CD3E", "mdv_cell_id"],
        });
    });

    test("with no obs to go by, every field stays with the viewer", () => {
        expect(splitTooltipFields(["cluster"], undefined)).toEqual({ viewer: ["cluster"], mdv: [] });
        expect(splitTooltipFields(undefined, new Set())).toEqual({ viewer: [], mdv: [] });
    });
});

describe("tooltipItemsFromDataStore", () => {
    const column = (name: string, values: (string | number)[] | null) => ({
        name,
        data: values,
        getValue: (row: number) => values?.[row] ?? "",
    });

    test("reads each field's value at the feature's row, labelled by column name", () => {
        const dataStore = {
            columnIndex: {
                gene_CD3E: column("CD3E", [0.5, 1.5]),
                unloaded: column("Unloaded", null),
            },
        };
        const items = tooltipItemsFromDataStore(dataStore, 1, ["gene_CD3E", "unloaded", "missing"]);
        expect(items).toEqual([
            { label: "CD3E", value: "1.5" },
            { label: "Unloaded", value: "loading…" },
            { label: "missing", value: "" },
        ]);
    });
});
