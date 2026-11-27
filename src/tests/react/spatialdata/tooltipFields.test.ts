import type { LayerConfig } from "@spatialdata/vis";
import { describe, expect, test } from "vitest";
import {
    splitTooltipFields,
    tooltipItemsFromDataStore,
    withViewerReadableColumns,
} from "@/react/spatialdata/table_association";

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

describe("withViewerReadableColumns", () => {
    const layer: LayerConfig = {
        type: "shapes",
        id: "cells",
        elementKey: "cells",
        visible: true,
        opacity: 1,
        tooltipFields: ["cluster", "gs|CD3E (gs)|0"],
        fillColorByColumn: { columnName: "gs|CD3E (gs)|0", mode: "auto" },
    };

    test("keeps only the columns in obs, and none when obs is unknown", () => {
        expect(withViewerReadableColumns(layer, new Set(["cluster"]))).toEqual({
            type: "shapes",
            id: "cells",
            elementKey: "cells",
            visible: true,
            opacity: 1,
            tooltipFields: ["cluster"],
        });
        expect(withViewerReadableColumns(layer, undefined)).toMatchObject({ tooltipFields: [] });
        expect(withViewerReadableColumns(layer, undefined)).not.toHaveProperty("fillColorByColumn");
    });

    test("a layer the viewer can read all of is returned as is", () => {
        const readable = new Set(["cluster", "gs|CD3E (gs)|0"]);
        expect(withViewerReadableColumns(layer, readable)).toBe(layer);
    });
});
