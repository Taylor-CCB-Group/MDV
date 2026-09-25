import { type TableRowSource, getTableRowResolver, planTableRows } from "@/react/spatialdata/table_rows";
import { describe, expect, test } from "vitest";

type Row = Record<string, string | number>;

/** A DataStore stand-in: `getValue` per column, as `DataStore.columnIndex` has it. */
function rowSource(rows: Row[], loaded?: string[]): TableRowSource {
    const fields = Object.keys(rows[0] ?? {});
    return {
        size: rows.length,
        columnsWithData: loaded ?? fields,
        columnIndex: Object.fromEntries(
            fields.map((field) => [field, { getValue: (i: number) => rows[i]?.[field] ?? "missing" }]),
        ),
    };
}

// What the converter writes for a `by-region` group holding the same region from two
// SpatialData objects: `anndata.concat` in record order, provenance on every row.
const twoTableProvenance = {
    table_handling: "by-region",
    tables: [
        {
            spatialdata_name: "sample-a.zarr",
            table_name: "cells",
            table_id: "sample-a.zarr/cells",
            region: "cell_boundaries",
            region_key: "region",
            instance_key: "cell_id",
        },
        {
            spatialdata_name: "sample-b.zarr",
            table_name: "cells",
            table_id: "sample-b.zarr/cells",
            region: "cell_boundaries",
            region_key: "region",
            instance_key: "cell_id",
        },
    ],
};

const twoTableRows: Row[] = [
    { spatialdata_table_id: "sample-a.zarr/cells", region: "cell_boundaries", cell_id: "a0" },
    { spatialdata_table_id: "sample-a.zarr/cells", region: "cell_boundaries", cell_id: "a1" },
    { spatialdata_table_id: "sample-b.zarr/cells", region: "cell_boundaries", cell_id: "b0" },
    { spatialdata_table_id: "sample-b.zarr/cells", region: "cell_boundaries", cell_id: "b1" },
    { spatialdata_table_id: "sample-b.zarr/cells", region: "cell_boundaries", cell_id: "b2" },
];

describe("SpatialData table rows → datasource rows", () => {
    test("maps the second table of a merged datasource to its own rows, not the first table's", () => {
        const plan = planTableRows(twoTableProvenance, "cells", "sample-b.zarr");
        const resolver = getTableRowResolver(rowSource(twoTableRows), plan);

        // sd.js reports rows within sample-b's own zarr table: 0, 1, 2. Used as
        // DataStore rows directly they would colour a0, a1 and b0.
        const rows = resolver?.resolveRows("cell_boundaries", ["b1", "b0", "b2", "nope"], [1, 0, 2, 3]);

        expect(Array.from(rows ?? [])).toEqual([3, 2, 4, -1]);
    });

    test("joins on (region, instance id), never on id alone", () => {
        const provenance = {
            tables: [
                {
                    spatialdata_name: "s.zarr",
                    table_name: "t",
                    table_id: "s.zarr/t",
                    region: ["left", "right"],
                    region_key: "region",
                    instance_key: "id",
                },
                {
                    spatialdata_name: "s.zarr",
                    table_name: "other",
                    table_id: "s.zarr/other",
                    region: "elsewhere",
                    region_key: "region",
                    instance_key: "id",
                },
            ],
        };
        const rows: Row[] = [
            { spatialdata_table_id: "s.zarr/other", region: "elsewhere", id: 1 },
            { spatialdata_table_id: "s.zarr/t", region: "left", id: 1 },
            { spatialdata_table_id: "s.zarr/t", region: "right", id: 1 },
        ];
        const resolver = getTableRowResolver(rowSource(rows), planTableRows(provenance, "t", "s.zarr"));

        // Label values arrive as strings; the datasource stored them as numbers.
        expect(Array.from(resolver?.resolveRows("right", ["1"], [-1]) ?? [])).toEqual([2]);
        expect(Array.from(resolver?.resolveRows("left", ["1"], [-1]) ?? [])).toEqual([1]);
        expect(Array.from(resolver?.resolveRows("elsewhere", ["1"], [-1]) ?? [])).toEqual([-1]);
    });

    test("uses table rows positionally when the datasource is exactly one table", () => {
        const plan = planTableRows({ tables: [twoTableProvenance.tables[0]] }, "cells", "sample-a.zarr");
        expect(plan).toEqual({ kind: "positional" });

        const resolver = getTableRowResolver(rowSource(twoTableRows.slice(0, 2)), plan);
        expect(Array.from(resolver?.resolveRows("cell_boundaries", ["a1", "x", "y"], [1, -1, 7]) ?? [])).toEqual([
            1, -1, -1,
        ]);
    });

    test("uses table rows positionally for a datasource with no SpatialData provenance", () => {
        expect(planTableRows(undefined, "cells", undefined)).toEqual({ kind: "positional" });
    });

    test("waits for the key columns to load before resolving", () => {
        const plan = planTableRows(twoTableProvenance, "cells", "sample-b.zarr");
        expect(plan.kind === "value" && plan.columns).toEqual(["spatialdata_table_id", "cell_id"]);

        expect(getTableRowResolver(rowSource(twoTableRows, ["spatialdata_table_id"]), plan)).toBeUndefined();
    });

    test("refuses to guess which of several same-named tables an element belongs to", () => {
        // No SpatialData name to tell sample-a's `cells` from sample-b's.
        expect(planTableRows(twoTableProvenance, "cells", undefined).kind).toBe("unresolved");
        expect(planTableRows(twoTableProvenance, "cells", "sample-c.zarr").kind).toBe("unresolved");
    });

    test("builds the join index once per datasource and key columns", () => {
        const source = rowSource(twoTableRows);
        const plan = planTableRows(twoTableProvenance, "cells", "sample-b.zarr");
        expect(getTableRowResolver(source, plan)).toBe(getTableRowResolver(source, { ...plan }));
    });
});
