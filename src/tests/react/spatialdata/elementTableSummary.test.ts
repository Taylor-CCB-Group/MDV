import { describe, expect, test } from "vitest";
import { elementTableDescription, summarizeElementTable } from "@/react/components/spatialLayers/ElementTableBadge";
import type { DataSourceAssociationCandidate } from "@/react/spatialdata/table_association";

function obsDataSource(name: string, tableName: string): DataSourceAssociationCandidate {
    return {
        name,
        dataStore: {
            config: {
                columns: [{ field: "spatialdata_table_id" }],
                spatialdata_tables: { tables: [{ table_name: tableName }] },
            },
        },
    };
}

/** A store where `cell_circles` is annotated by the given tables and nothing else is. */
function storeAnnotating(tableNames: string[]) {
    return {
        getAssociatedTables: (kind: string, key: string): Array<[string, object]> =>
            kind === "shapes" && key === "cell_circles" ? tableNames.map((name) => [name, {}]) : [],
    };
}

describe("summarizeElementTable", () => {
    test("names the annotating datasource", () => {
        const summary = summarizeElementTable(
            storeAnnotating(["table"]),
            [obsDataSource("cells", "table")],
            "shapes",
            "cell_circles",
        );
        expect(summary).toEqual({ status: "resolved", tableName: "table", dataSourceName: "cells" });
    });

    test("reports an element no table annotates", () => {
        const summary = summarizeElementTable(
            storeAnnotating(["table"]),
            [obsDataSource("cells", "table")],
            "shapes",
            "cell_boundaries",
        );
        expect(summary).toEqual({ status: "none" });
    });

    test("flags several tables on one element as ambiguous", () => {
        const summary = summarizeElementTable(storeAnnotating(["a", "b"]), [], "shapes", "cell_circles");
        expect(summary).toEqual({ status: "ambiguous", detail: "several tables annotate it (a, b)" });
    });

    test("flags a table that matches several datasources as ambiguous", () => {
        const summary = summarizeElementTable(
            storeAnnotating(["table"]),
            [obsDataSource("cells", "table"), obsDataSource("cells_copy", "table")],
            "shapes",
            "cell_circles",
        );
        expect(summary?.status).toBe("ambiguous");
        if (summary) expect(elementTableDescription(summary)).toContain("cells, cells_copy");
    });

    test("says nothing for images and points, or before the store is read", () => {
        const store = storeAnnotating(["table"]);
        expect(summarizeElementTable(store, [], "image", "cell_circles")).toBeNull();
        expect(summarizeElementTable(store, [], "points", "cell_circles")).toBeNull();
        expect(summarizeElementTable(undefined, [], "shapes", "cell_circles")).toBeNull();
    });
});
