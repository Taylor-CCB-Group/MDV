import { describe, expect, test } from "vitest";
import { ColumnQueryMapper } from "@/datastore/decorateColumnMethod";
import { RowsAsColsQuery, type RowsAsColslink } from "@/links/link_utils";
import type BaseChart from "@/charts/BaseChart";

const link: RowsAsColslink = {
    name_column: "gene",
    name: "genes",
    subgroups: { gs: { name: "gs", label: "Main", type: "double" } },
    observableFields: [],
    initPromise: Promise.resolve(),
};

// #297 - settings UI must see the active link, not the concrete column it resolved to
describe("ColumnQueryMapper.activeValue", () => {
    const mapper = new ColumnQueryMapper({} as BaseChart<any>, { colorByColumn: "color_by" });

    test("returns the active RowsAsColsQuery rather than the concrete column", () => {
        const q = new RowsAsColsQuery(link, "genes", 1);
        mapper.userValues.colorByColumn = [q];
        expect(mapper.activeValue("colorByColumn", "gene1")).toBe(q);
    });

    test("returns the concrete value when the user value is a plain column", () => {
        mapper.userValues.colorByColumn = ["colA"];
        expect(mapper.activeValue("colorByColumn", "colA")).toBe("colA");
    });

    test("ignores a stale query once the config value has been cleared", () => {
        mapper.userValues.colorByColumn = [new RowsAsColsQuery(link, "genes", 1)];
        expect(mapper.activeValue("colorByColumn", undefined)).toBeUndefined();
    });
});
