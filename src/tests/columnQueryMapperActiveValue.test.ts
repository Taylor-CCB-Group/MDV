import { describe, expect, test, vi } from "vitest";
import { ColumnQueryMapper } from "@/datastore/decorateColumnMethod";
import { RowsAsColsQuery, type RowsAsColslink } from "@/links/link_utils";
import { $mobx, type IReactionDisposer } from "mobx";
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

    test("keeps an unresolved query distinct from a cleared setting", () => {
        const q = new RowsAsColsQuery(link, "genes", 1);
        mapper.userValues.setToolTipColumn = [q];
        // query resolves to no fields -> concrete config value is undefined, but it is still active
        expect(mapper.activeValue("setToolTipColumn", undefined)).toBe(q);
    });

    test("clearActiveQuery makes the setting count as cleared", () => {
        const spy = vi.fn();
        const dispose: IReactionDisposer = Object.assign(() => spy(), { [$mobx]: {} as IReactionDisposer[typeof $mobx] });
        const chart = { reactionDisposers: [dispose] } as unknown as BaseChart<any>;
        const m = new ColumnQueryMapper(chart, { setToolTipColumn: "tooltip.column" });
        m.userValues.setToolTipColumn = [new RowsAsColsQuery(link, "genes", 1)];
        m.reactionDisposers.set("setToolTipColumn", dispose);
        m.clearActiveQuery("setToolTipColumn");
        expect(spy).toHaveBeenCalledOnce();
        expect(chart.reactionDisposers).toEqual([]);
        expect(m.activeValue("setToolTipColumn", undefined)).toBeUndefined();
    });
});
