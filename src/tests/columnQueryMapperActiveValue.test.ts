import { $mobx, type IReactionDisposer } from "mobx";
import { describe, expect, test, vi } from "vitest";
import BaseChart from "@/charts/BaseChart";
import { ColumnQueryMapper } from "@/datastore/decorateColumnMethod";
import { type RowsAsColslink, RowsAsColsQuery } from "@/links/link_utils";

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
        expect(mapper.activeValue("colorByColumn", "colB")).toBe("colB");
    });

    test("keeps an unresolved query distinct from a cleared setting", () => {
        const q = new RowsAsColsQuery(link, "genes", 1);
        mapper.userValues.setToolTipColumn = [q];
        // query resolves to no fields -> concrete config value is undefined, but it is still active
        expect(mapper.activeValue("setToolTipColumn", undefined)).toBe(q);
    });

    test("clearActiveQuery makes the setting count as cleared", () => {
        const spy = vi.fn();
        const dispose: IReactionDisposer = Object.assign(() => spy(), {
            [$mobx]: {} as IReactionDisposer[typeof $mobx],
        });
        const chart = { reactionDisposers: [dispose] } as unknown as BaseChart<any>;
        const m = new ColumnQueryMapper(chart, { setToolTipColumn: "tooltip.column" });
        m.userValues.setToolTipColumn = [new RowsAsColsQuery(link, "genes", 1)];
        m.reactionDisposers.set("setToolTipColumn", dispose);
        m.clearActiveQuery("setToolTipColumn");
        expect(spy).toHaveBeenCalledOnce();
        expect(chart.reactionDisposers).toEqual([]);
        expect(m.activeValue("setToolTipColumn", undefined)).toBeUndefined();
    });

    test("clearActiveQuery invalidates in-flight column loads", () => {
        const m = new ColumnQueryMapper({ reactionDisposers: [] } as unknown as BaseChart<any>, {});
        const before = m.generation("colorByColumn");
        m.clearActiveQuery("colorByColumn");
        expect(m.generation("colorByColumn")).not.toBe(before);
    });
});

// the colour toggles (overlay, log scale, etc.) all re-apply colour-by through this helper
describe("BaseChart.reapplyColorBy", () => {
    function stubChart(userValue: unknown[] | undefined, color_by: string | undefined) {
        const mapper = new ColumnQueryMapper({} as BaseChart<any>, { colorByColumn: "color_by" });
        if (userValue) mapper.userValues.colorByColumn = userValue;
        const colorByColumn = vi.fn();
        const chart = { config: { color_by }, activeQueries: mapper, colorByColumn };
        return { chart, colorByColumn };
    }

    test("passes the active query, not the concrete column", () => {
        const q = new RowsAsColsQuery(link, "genes", 1);
        const { chart, colorByColumn } = stubChart([q], "gene1");
        BaseChart.prototype.reapplyColorBy.call(chart as unknown as BaseChart<any>);
        expect(colorByColumn).toHaveBeenCalledWith(q);
    });

    test("passes the concrete column when there is no active query", () => {
        const { chart, colorByColumn } = stubChart(["colA"], "colA");
        BaseChart.prototype.reapplyColorBy.call(chart as unknown as BaseChart<any>);
        expect(colorByColumn).toHaveBeenCalledWith("colA");
    });
});
