import { autorun, makeObservable, observable, runInAction } from "mobx";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type BaseChart from "@/charts/BaseChart";
import { ColumnQueryMapper, loadColumnData } from "@/datastore/decorateColumnMethod";
import type { FieldSpec } from "@/lib/columnTypeHelpers";
import { getFieldName, type RowsAsColslink, RowsAsColsQuery } from "@/links/link_utils";

/**
 * A `@loadColumnData` method called with an active link re-runs whenever the link's fields change: each run loads
 * the columns for the fields as they were then, and calls the method once they have loaded. Loads can finish in a
 * different order from the one they started in (e.g. as a view loads, ten unfiltered genes take longer to fetch than
 * the one gene left after the filter lands), and the chart must end up showing the link's current fields.
 */

type Row = RowsAsColslink["observableFields"][number];
function row(index: number, name: string): Row {
    return {
        index,
        name,
        fieldName: getFieldName("gs", name, index),
        fieldNameForSubgroup: (sg: string) => getFieldName(sg, name, index),
        // these tests are about which fields get applied, not the columns behind them
        get column(): never {
            throw new Error("column not available in this test");
        },
        columnForSubgroup: () => {
            throw new Error("column not available in this test");
        },
    };
}

/** column loads requested through `ChartManager._getColumnsThen`, finished by hand */
let loads: { columns: string[]; finish: () => void }[];
let link: RowsAsColslink;
/** what the decorated method was called with: a single-column method gets the first of the link's fields */
let applied: FieldSpec[];
let chart: any;
const setParam = loadColumnData(function setParam(this: BaseChart<any>, column: FieldSpec) {
    applied.push(column);
});

function setLinkRows(...rows: Row[]) {
    runInAction(() => {
        link.observableFields = rows;
    });
}

beforeEach(() => {
    vi.spyOn(console, "log").mockImplementation(() => {});
    loads = [];
    applied = [];
    vi.stubGlobal("mdv", {
        chartManager: {
            _getColumnsThen: (_ds: string, columns: string[], finish: () => void) => loads.push({ columns, finish }),
        },
    });
    link = makeObservable(
        {
            name_column: "gene",
            name: "genes",
            subgroups: { gs: { name: "gs", label: "Main", type: "double" } },
            observableFields: [],
            initPromise: Promise.resolve(),
        },
        { observableFields: observable.ref },
    );
    chart = {
        config: { type: "heat_map" },
        dataSource: { name: "cells" },
        reactionDisposers: [],
        mobxAutorun(fn: () => void) {
            const disposer = autorun(fn);
            this.reactionDisposers.push(disposer);
            return disposer;
        },
    };
    chart.activeQueries = new ColumnQueryMapper(chart, { setParam: "param" });
});
afterEach(() => {
    for (const dispose of chart.reactionDisposers) dispose();
    vi.restoreAllMocks();
    vi.unstubAllGlobals();
});

const A = row(0, "A");
const B = row(1, "B");
const TNF = row(96, "TNFRSF1B");

describe("@loadColumnData with an active link", () => {
    test("applies the fields once their columns have loaded", () => {
        setLinkRows(TNF);
        setParam.call(chart, new RowsAsColsQuery(link, "genes", 10));
        expect(applied).toEqual([]);
        loads[0].finish();
        expect(applied).toEqual([TNF.fieldName]);
    });

    test("a load for superseded fields that finishes last does not overwrite the current ones", () => {
        setLinkRows(A, B); // before the gene filter has been applied
        setParam.call(chart, new RowsAsColsQuery(link, "genes", 10));
        setLinkRows(TNF); // the filter lands while those columns are still loading
        expect(loads.map((l) => l.columns)).toEqual([[A.fieldName, B.fieldName], [TNF.fieldName]]);

        loads[1].finish();
        loads[0].finish();

        expect(applied.at(-1)).toBe(TNF.fieldName);
    });

    test("a load for superseded fields is not applied at all, even when it finishes first", () => {
        setLinkRows(A, B);
        setParam.call(chart, new RowsAsColsQuery(link, "genes", 10));
        setLinkRows(TNF);

        loads[0].finish();
        loads[1].finish();

        expect(applied).toEqual([TNF.fieldName]);
    });
});
