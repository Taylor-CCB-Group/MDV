import type { ElementIssue, EntryNotice } from "@spatialdata/core";
import type { EntryStatus } from "@spatialdata/vis";
import { afterEach, describe, expect, test } from "vitest";

import {
    dataQualityKey,
    elementIssueId,
    elementPath,
    isDataQualityDismissed,
    issuesByElement,
    projectSpatialStores,
    setDataQualityDismissed,
    withoutDismissed,
} from "@/react/spatialdata/data_quality";

const STORE = "https://example.org/store.zarr";
const noPyramid: ElementIssue = {
    kind: "no-pyramid",
    width: 4000,
    height: 3000,
    remedy: { op: "build-pyramid", factor: 2, coarsestAtMostPixels: 2048 },
};
const unit = (coordinateSystem: string): ElementIssue => ({
    kind: "no-length-unit",
    coordinateSystem,
    remedy: { op: "declare-length-unit", coordinateSystem },
});
const elementIssues = (issues: ElementIssue[]): EntryNotice => ({ kind: "element-issues", message: "m", issues });
const unmatched: EntryNotice = { kind: "table-rows-unmatched", message: "no rows", tableName: "t", instanceKey: "k" };

const touched: string[] = [];
function dismiss(path: string, id: string) {
    const key = dataQualityKey(STORE, path, id);
    touched.push(key);
    setDataQualityDismissed(key, true);
}

afterEach(() => {
    for (const key of touched.splice(0)) setDataQualityDismissed(key, false);
});

describe("withoutDismissed", () => {
    const status: EntryStatus = { errors: [], notices: [elementIssues([noPyramid, unit("global")]), unmatched] };

    test("takes out dismissed issues one by one, and other notices whole", () => {
        dismiss("images/he", elementIssueId(noPyramid));
        const { status: left, dismissed } = withoutDismissed(status, STORE, "images/he");
        expect(dismissed).toBe(1);
        expect(left?.notices).toEqual([elementIssues([unit("global")]), unmatched]);

        dismiss("images/he", "no-length-unit:global");
        dismiss("images/he", "table-rows-unmatched");
        expect(withoutDismissed(status, STORE, "images/he")).toEqual({ status: undefined, dismissed: 3 });
    });

    test("is keyed by store and element, and keeps errors", () => {
        dismiss("images/he", "table-rows-unmatched");
        expect(withoutDismissed(status, STORE, "images/other").dismissed).toBe(0);
        expect(withoutDismissed(status, "https://elsewhere/store.zarr", "images/he").dismissed).toBe(0);
        const error = {
            resource: "element",
            error: { kind: "element-not-found", message: "gone", retryable: false, elementKey: "he" },
        } satisfies EntryStatus["errors"][number];
        dismiss("images/he", elementIssueId(noPyramid));
        dismiss("images/he", "no-length-unit:global");
        const { status: left } = withoutDismissed({ errors: [error], notices: status.notices }, STORE, "images/he");
        expect(left).toEqual({ errors: [error], notices: [] });
    });

    test("persists in this browser", () => {
        dismiss("images/he", "table-rows-unmatched");
        expect(isDataQualityDismissed(dataQualityKey(STORE, "images/he", "table-rows-unmatched"))).toBe(true);
        expect(localStorage.getItem("mdv.dataQuality.dismissed")).toContain("table-rows-unmatched");
    });
});

describe("store issues", () => {
    test("names elements as sd.js does and groups issues by element", () => {
        expect(elementPath("image", "he")).toBe("images/he");
        expect(elementPath("labels", "cells")).toBe("labels/cells");
        const grouped = issuesByElement([
            { path: "images/he", issue: noPyramid },
            { path: "shapes/cells", issue: unit("global") },
            { path: "images/he", issue: unit("global") },
        ]);
        expect([...grouped.keys()]).toEqual(["images/he", "shapes/cells"]);
        expect(grouped.get("images/he")).toEqual([noPyramid, unit("global")]);
    });

    test("lists each store the project's regions point at once", () => {
        const region = (url: string) => ({ spatial: { file: "store.zarr", url } });
        const stores = projectSpatialStores([
            { name: "cells", dataStore: { regions: { all_regions: { a: region(STORE), b: region(STORE) } } } },
            { name: "genes", dataStore: { regions: { all_regions: { a: region(STORE), c: { spatial: {} } } } } },
            { name: "plain", dataStore: {} },
        ]);
        expect(stores).toEqual([{ url: STORE, dataSources: ["cells", "genes"] }]);
    });
});
