import type { ElementIssue, EntryNotice } from "@spatialdata/core";
import type { EntryStatus } from "@spatialdata/vis";
import { describe, expect, test } from "vitest";

import { issuesByElement, projectSpatialStores, splitStorageIssues } from "@/react/spatialdata/data_quality";

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

describe("splitStorageIssues", () => {
    test("takes an element's storage issues out of a layer's status, keeping the rest", () => {
        const status: EntryStatus = { errors: [], notices: [elementIssues([noPyramid, unit("global")]), unmatched] };
        expect(splitStorageIssues(status)).toEqual({
            status: { errors: [], notices: [unmatched] },
            storageIssues: [noPyramid, unit("global")],
        });
    });

    test("leaves nothing when storage issues were all there was, and keeps errors", () => {
        expect(splitStorageIssues({ errors: [], notices: [elementIssues([noPyramid])] })).toEqual({
            status: undefined,
            storageIssues: [noPyramid],
        });
        const error = {
            resource: "element",
            error: { kind: "element-not-found", message: "gone", retryable: false, elementKey: "he" },
        } satisfies EntryStatus["errors"][number];
        expect(splitStorageIssues({ errors: [error], notices: [elementIssues([noPyramid])] }).status).toEqual({
            errors: [error],
            notices: [],
        });
        const plain: EntryStatus = { errors: [], notices: [unmatched] };
        expect(splitStorageIssues(plain)).toEqual({ status: plain, storageIssues: [] });
        expect(splitStorageIssues(undefined)).toEqual({ status: undefined, storageIssues: [] });
    });
});

describe("store issues", () => {
    test("groups issues by element", () => {
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
