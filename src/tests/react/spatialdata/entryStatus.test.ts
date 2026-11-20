import type { ElementIssue, EntryNotice } from "@spatialdata/core";
import type { EntryStatus } from "@spatialdata/vis";
import { describe, expect, test } from "vitest";

import {
    describeRemedy,
    entryStatusErrorKey,
    entryStatusForDisplay,
    remediesFor,
} from "@/react/spatialdata/entry_status";

const noPyramid: ElementIssue = {
    kind: "no-pyramid",
    width: 4000,
    height: 3000,
    remedy: { op: "build-pyramid", factor: 2, coarsestAtMostPixels: 2048 },
};
const oversized: ElementIssue = {
    kind: "oversized-chunks",
    chunkWidth: 4000,
    chunkHeight: 4096,
    decodedBytes: 4000 * 4096 * 2,
    remedy: { op: "rechunk", chunkPixels: { min: 512, max: 1024 }, decodedBytesAtMost: 16 * 1024 * 1024 },
};
const unitIssue = (coordinateSystem: string): ElementIssue => ({
    kind: "no-length-unit",
    coordinateSystem,
    remedy: { op: "declare-length-unit", coordinateSystem },
});

function elementIssues(issues: ElementIssue[]): EntryNotice {
    return { kind: "element-issues", message: "issues", issues };
}

const truncated: EntryNotice = { kind: "preload-truncated", message: "Loaded 10 of 20", loaded: 10, total: 20 };

describe("describeRemedy", () => {
    test("words each remedy, once per distinct one", () => {
        expect(describeRemedy(noPyramid.remedy)).toContain("downsample by 2");
        expect(describeRemedy(oversized.remedy)).toContain("512–1024 px chunks (at most 16 MB");
        expect(describeRemedy(unitIssue("global").remedy)).toContain("'global'");
        expect(remediesFor([noPyramid, { ...noPyramid, kind: "no-pyramid" }, oversized])).toHaveLength(2);
    });
});

describe("entryStatusForDisplay", () => {
    test("keeps only the shown coordinate system's unit issue", () => {
        const status: EntryStatus = {
            errors: [],
            notices: [elementIssues([noPyramid, unitIssue("global"), unitIssue("aligned")])],
        };
        const shown = entryStatusForDisplay(status, "aligned", "image");
        const notice = shown?.notices[0];
        expect(notice?.kind === "element-issues" && notice.issues).toEqual([noPyramid, unitIssue("aligned")]);
        expect(notice?.message).not.toBe("issues");
    });

    test("drops a notice whose issues are all about other coordinate systems, and empty statuses", () => {
        const status: EntryStatus = { errors: [], notices: [elementIssues([unitIssue("global")])] };
        expect(entryStatusForDisplay(status, "aligned", "image")).toBeUndefined();
        expect(entryStatusForDisplay(undefined, "aligned")).toBeUndefined();
    });

    test("leaves out points' notices, which the points panel shows, but not their errors", () => {
        const error = {
            resource: "element",
            error: { kind: "element-not-found", message: "gone", retryable: false, elementKey: "p" },
        } satisfies EntryStatus["errors"][number];
        expect(entryStatusForDisplay({ errors: [], notices: [truncated] }, "global", "points")).toBeUndefined();
        expect(entryStatusForDisplay({ errors: [error], notices: [truncated] }, "global", "points")).toEqual({
            errors: [error],
            notices: [],
        });
        expect(entryStatusErrorKey("p", error)).toBe("p\0element\0element-not-found\0gone");
    });
});
