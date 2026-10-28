import { cpuEditablePicks } from "@/editable-layers/cpuPicks";
import { describe, expect, it } from "vitest";

// Screen = data * 2, so distances in screen pixels are twice those in data units.
const project = (p: readonly number[]) => [p[0] * 2, p[1] * 2];
const square = {
    type: "Feature" as const,
    properties: {},
    geometry: {
        type: "Polygon",
        coordinates: [
            [
                [10, 10],
                [20, 10],
                [20, 20],
                [10, 20],
                [10, 10],
            ],
        ],
    },
};
const handle = (x: number, y: number, positionIndex: number) => ({
    type: "Feature" as const,
    properties: {
        guideType: "editHandle",
        editHandleType: "existing",
        featureIndex: 0,
        positionIndexes: [0, positionIndex],
    },
    geometry: { type: "Point", coordinates: [x, y] },
});
const guides = [handle(10, 10, 0), handle(20, 10, 1), handle(20, 20, 2), handle(10, 20, 3)];

const picksAt = (x: number, y: number, radius = 4) =>
    cpuEditablePicks({ features: [square], guides, point: [x, y], project, radius });

describe("cpuEditablePicks", () => {
    it("picks the feature inside it", () => {
        expect(picksAt(30, 30).map((p) => [p.isGuide, p.index])).toEqual([[false, 0]]);
    });

    it("picks nothing outside it", () => {
        expect(picksAt(80, 80)).toEqual([]);
    });

    it("counts a point within the radius of an edge as on the feature", () => {
        expect(picksAt(30, 18).map((p) => p.index)).toEqual([0]);
        expect(picksAt(30, 14)).toEqual([]);
    });

    it("puts the nearest edit handle first, ahead of the feature", () => {
        const picks = picksAt(41, 41);
        expect(picks[0]).toMatchObject({ isGuide: true, index: 2 });
        expect(picks[0].object.properties).toMatchObject({ guideType: "editHandle", positionIndexes: [0, 2] });
        expect(picks.at(-1)).toMatchObject({ isGuide: false, index: 0 });
    });

    it("ignores guides that are not edit handles", () => {
        const line = {
            type: "Feature" as const,
            properties: { guideType: "tentative" },
            geometry: { type: "Point", coordinates: [20, 20] },
        };
        const picks = cpuEditablePicks({ features: [], guides: [line], point: [40, 40], project, radius: 4 });
        expect(picks).toEqual([]);
    });

    it("treats holes as outside", () => {
        const withHole = {
            ...square,
            geometry: {
                type: "Polygon",
                coordinates: [
                    square.geometry.coordinates[0],
                    [
                        [13, 13],
                        [17, 13],
                        [17, 17],
                        [13, 17],
                        [13, 13],
                    ],
                ],
            },
        };
        expect(cpuEditablePicks({ features: [withHole], guides: [], point: [30, 30], project, radius: 1 })).toEqual([]);
        expect(
            cpuEditablePicks({ features: [withHole], guides: [], point: [23, 30], project, radius: 1 }),
        ).toHaveLength(1);
    });
});
