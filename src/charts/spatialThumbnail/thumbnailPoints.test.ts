import { describe, expect, it } from "vitest";
import { layoutCellPositions } from "./thumbnailPoints";

describe("layoutCellPositions", () => {
    const center = { x: 64, y: -64 };

    it("puts the smaller y at the top of the thumbnail", () => {
        const positions = layoutCellPositions({
            x: [0, 10],
            y: [0, 20],
            regionCode: null,
            slotCodes: [],
            filter: null,
            centers: [center],
            size: 100,
        });
        expect(Array.from(positions)).toEqual([64 - 50, -64 + 50, 64 + 50, -64 - 50]);
    });

    it("drops filtered rows and rows from other regions", () => {
        const positions = layoutCellPositions({
            x: [0, 10, 5],
            y: [0, 10, 5],
            regionCode: [0, 1, 0],
            slotCodes: [0],
            filter: [0, 0, 1],
            centers: [center],
            size: 10,
        });
        expect(positions.length).toBe(2);
        expect(positions[0]).toBeCloseTo(center.x - 5);
        expect(positions[1]).toBeCloseTo(center.y + 5);
    });

    it("keeps each region's cells inside that thumbnail", () => {
        const positions = layoutCellPositions({
            x: [0, 10, 0, 10],
            y: [0, 0, 0, 0],
            regionCode: [3, 3, 7, 7],
            slotCodes: [3, 7],
            filter: null,
            centers: [
                { x: 0, y: 0 },
                { x: 200, y: 0 },
            ],
            size: 10,
        });
        expect(positions[0]).toBeCloseTo(-5);
        expect(positions[2]).toBeCloseTo(5);
        expect(positions[4]).toBeCloseTo(195);
        expect(positions[6]).toBeCloseTo(205);
    });
});
