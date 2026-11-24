import { describe, expect, test } from "vitest";
import { isRegionGridShown, partitionRowsByRegion, spatialRegionCells } from "./regionGridUtils";

describe("partitionRowsByRegion", () => {
    test("splits rows by region, in order, with their positions", () => {
        // row -> region code
        const codes = [0, 1, 0, 2, 1, 0];
        const values = ["a", "b", "c"];
        const rows = new Uint32Array([5, 1, 2, 4, 0]);
        const parts = partitionRowsByRegion(rows, codes, values);
        expect([...parts.keys()]).toEqual(["a", "b"]);
        expect([...(parts.get("a")?.rows ?? [])]).toEqual([5, 2, 0]);
        expect([...(parts.get("a")?.positions ?? [])]).toEqual([0, 2, 4]);
        expect([...(parts.get("b")?.rows ?? [])]).toEqual([1, 4]);
        expect([...(parts.get("b")?.positions ?? [])]).toEqual([1, 3]);
        expect(parts.has("c")).toBe(false);
    });
});

describe("spatialRegionCells", () => {
    test("a cell per region with a store and coordinate system", () => {
        const cells = spatialRegionCells({
            r1: { spatial: { file: "one.zarr", coordinate_system: "global" } },
            r2: { viv_image: { file: "x" } },
            r3: { spatial: { file: "three.zarr", coordinate_system: "global", url: "https://example.org/three.zarr" } },
        });
        expect(cells.map((cell) => cell.key)).toEqual(["r1", "r3"]);
        expect(cells[1]).toMatchObject({
            file: "three.zarr",
            coordinateSystem: "global",
            url: "https://example.org/three.zarr",
        });
    });
});

describe("isRegionGridShown", () => {
    test("only the SpatialData chart, with grid by region", () => {
        expect(isRegionGridShown({ type: "SpatialDataMdvRegionReact", grid: { by: "region" } })).toBe(true);
        expect(isRegionGridShown({ type: "SpatialDataMdvRegionReact", grid: { by: "densityFields" } })).toBe(false);
        expect(isRegionGridShown({ type: "DeckScatter", grid: { by: "region" } })).toBe(false);
    });
});
