import DataStore from "@/datastore/DataStore";
import RangeDimension from "@/datastore/RangeDimension";
import { useRangeDimension2D } from "@/react/hooks";
import { act, renderHook } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";

const mockUseDataStore = vi.fn();

vi.mock("@/react/context", () => ({
    useChart: () => ({ config: { param: ["x", "y"] } }),
    useDataStore: () => mockUseDataStore(),
}));

vi.mock("@/react/scatter_state", () => ({
    useRegionScale: () => 1,
}));

// x/y on a line: rows 0..4 at (i, i)
const SIZE = 5;

function createDataStore() {
    const ds = {
        size: SIZE,
        filterSize: SIZE,
        filterArray: new Uint8Array(SIZE),
        dimensions: [] as RangeDimension[],
        columnIndex: {
            x: { data: Float32Array.from({ length: SIZE }, (_, i) => i) },
            y: { data: Float32Array.from({ length: SIZE }, (_, i) => i) },
        },
        _callListeners: vi.fn(),
        getDimension(type: string) {
            expect(type).toBe("range_dimension");
            const dim = new RangeDimension(ds);
            ds.dimensions.push(dim);
            return dim;
        },
        removeAllFilters() {
            DataStore.prototype.removeAllFilters.call(ds);
        },
    };
    return ds;
}

const square = (lo: number, hi: number): [number, number][] => [
    [lo, lo],
    [hi, lo],
    [hi, hi],
    [lo, hi],
    [lo, lo],
];

describe("useRangeDimension2D", () => {
    beforeEach(() => {
        vi.stubGlobal(
            "Worker",
            vi.fn().mockImplementation(function MockWorker() {
                return { addEventListener: vi.fn(), postMessage: vi.fn(), terminate: vi.fn() };
            }),
        );
        vi.spyOn(console, "log").mockImplementation(() => {});
    });

    afterEach(() => {
        vi.unstubAllGlobals();
        vi.restoreAllMocks();
    });

    test("keeps a registered dimension after StrictMode's simulated unmount", () => {
        const ds = createDataStore();
        mockUseDataStore.mockReturnValue(ds);

        const { result, unmount } = renderHook(() => useRangeDimension2D(), { reactStrictMode: true });

        const dim = result.current.rangeDimension;
        expect(dim).toBeDefined();
        expect(ds.dimensions).toEqual([dim]);

        unmount();
        expect(ds.dimensions).toEqual([]);
    });

    test("a new filter after removeAllFilters only reflects the new polygon", () => {
        const ds = createDataStore();
        mockUseDataStore.mockReturnValue(ds);

        const { result } = renderHook(() => useRangeDimension2D(), { reactStrictMode: true });

        act(() => result.current.filterPoly(square(-0.5, 1.5)));
        expect(ds.filterSize).toBe(2);

        act(() => ds.removeAllFilters());
        expect(ds.filterSize).toBe(SIZE);

        act(() => result.current.filterPoly(square(2.5, 4.5)));
        expect(ds.filterSize).toBe(2);
        expect(Array.from(ds.filterArray)).toEqual([1, 1, 1, 0, 0]);
    });

    test("destroying a dimension twice does not remove another one", () => {
        const ds = createDataStore();
        const a = ds.getDimension("range_dimension");
        const b = ds.getDimension("range_dimension");
        a.destroy();
        a.destroy();
        expect(ds.dimensions).toEqual([b]);
    });
});
