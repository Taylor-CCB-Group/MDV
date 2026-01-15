import useSortedFilteredIndices from "@/react/hooks/useSortedFilteredIndices";
import { renderHook, waitFor, act } from "@testing-library/react";
import { describe, test, expect, beforeEach, vi } from "vitest";
import { observable, runInAction } from "mobx";
import type { TableChartReactConfig } from "@/react/components/TableChartReactWrapper";

// Mock the context hooks
vi.mock("@/react/context", () => ({
    useDataStore: () => mockDataStore,
}));

vi.mock("@/react/hooks", () => ({
    useConfig: () => mockConfig,
    useSimplerFilteredIndices: () => mockFilteredIndices,
}));

let mockDataStore: any;
let mockConfig: TableChartReactConfig;
let mockFilteredIndices: Uint32Array;

describe("useSortedFilteredIndices", () => {
    beforeEach(() => {
        vi.clearAllMocks();

        mockFilteredIndices = new Uint32Array([0, 1, 2, 3, 4]);

        mockConfig = observable({
            sort: undefined as { columnId: string; ascending: boolean } | undefined,
        }) as any;

        mockDataStore = {
            columnIndex: {
                age: { datatype: "integer" },
            },
            getRawColumn: vi.fn((columnId: string) => {
                if (columnId === "age") {
                    return new Uint32Array([25, 30, 20, 35, 28]);
                }
                return null;
            }),
        };
    });

    test("should return unsorted indices by default", async () => {
        const { result } = renderHook(() => useSortedFilteredIndices());

        await waitFor(() => {
            expect(result.current).toEqual(new Uint32Array([0, 1, 2, 3, 4]));
        });
    });

    test("should sort ascending", async () => {
        const { result } = renderHook(() => useSortedFilteredIndices());

        act(() => {
            runInAction(() => {
                mockConfig.sort = { columnId: "age", ascending: true };
            });
        });

        await waitFor(() => {
            // [25, 30, 20, 35, 28] -> [20, 25, 28, 30, 35]
            expect(result.current).toEqual(new Uint32Array([2, 0, 4, 1, 3]));
        });
    });

    test("should sort descending", async () => {
        const { result } = renderHook(() => useSortedFilteredIndices());

        act(() => {
            runInAction(() => {
                mockConfig.sort = { columnId: "age", ascending: false };
            });
        });

        await waitFor(() => {
            // [25, 30, 20, 35, 28] -> [35, 30, 28, 25, 20]
            expect(result.current).toEqual(new Uint32Array([3, 1, 4, 0, 2]));
        });
    });
});
