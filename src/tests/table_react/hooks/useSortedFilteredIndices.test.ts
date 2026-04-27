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

        // Unique column: 4 bytes each
        const namesArr = ["bb", "aa", "dd", "cc", "ee"];
        const namesBytes = new Uint8Array(namesArr.length * 2);
        const encoder = new TextEncoder();
        namesArr.forEach(
            (str, index) => namesBytes.set(encoder.encode(str), index * 2)
        );

        mockDataStore = {
            columnIndex: {
                age: { datatype: "integer" },
                name: {
                    datatype: "unique",
                    stringLength: 2,
                },
                tags: {
                    datatype: "multitext",
                    stringLength: 2,
                    values: ["N/A", "B", "C", "D", "E"],
                },
                cell_type: {
                    datatype: "text",
                    values: ["T cells", "B cells", "CD4 cells", "B2 cells", ""],
                },
            },
            getRawColumn: vi.fn((columnId: string) => {
                if (columnId === "age") {
                    return new Uint32Array([25, 30, 20, 35, 28]);
                }
                if (columnId === "name") {
                    return namesBytes;
                }
                if (columnId === "tags") {
                    // String length is 2
                    // Will be displayed as Row 0: "E, E", Row 1: "A, A"....
                    return new Uint16Array([
                        4, 4, 0, 0, 2, 2, 1, 1, 3, 3,
                    ]);
                }
                if (columnId === "cell_type") {
                    return new Uint8Array([0, 1, 2, 3, 4]);
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

    test("should sort numeric column ascending", async () => {
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

    test("should sort numeric column descending", async () => {
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

    test("should sort multitext column and put N/A at the end", async () => {
        const { result } = renderHook(() => useSortedFilteredIndices());

        act(() => {
            runInAction(() => {
                mockConfig.sort = { columnId: "tags", ascending: true };
            });
        });

        // Row indices - 0: "E, E", 1: "A, A", 2: "C, C", 3: "B, B", 4: "D, D"
        // Sorted indices - 1, 3, 2, 4, 0

        await waitFor(() => {
            expect(result.current).toEqual(new Uint32Array([3, 2, 4, 0, 1]));
        });
    });

    test("should sort text column and move empty value to the end", async () => {
        const { result } = renderHook(() => useSortedFilteredIndices());

        act(() => {
            runInAction(() => {
                mockConfig.sort = { columnId: "cell_type", ascending: true };
            });
        });

        await waitFor(() => {
            // Rows: 0: "T cells", 1: "B cells", 2: "CD4 cells", 3: "B2 cells", 4: ""
            // Sorted: B cells, B2 cells, CD4 cells, T cells, "" -> indices 1, 3, 2, 0, 4
            // Moving the empty value to the end
            expect(result.current).toEqual(new Uint32Array([1, 3, 2, 0, 4]));
        });
    });

    test("should sort unique column", async () => {
        const { result } = renderHook(() => useSortedFilteredIndices());

        act(() => {
            runInAction(() => {
                mockConfig.sort = { columnId: "name", ascending: true };
            });
        });

        // Row indices - 0: "bb", 1: "aa", 2: "dd", 3: "cc", 4: "ee"
        // Sorted indices - 1, 0, 3, 2, 4

        await waitFor(() => {
            expect(result.current).toEqual(new Uint32Array([1, 0, 3, 2, 4]));
        });
    });
});
