import type { DataType, LoadedDataColumn } from "@/charts/charts";
import type DataStore from "@/datastore/DataStore";
import type { FeedbackAlert } from "@/react/components/TableChartReactComponent";
import useFindReplace from "@/react/hooks/useFindReplace";
import { renderHook, act } from "@testing-library/react";
import { describe, test, expect, beforeEach, vi } from "vitest";

describe("useFindReplace", () => {
    let orderedParamColumns: LoadedDataColumn<DataType>[];
    let sortedIndices: Uint32Array;
    let dataStore: DataStore;
    let searchColumn: string | null;
    let config: { include_index?: boolean };
    let gridRef: React.MutableRefObject<any>;
    let selectionSourceRef: React.MutableRefObject<'user' | 'programmatic' | null>;
    let setFeedbackAlert: ReturnType<typeof vi.fn>;
    let mockGrid: any;

    beforeEach(() => {
        // Reset all mocks before each test
        vi.clearAllMocks();

        // Setup mock grid
        mockGrid = {
            gotoCell: vi.fn(),
            invalidate: vi.fn(),
            render: vi.fn(),
        };

        // Setup basic columns
        orderedParamColumns = [
            {
                field: "sample_count",
                datatype: "integer",
                data: new Uint32Array([10, 20, 30]),
                editable: true,
                getValue: vi.fn((dataIndex: number) => {
                    return orderedParamColumns[0].data[dataIndex];
                })
            },
            {
                field: "cell_type",
                datatype: "text",
                data: new Uint8Array([0, 1, 2]),
                values: ["T cells", "B cells", "NK cells"],
                editable: true,
                getValue: vi.fn((dataIndex: number) => {
                    const column = orderedParamColumns[1];
                    const valueIndex = column.data[dataIndex];
                    return column.values[valueIndex] || null;
                })
            },
        ] as any;

        sortedIndices = new Uint32Array([0, 1, 2]);
        searchColumn = "cell_type";
        config = { include_index: false };

        // Mock DataStore
        dataStore = {
            dataChanged: vi.fn(),
        } as any;

        gridRef = { current: { slickGrid: mockGrid } } as any;
        selectionSourceRef = { current: null };
        setFeedbackAlert = vi.fn();
    });

    describe("handleFind", () => {
        test("should find matches in text column", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("cells");
            });

            expect(result.current.foundMatches).toHaveLength(3);
            expect(result.current.matchCount).toBe(3);
            expect(result.current.currentMatchIndex).toBe(0);
            expect(mockGrid.gotoCell).toHaveBeenCalledWith(0, 1, false);
        });

        test("should find matches case-insensitively", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );
        
            act(() => {
                result.current.handleFind("B CELLS");  // uppercase search
            });
        
            expect(result.current.foundMatches).toHaveLength(1);
            expect(result.current.foundMatches[0].value).toBe("B cells");
        });

        test("should find no matches for non-existent text", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("nonexistent");
            });

            expect(result.current.foundMatches).toHaveLength(0);
            expect(result.current.matchCount).toBe(0);
            expect(result.current.currentMatchIndex).toBe(-1);
        });

        test("should clear matches when find text is empty", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("");
            });

            expect(result.current.foundMatches).toHaveLength(0);
            expect(result.current.matchCount).toBeNull();
            expect(result.current.currentMatchIndex).toBe(-1);
        });

        test("should show error feedback if column not found", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    "invalid_column",
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("test");
            });

            expect(setFeedbackAlert).toHaveBeenCalledWith(
                expect.objectContaining({
                    type: "error",
                    title: "Find Error",
                    message: expect.stringContaining("Column index not found"),
                }),
            );
        });

        test("should find matches in numeric column", () => {
            const testSearchColumn = "sample_count";
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    testSearchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("20");
            });

            expect(result.current.foundMatches).toHaveLength(1);
            expect(result.current.foundMatches[0].value).toBe("20");
        });
    });

    describe("handleFindNext", () => {
        test("should navigate to next match", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("cells");
            });

            expect(result.current.currentMatchIndex).toBe(0);

            act(() => {
                result.current.handleFindNext();
            });

            expect(result.current.currentMatchIndex).toBe(1);
            expect(mockGrid.gotoCell).toHaveBeenCalledTimes(2); // Once for find, once for findNext
        });

        test("should not navigate past last match", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("B cells");
            });

            expect(result.current.currentMatchIndex).toBe(0);

            act(() => {
                result.current.handleFindNext();
            });

            // Should still be at 0 since there's only 1 match
            expect(result.current.currentMatchIndex).toBe(0);
        });
    });

    describe("handleFindPrev", () => {
        test("should navigate to previous match", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("cells");
            });

            act(() => {
                result.current.handleFindNext();
            });

            expect(result.current.currentMatchIndex).toBe(1);

            act(() => {
                result.current.handleFindPrev();
            });

            expect(result.current.currentMatchIndex).toBe(0);
        });

        test("should not navigate before first match", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("cells");
            });

            expect(result.current.currentMatchIndex).toBe(0);

            act(() => {
                result.current.handleFindPrev();
            });

            // Should still be at 0
            expect(result.current.currentMatchIndex).toBe(0);
        });
    });

    describe("handleReplace", () => {
        test("should replace current match successfully", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("T cells");
            });

            act(() => {
                result.current.handleReplace("T cells", "T lymphocytes");
            });

            expect(setFeedbackAlert).toHaveBeenCalledWith(
                expect.objectContaining({
                    type: "success",
                    title: "Replace Successful",
                    message: expect.stringContaining(`Replaced T cells with T lymphocytes in column: ${searchColumn}`)
                })
            );
        });

        test("should show error if column not editable", () => {
            orderedParamColumns[1].editable = false;

            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("cells");
                result.current.handleReplace("cells", "new");
            });

            expect(setFeedbackAlert).toHaveBeenCalledWith(
                expect.objectContaining({
                    type: "error",
                    title: "Replace Error",
                    message: expect.stringContaining("not editable"),
                }),
            );
        });

        test("should show error if no match found", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                // Don't call handleFind first
                result.current.handleReplace("find", "replace");
            });

            expect(setFeedbackAlert).toHaveBeenCalledWith(
                expect.objectContaining({
                    type: "error",
                    title: "Replace Error",
                    message: expect.stringContaining("No current match"),
                }),
            );
        });
    });

    describe("handleReplaceAll", () => {
        test("should replace all matches successfully", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("cells");
            });

            act(() => {
                result.current.handleReplaceAll("cells", "lymphocytes");
            });

            expect(setFeedbackAlert).toHaveBeenCalledWith(
                expect.objectContaining({
                    type: "success",
                    title: "Replace All Successful",
                    message: expect.stringContaining(`Replaced 3 occurrences of cells with lymphocytes in column: ${searchColumn}`)
                })
            );
        });

        test("should show error if no matches found", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleReplaceAll("nonexistent", "new");
            });

            expect(setFeedbackAlert).toHaveBeenCalledWith(
                expect.objectContaining({
                    type: "error",
                    title: "Replace All Error",
                    message: expect.stringContaining("No match found"),
                }),
            );
        });
    });

    describe("onReset", () => {
        test("should reset all state", () => {
            const { result } = renderHook(() =>
                useFindReplace(
                    orderedParamColumns,
                    sortedIndices,
                    dataStore,
                    searchColumn,
                    config,
                    gridRef,
                    selectionSourceRef,
                    setFeedbackAlert,
                ),
            );

            act(() => {
                result.current.handleFind("cells");
            });

            expect(result.current.foundMatches).toHaveLength(3);
            expect(result.current.matchCount).toBe(3);

            act(() => {
                result.current.onReset();
            });

            expect(result.current.foundMatches).toHaveLength(0);
            expect(result.current.matchCount).toBeNull();
        });
    });
});