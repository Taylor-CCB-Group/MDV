import type { DataType, LoadedDataColumn } from "@/charts/charts";
import type DataStore from "@/datastore/DataStore";
import useEditCell from "@/react/hooks/useEditCell";
import { renderHook, act } from "@testing-library/react";
import { describe, test, expect, beforeEach, vi } from "vitest";
import type { OnBeforeEditCellEventArgs, OnCellChangeEventArgs } from "slickgrid-react";

describe("useEditCell", () => {
    let orderedParamColumns: LoadedDataColumn<DataType>[];
    let sortedIndices: Uint32Array;
    let dataStore: DataStore;
    let gridRef: React.MutableRefObject<any>;
    let setFeedbackAlert: ReturnType<typeof vi.fn>;
    let mockGrid: any;
    let orderedParamColumnsRef: React.MutableRefObject<LoadedDataColumn<DataType>[]>;
    let sortedIndicesRef: React.MutableRefObject<Uint32Array>;

    beforeEach(() => {
        // Reset all mocks before each test
        vi.clearAllMocks();

        // Setup mock grid
        mockGrid = {
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
            },
            {
                field: "cell_type",
                datatype: "text",
                data: new Uint8Array([0, 1, 2]),
                values: ["T cells", "B cells", "NK cells"],
                editable: true,
            },
            {
                field: "readonly",
                datatype: "text",
                data: [],
                values: [],
                editable: false,
            }
        ] as any;

        sortedIndices = new Uint32Array([0, 1, 2]);

        // Mock DataStore
        dataStore = {
            dataChanged: vi.fn(),
        } as any;

        orderedParamColumnsRef = { current: orderedParamColumns };
        sortedIndicesRef = { current: sortedIndices };
        gridRef = { current: { slickGrid: mockGrid } } as any;
        setFeedbackAlert = vi.fn();
    });

    test("should handle null grid reference gracefully", () => {
        const gridRefWithNull = { current: null } as any;
        
        const { result } = renderHook(() =>
            useEditCell(
                orderedParamColumnsRef,
                sortedIndicesRef,
                dataStore,
                gridRefWithNull,
                setFeedbackAlert,
            ),
        );
    
        // Should not crash when grid is null
        expect(result.current.handleBeforeEditCell).toBeDefined();
        expect(result.current.handleCellChange).toBeDefined();
    });

    describe("handleBeforeEditCell", () => {
        test("should store old cell value", () => {
            const { result } = renderHook(() =>
                useEditCell(
                    orderedParamColumnsRef,
                    sortedIndicesRef,
                    dataStore,
                    gridRef,
                    setFeedbackAlert,
                ),
            );

            const mockEvent = new CustomEvent("beforeEditCell", {
                detail: {
                    eventData: {},
                    args: {
                        item: { sample_count: 10 },
                        column: { field: "sample_count" },
                    } as OnBeforeEditCellEventArgs,
                },
            }) as any;

            act(() => {
                result.current.handleBeforeEditCell(mockEvent);
            });
        });
    });

    describe("handleCellChange", () => {
        test("should handle successful cell edit", () => {
            const { result } = renderHook(() =>
                useEditCell(
                    orderedParamColumnsRef,
                    sortedIndicesRef,
                    dataStore,
                    gridRef,
                    setFeedbackAlert,
                ),
            );

            // Storing old cell value
            const beforeEditEvent = new CustomEvent("beforeEditCell", {
                detail: {
                    eventData: {},
                    args: {
                        item: { sample_count: 10 },
                        column: { field: "sample_count" },
                    } as OnBeforeEditCellEventArgs,
                },
            }) as any;

            act(() => {
                result.current.handleBeforeEditCell(beforeEditEvent);
            });

            // Cell change event
            const changeEvent = new CustomEvent("cellChange", {
                detail: {
                    eventData: {},
                    args: {
                        row: 0,
                        column: { field: "sample_count" },
                        item: { sample_count: 15 },
                    } as OnCellChangeEventArgs,
                },
            }) as any;

            act(() => {
                result.current.handleCellChange(changeEvent);
            });

            expect(setFeedbackAlert).toHaveBeenCalledWith(
                expect.objectContaining({
                    type: "success",
                    title: "Edit Successful",
                }),
            );
            expect(dataStore.dataChanged).toHaveBeenCalledWith(["sample_count"]);
            expect(mockGrid.invalidate).toHaveBeenCalled();
            expect(mockGrid.render).toHaveBeenCalled();
        });

        test("should show error when column not found", () => {
            const { result } = renderHook(() =>
                useEditCell(
                    orderedParamColumnsRef,
                    sortedIndicesRef,
                    dataStore,
                    gridRef,
                    setFeedbackAlert,
                ),
            );

            const changeEvent = new CustomEvent("cellChange", {
                detail: {
                    eventData: {},
                    args: {
                        row: 0,
                        column: { field: "nonexistent" },
                        item: { nonexistent: "value" },
                    } as OnCellChangeEventArgs,
                },
            }) as any;

            act(() => {
                result.current.handleCellChange(changeEvent);
            });

            expect(setFeedbackAlert).toHaveBeenCalledWith(
                expect.objectContaining({
                    type: "error",
                    title: "Edit Error",
                    message: "Column not found",
                }),
            );
        });

        test("should show error when column is not editable", () => {
            const { result } = renderHook(() =>
                useEditCell(
                    orderedParamColumnsRef,
                    sortedIndicesRef,
                    dataStore,
                    gridRef,
                    setFeedbackAlert,
                ),
            );

            const changeEvent = new CustomEvent("cellChange", {
                detail: {
                    eventData: {},
                    args: {
                        row: 0,
                        column: { field: "readonly" },
                        item: { readonly: "" },
                    } as OnCellChangeEventArgs,
                },
            }) as any;

            act(() => {
                result.current.handleCellChange(changeEvent);
            });

            expect(setFeedbackAlert).toHaveBeenCalledWith(
                expect.objectContaining({
                    type: "error",
                    title: "Edit Error",
                    message: "Column readonly not editable",
                }),
            );
        });
    });
});

