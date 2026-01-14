import useSlickGridReact from "@/react/hooks/useSlickGridReact";
import { renderHook, act, waitFor } from "@testing-library/react";
import { describe, test, expect, beforeEach, vi } from "vitest";
import { observable, runInAction } from "mobx";
import type { TableChartReactConfig } from "@/react/components/TableChartReactWrapper";
import type { DataType, LoadedDataColumn } from "@/charts/charts";

// Mock all the context hooks
vi.mock("@/react/context", () => ({
    useDataStore: () => mockDataStore,
    useChart: () => mockChart,
}));

vi.mock("@/react/hooks", () => ({
    useConfig: () => mockConfig,
    useChartID: () => "test-chart-id",
    useOrderedParamColumns: () => mockOrderedParamColumns,
}));

vi.mock("@/react/selectionHooks", () => ({
    useHighlightedIndex: () => mockHighlightedIndex,
    useHighlightedIndices: () => mockHighlightedIndex,
}));

vi.mock("@/react/hooks/useSortedIndices", () => ({
    default: () => mockSortedIndices,
}));

let mockDataStore: any;
let mockConfig: TableChartReactConfig;
let mockChart: any;
let mockOrderedParamColumns: LoadedDataColumn<DataType>[];
let mockSortedIndices: Uint32Array;
let mockHighlightedIndex: number;
let mockHighlightedIndices: number[];

describe("useSlickGridReact", () => {
    beforeEach(() => {
        vi.clearAllMocks();

        // Setup mock columns
        mockOrderedParamColumns = [
            {
                field: "age",
                name: "Age",
                datatype: "integer",
                data: new Uint32Array([25, 30, 20]),
                editable: true,
            },
            {
                field: "name",
                name: "Name",
                datatype: "text",
                data: new Uint8Array([0, 1, 2]),
                values: ["Alice", "Bob", "Charlie"],
                editable: false,
            },
        ] as any;

        mockSortedIndices = new Uint32Array([0, 1, 2]);
        mockHighlightedIndex = -1;
        mockHighlightedIndices = [];

        // Setup mock config with observable properties
        mockConfig = observable({
            include_index: false,
            column_widths: {},
            order: {},
            sort: undefined as { columnId: string; ascending: boolean } | undefined,
        }) as any;

        mockDataStore = {
            dataHighlighted: vi.fn(),
        };

        mockChart = {};
    });

    describe("initialization", () => {
        test("should initialize with default state", () => {
            const { result } = renderHook(() => useSlickGridReact());

            expect(result.current.isFindReplaceOpen).toBe(false);
            expect(result.current.searchColumn).toBeNull();
            expect(result.current.chartId).toBe("test-chart-id");
        });

        test("should create column definitions without index column", () => {
            const { result } = renderHook(() => useSlickGridReact());

            expect(result.current.columnDefs).toHaveLength(2);
            expect(result.current.columnDefs[0].field).toBe("age");
            expect(result.current.columnDefs[1].field).toBe("name");
        });

        test("should create column definitions with index column", () => {
            act(() => {
                runInAction(() => {
                    mockConfig.include_index = true;
                });
            });

            const { result } = renderHook(() => useSlickGridReact());

            expect(result.current.columnDefs).toHaveLength(3);
            expect(result.current.columnDefs[0].field).toBe("__index__");
            expect(result.current.columnDefs[1].field).toBe("age");
            expect(result.current.columnDefs[2].field).toBe("name");
        });

        test("should include find-replace menu item in header", () => {
            const { result } = renderHook(() => useSlickGridReact());

            const ageColumn = result.current.columnDefs.find((col) => col.field === "age");
            expect(ageColumn?.header?.menu?.commandItems).toEqual([
                {
                    command: "find-replace",
                    title: "Find & Replace",
                    iconCssClass: "mdi mdi-magnify",
                },
            ]);
        });

        test("should apply column widths from config", () => {
            act(() => {
                runInAction(() => {
                    mockConfig.column_widths = { age: 150, name: 200 };
                });
            });

            const { result } = renderHook(() => useSlickGridReact());

            const ageColumn = result.current.columnDefs.find((col) => col.field === "age");
            const nameColumn = result.current.columnDefs.find((col) => col.field === "name");

            expect(ageColumn?.width).toBe(150);
            expect(nameColumn?.width).toBe(200);
        });
    });

    describe("isColumnEditable", () => {
        test("should return false when no column selected", () => {
            const { result } = renderHook(() => useSlickGridReact());

            expect(result.current.isColumnEditable).toBe(false);
        });

        test("should return true for editable column", () => {
            const { result } = renderHook(() => useSlickGridReact());

            act(() => {
                result.current.setSearchColumn("age");
            });

            expect(result.current.isColumnEditable).toBe(true);
        });
    });

    describe("handleGridCreated", () => {
        test("should set grid ref and initialize grid", () => {
            const { result } = renderHook(() => useSlickGridReact());

            const mockGrid = {
                setData: vi.fn(),
                render: vi.fn(),
                onSelectedRowsChanged: { subscribe: vi.fn(), unsubscribe: vi.fn() },
                onSort: { subscribe: vi.fn(), unsubscribe: vi.fn()  },
                getPubSubService: vi.fn(() => ({
                    subscribe: vi.fn(),
                    unsubscribe: vi.fn() 
                })),
                onColumnsResized: { subscribe: vi.fn(), unsubscribe: vi.fn() },
                onColumnsReordered: { subscribe: vi.fn(), unsubscribe: vi.fn() },
            };

            const mockGridInstance = {
                slickGrid: mockGrid,
            } as any;

            const mockEvent = new CustomEvent("gridCreated", {
                detail: mockGridInstance,
            }) as any;

            act(() => {
                result.current.handleGridCreated(mockEvent);
            });

            expect(mockGrid.setData).toHaveBeenCalled();
            expect(mockGrid.render).toHaveBeenCalled();
        });
    });

    describe("find and replace dialog", () => {
        test("should open find replace dialog", () => {
            const { result } = renderHook(() => useSlickGridReact());

            act(() => {
                result.current.setSearchColumn("age");
                result.current.setIsFindReplaceOpen(true);
            });

            expect(result.current.isFindReplaceOpen).toBe(true);
            expect(result.current.searchColumn).toBe("age");
        });

        test("should close dialog and reset search column", () => {
            const { result } = renderHook(() => useSlickGridReact());

            act(() => {
                result.current.setSearchColumn("age");
                result.current.setIsFindReplaceOpen(true);
            });

            expect(result.current.isFindReplaceOpen).toBe(true);

            act(() => {
                result.current.onDialogClose();
            });

            expect(result.current.isFindReplaceOpen).toBe(false);
            expect(result.current.searchColumn).toBeNull();
        });
    });

    describe("column config updates", () => {
        test("should update column defs when include_index changes", () => {
            const { result, rerender } = renderHook(() => useSlickGridReact());

            expect(result.current.columnDefs).toHaveLength(2);

            act(() => {
                runInAction(() => {
                    mockConfig.include_index = true;
                });
            });

            rerender();

            expect(result.current.columnDefs).toHaveLength(3);
            expect(result.current.columnDefs[0].field).toBe("__index__");
        });

        test("should update column widths when config changes", () => {
            const { result, rerender } = renderHook(() => useSlickGridReact());

            const ageColumn = result.current.columnDefs.find((col) => col.field === "age");
            expect(ageColumn?.width).toBe(100); // default

            act(() => {
                runInAction(() => {
                    mockConfig.column_widths = { age: 250 };
                });
            });

            rerender();

            const updatedAgeColumn = result.current.columnDefs.find((col) => col.field === "age");
            expect(updatedAgeColumn?.width).toBe(250);
        });
    });
});

