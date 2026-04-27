import useSlickGridReact from "@/react/hooks/useSlickGridReact";
import { renderHook, act } from "@testing-library/react";
import { describe, test, expect, beforeEach, vi } from "vitest";
import { observable, runInAction } from "mobx";
import type { TableChartReactConfig } from "@/react/components/TableChartReactWrapper";
import type { DataType, LoadedDataColumn } from "@/charts/charts";
import { createSlickGridMock } from "./testUtils/createSlickGridMock";

// Mock all the context hooks
vi.mock("@/react/context", () => ({
    useDataStore: () => mockDataStore,
    useChart: () => mockChart,
}));

vi.mock("@/react/hooks", () => ({
    useConfig: () => mockConfig,
    useChartID: () => "test-chart-id",
    useTheme: () => "dark",
    useOrderedParamColumns: () => mockOrderedParamColumns,
    useChartManager: () => mockChartManager,
}));

vi.mock("@/react/selectionHooks", () => ({
    useHighlightedIndices: () => mockHighlightedIndices,
}));

vi.mock("@/react/hooks/useSortedFilteredIndices", () => ({
    default: () => mockSortedIndices,
}));

let mockDataStore: any;
let mockConfig: TableChartReactConfig;
let mockChart: any;
let mockOrderedParamColumns: LoadedDataColumn<DataType>[];
let mockSortedIndices: Uint32Array;
let mockHighlightedIndices: number[];
let mockPermission: "edit" | "view";
let mockChartManager: any;

function setupGrid(result: { current: ReturnType<typeof useSlickGridReact> }) {
    const mockGridInstance = createSlickGridMock();
    const mockEvent = new CustomEvent("gridCreated", {
        detail: mockGridInstance,
    }) as any;

    act(() => {
        result.current.handleGridCreated(mockEvent);
    });

    const pubSub = mockGridInstance.slickGrid.getPubSubService() as any;
    const headerMenuCall = pubSub.subscribe.mock.calls.find(
        ([eventName]: [string]) => eventName === "onHeaderMenuCommand",
    );
    const headerMenuHandler = headerMenuCall?.[1];

    expect(headerMenuHandler).toBeTruthy();

    return {
        gridInstance: mockGridInstance,
        headerMenuHandler,
        columnsReorderedHandler: (mockGridInstance.slickGrid.onColumnsReordered as any).subscribe.mock.calls[0]?.[0],
        columnsResizedHandler: (mockGridInstance.slickGrid.onColumnsResized as any).subscribe.mock.calls[0]?.[0],
    };
}

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
            {
                field: "id",
                name: "Id",
                datatype: "unique",
                values: ["AB", "A"],
                data: new Uint32Array([65, 66, 65, 0]),
                stringLength: 2,
                editable: true,
            }
        ] as any;

        mockSortedIndices = new Uint32Array([0, 1, 2]);
        mockHighlightedIndices = [];
        mockPermission = "edit";
        mockChartManager = {
            config: { permission: mockPermission },
            analyzeColumnRemoval: vi.fn().mockResolvedValue({
                dataSourceName: "test-ds",
                columnName: "age",
                currentViewCharts: [],
                savedViews: [],
            }),
            saveState: vi.fn(),
            viewManager: {
                saveView: vi.fn().mockResolvedValue(undefined),
            },
            loadColumnSet: vi.fn((_columns: string[], _dsName: string, callback: () => void) => {
                callback();
            }),
        };

        // Setup mock config with observable properties
        mockConfig = observable({
            include_index: false,
            column_widths: {},
            order: {},
            sort: undefined as { columnId: string; ascending: boolean } | undefined,
            param: ["age", "name", "id"],
        }) as any;

        mockDataStore = {
            dataHighlighted: vi.fn(),
            dataChanged: vi.fn(),
            addColumn: vi.fn(),
            cleanColumnData: vi.fn(),
            removeColumn: vi.fn(),
            addListener: vi.fn(),
            columns: mockOrderedParamColumns,
            size: 3,
            name: "test-ds",
            columnIndex: Object.fromEntries(
                mockOrderedParamColumns.map((column) => [column.field, column]),
            ),
        };

        mockChart = {
            setAddColumnDialogOpener: vi.fn(),
            setGridRef: vi.fn(),
            activeQueries: {
                activeParams: vi.fn(() => mockConfig.param),
            },
            setParams: vi.fn((nextParam: Array<string | Record<string, unknown>>) => {
                runInAction(() => {
                    mockConfig.param = nextParam as typeof mockConfig.param;
                });
            }),
        };
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

            expect(result.current.columnDefs).toHaveLength(3);
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

            expect(result.current.columnDefs).toHaveLength(4);
            expect(result.current.columnDefs[0].field).toBe("__index__");
            expect(result.current.columnDefs[1].field).toBe("age");
            expect(result.current.columnDefs[2].field).toBe("name");
        });

        test("should include remove column in editable column header menu", () => {
            const { result } = renderHook(() => useSlickGridReact());

            const ageColumn = result.current.columnDefs.find((col) => col.field === "age");
            expect(ageColumn?.header?.menu?.commandItems).toEqual([
                {
                    command: "find-replace",
                    title: "Find & Replace",
                    iconCssClass: "mdi mdi-magnify",
                },
                {
                    command: "bulk-edit",
                    title: "Bulk Edit",
                    iconCssClass: "mdi mdi-table-edit",
                },
                {
                    command: "remove-column",
                    title: "Remove Column",
                    iconCssClass: "mdi mdi-delete",
                },
            ]);
        });

        test("should not include remove column for non-editable column", () => {
            const { result } = renderHook(() => useSlickGridReact());

            const nameColumn = result.current.columnDefs.find((col) => col.field === "name");
            expect(nameColumn?.header?.menu?.commandItems).toEqual([
                {
                    command: "find-replace",
                    title: "Find",
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

        test("should not have editor for non-editable column", () => {
            const { result } = renderHook(() => useSlickGridReact());

            const uniqueCol = result.current.columnDefs.find((col) => col.field === "name");
            expect(uniqueCol?.editor).toBeNull();
        });

        test("should have editor for editable integer column", () => {
            const { result } = renderHook(() => useSlickGridReact());

            const uniqueCol = result.current.columnDefs.find((col) => col.field === "age");
            expect(uniqueCol?.editor).toBeTruthy();
        });

        test("should have the custom editor for the unique editable column", () => {
            const { result } = renderHook(() => useSlickGridReact());

            const uniqueCol = result.current.columnDefs.find((col) => col.field === "id");
            expect(uniqueCol?.editor?.model?.name).toBe("TextEditorWithMaxLength");
            expect(uniqueCol?.editor?.maxLength).toBe(2);
        });
    });

    describe("handleGridCreated", () => {
        test("should set grid ref and initialize grid", () => {
            const { result } = renderHook(() => useSlickGridReact());

            const { gridInstance: mockGridInstance } = setupGrid(result);

            expect(mockGridInstance.slickGrid.setData).toHaveBeenCalled();
            expect(mockGridInstance.slickGrid.render).toHaveBeenCalled();
        });

        test("should open the removal impact dialog when charts use the column", async () => {
            mockChartManager.analyzeColumnRemoval.mockResolvedValueOnce({
                dataSourceName: "test-ds",
                columnName: "age",
                currentViewCharts: [
                    {
                        chartTitle: "Scatter",
                        chartType: "scatter_plot",
                        chartTypeLabel: "scatter_plot",
                        action: "remove_chart",
                        usage: "param",
                        usagePaths: ["param"],
                    },
                ],
                savedViews: [],
            });
            const { result } = renderHook(() => useSlickGridReact());
            const { headerMenuHandler } = setupGrid(result);

            await act(async () => {
                headerMenuHandler({
                    column: { field: "age" },
                    command: "remove-column",
                });
                await Promise.resolve();
            });

            expect(mockChartManager.analyzeColumnRemoval).toHaveBeenCalledWith(
                "test-ds",
                "age",
                "test-chart-id",
            );
            expect(mockDataStore.removeColumn).not.toHaveBeenCalled();
            expect(result.current.pendingColumnRemoval).toEqual(
                expect.objectContaining({
                    columnName: "age",
                }),
            );
        });

        test("should respect updated permissions in existing header menu handlers", async () => {
            const { result, rerender } = renderHook(() => useSlickGridReact());
            const { headerMenuHandler } = setupGrid(result);

            mockPermission = "view";
            mockChartManager.config.permission = "view";
            rerender();

            await act(async () => {
                headerMenuHandler({
                    column: { field: "age" },
                    command: "remove-column",
                });
                await Promise.resolve();
            });

            expect(mockDataStore.removeColumn).not.toHaveBeenCalled();
            expect(mockChartManager.analyzeColumnRemoval).not.toHaveBeenCalled();
            expect(result.current.pendingColumnRemoval).toBeNull();
        });
    });

    describe("column removal", () => {
        test("should open a confirmation dialog when there are no impacts", async () => {
            const { result } = renderHook(() => useSlickGridReact());
            const { headerMenuHandler } = setupGrid(result);

            await act(async () => {
                mockChartManager.analyzeColumnRemoval.mockResolvedValueOnce({
                    dataSourceName: "test-ds",
                    columnName: "age",
                    currentViewCharts: [],
                    savedViews: [],
                });
                headerMenuHandler({
                    column: { field: "age" },
                    command: "remove-column",
                });
                await Promise.resolve();
            });

            expect(mockDataStore.removeColumn).not.toHaveBeenCalled();
            expect(result.current.pendingColumnRemoval).toEqual(
                expect.objectContaining({
                    columnName: "age",
                }),
            );
        });

        test("should remove the column after confirmation when there are no impacts", async () => {
            const { result } = renderHook(() => useSlickGridReact());
            const { headerMenuHandler } = setupGrid(result);

            await act(async () => {
                mockChartManager.analyzeColumnRemoval.mockResolvedValueOnce({
                    dataSourceName: "test-ds",
                    columnName: "age",
                    currentViewCharts: [],
                    savedViews: [],
                });
                headerMenuHandler({
                    column: { field: "age" },
                    command: "remove-column",
                });
                await Promise.resolve();
            });

            await act(async () => {
                await result.current.confirmColumnRemoval();
            });

            expect(mockDataStore.removeColumn).toHaveBeenCalledWith("age", true, true);
            expect(mockChartManager.viewManager.saveView).toHaveBeenCalledTimes(1);
            expect(mockChartManager.saveState).not.toHaveBeenCalled();
            expect(result.current.pendingColumnRemoval).toBeNull();
        });

        test("should block removal when impacts exist in saved views", async () => {
            mockChartManager.analyzeColumnRemoval.mockResolvedValueOnce({
                dataSourceName: "test-ds",
                columnName: "age",
                currentViewCharts: [],
                savedViews: [
                    {
                        viewName: "Other View",
                        charts: [
                            {
                                chartTitle: "Scatter",
                                chartType: "scatter_plot",
                                chartTypeLabel: "scatter_plot",
                                action: "remove_chart",
                                usage: "param",
                                usagePaths: ["param"],
                            },
                        ],
                    },
                ],
            });

            const { result } = renderHook(() => useSlickGridReact());
            const { headerMenuHandler } = setupGrid(result);

            await act(async () => {
                headerMenuHandler({
                    column: { field: "age" },
                    command: "remove-column",
                });
                await Promise.resolve();
            });

            expect(result.current.pendingColumnRemoval?.impact.savedViews).toHaveLength(1);
            expect(mockDataStore.removeColumn).not.toHaveBeenCalled();
        });

        test("should surface an error and block removal when column analysis fails", async () => {
            mockChartManager.analyzeColumnRemoval.mockRejectedValueOnce(
                new Error("Failed to check column usage in saved views: Other View"),
            );

            const { result } = renderHook(() => useSlickGridReact());
            const { headerMenuHandler } = setupGrid(result);

            await act(async () => {
                headerMenuHandler({
                    column: { field: "age" },
                    command: "remove-column",
                });
                await Promise.resolve();
            });

            expect(result.current.pendingColumnRemoval).toBeNull();
            expect(mockDataStore.removeColumn).not.toHaveBeenCalled();
            expect(result.current.feedbackAlert).toEqual(
                expect.objectContaining({
                    type: "error",
                    title: "Remove Column Error",
                    message: "Failed to check column usage in saved views: Other View",
                }),
            );
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

            expect(result.current.columnDefs).toHaveLength(3);

            act(() => {
                runInAction(() => {
                    mockConfig.include_index = true;
                });
            });

            rerender();

            expect(result.current.columnDefs).toHaveLength(4);
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

    describe("add column", () => {
        test("should expose cloneable text columns", () => {
            mockDataStore.columns = [
                ...mockOrderedParamColumns,
                {
                    field: "hidden",
                    name: "Hidden",
                    datatype: "text",
                },
                {
                    field: "genes|TP53",
                    name: "TP53",
                    datatype: "double",
                    subgroup: "genes",
                },
            ];
            mockDataStore.columnIndex.hidden = {
                field: "hidden",
                name: "Hidden",
                datatype: "text",
            };
            const { result } = renderHook(() => useSlickGridReact());

            expect(result.current.cloneableColumns).toEqual([
                { field: "age", name: "Age", datatype: "integer" },
                { field: "hidden", name: "Hidden", datatype: "text" },
                { field: "id", name: "Id", datatype: "unique", stringLength: 2 },
                { field: "name", name: "Name", datatype: "text" },
            ]);
        });

        test("should create a new empty typed column and insert it at the requested position", () => {
            const { result } = renderHook(() => useSlickGridReact());

            act(() => {
                result.current.handleAddColumn({
                    name: "annotation",
                    datatype: "double",
                    position: 2,
                });
            });

            expect(mockDataStore.addColumn).toHaveBeenCalledTimes(1);
            expect(mockChart.setParams).toHaveBeenCalledWith(["age", "annotation", "name", "id"]);
            expect(mockConfig.param).toEqual(["age", "annotation", "name", "id"]);
            expect(mockConfig.order).toEqual({
                age: 0,
                annotation: 1,
                name: 2,
                id: 3,
            });
            expect(mockDataStore.dataChanged).toHaveBeenCalledWith(["annotation"]);
        });

        test("should insert a new column using the current unsaved grid order", () => {
            const { result } = renderHook(() => useSlickGridReact());
            const { gridInstance } = setupGrid(result);

            gridInstance.slickGrid.getColumns = vi.fn(() => [
                { field: "id", width: 100 },
                { field: "age", width: 100 },
                { field: "name", width: 100 },
            ] as any);

            act(() => {
                result.current.handleAddColumn({
                    name: "annotation",
                    datatype: "double",
                    position: 2,
                });
            });

            expect(mockChart.setParams).toHaveBeenCalledWith(["id", "annotation", "age", "name"]);
            expect(mockConfig.order).toEqual({
                id: 0,
                annotation: 1,
                age: 2,
                name: 3,
            });
        });

        test("should preserve existing column widths when adding a column", () => {
            mockConfig.column_widths = { age: 150, name: 220 };
            const { result } = renderHook(() => useSlickGridReact());

            act(() => {
                result.current.handleAddColumn({
                    name: "annotation",
                    datatype: "text",
                    position: 2,
                });
            });

            expect(mockConfig.column_widths).toEqual({ age: 150, name: 220 });
        });

        test("should keep an active-link param intact when inserting near its expanded columns", () => {
            const activeLinkParam = {
                columns: [],
                fields: ["gene_a", "gene_b"],
                initialize: async () => undefined,
            };
            mockConfig.param = ["gene_a", "gene_b", "age"] as typeof mockConfig.param;
            mockChart.activeQueries.activeParams.mockReturnValue([activeLinkParam, "age"]);
            mockOrderedParamColumns = [
                {
                    field: "gene_a",
                    name: "gene_a",
                    datatype: "double",
                    data: new Float32Array([1, 2, 3]),
                    editable: false,
                },
                {
                    field: "gene_b",
                    name: "gene_b",
                    datatype: "double",
                    data: new Float32Array([4, 5, 6]),
                    editable: false,
                },
                mockOrderedParamColumns[0],
            ] as any;
            mockDataStore.columns = mockOrderedParamColumns;
            mockDataStore.columnIndex = Object.fromEntries(
                mockOrderedParamColumns.map((column) => [column.field, column]),
            );

            const { result } = renderHook(() => useSlickGridReact());
            const { gridInstance } = setupGrid(result);
            vi.spyOn(gridInstance.slickGrid, "getColumns").mockReturnValue([
                { id: "gene_a", field: "gene_a", width: 100 },
                { id: "gene_b", field: "gene_b", width: 100 },
                { id: "age", field: "age", width: 100 },
            ]);

            act(() => {
                result.current.handleAddColumn({
                    name: "annotation",
                    datatype: "text",
                    position: 2,
                });
            });

            expect(mockChart.setParams).toHaveBeenCalledWith([
                expect.objectContaining({
                    fields: ["gene_a", "gene_b"],
                }),
                "annotation",
                "age",
            ]);
            expect(mockConfig.order).toEqual({
                gene_a: 0,
                gene_b: 1,
                annotation: 2,
                age: 3,
            });
        });

        test("should pass clone metadata through to the model", () => {
            const { result } = renderHook(() => useSlickGridReact());

            mockDataStore.columnIndex = {
                age: {
                    field: "age",
                    name: "Age",
                    datatype: "integer",
                    data: new Float32Array([1, 2, 3]),
                    getValue: (idx: number) => [1, 2, 3][idx],
                },
            };

            act(() => {
                result.current.handleAddColumn({
                    name: "copied_age",
                    datatype: "text",
                    cloneColumn: "age",
                    position: 3,
                });
            });

            expect(mockDataStore.addColumn).toHaveBeenCalledTimes(1);
            const [column] = mockDataStore.addColumn.mock.calls[0];
            expect(column).toMatchObject({
                name: "copied_age",
                field: "copied_age",
                datatype: "integer",
                editable: true,
            });
        });

        test("should load a hidden real column before cloning it", async () => {
            const { result } = renderHook(() => useSlickGridReact());

            mockDataStore.columns = [
                ...mockOrderedParamColumns,
                {
                    field: "n_genes_by_counts",
                    name: "n_genes_by_counts",
                    datatype: "double",
                },
            ];
            mockDataStore.columnIndex.n_genes_by_counts = {
                field: "n_genes_by_counts",
                name: "n_genes_by_counts",
                datatype: "double",
                data: null,
            };
            mockChartManager.loadColumnSet = vi.fn((_columns: string[], _dsName: string, callback: () => void) => {
                mockDataStore.columnIndex.n_genes_by_counts = {
                    field: "n_genes_by_counts",
                    name: "n_genes_by_counts",
                    datatype: "double",
                    data: new Float32Array([1, 2, 3]),
                };
                callback();
            });

            await act(async () => {
                await result.current.handleAddColumn({
                    name: "active link col",
                    datatype: "text",
                    cloneColumn: "n_genes_by_counts",
                    position: 2,
                });
            });

            expect(mockChartManager.loadColumnSet).toHaveBeenCalledWith(
                ["n_genes_by_counts"],
                "test-ds",
                expect.any(Function),
            );
            expect(mockDataStore.addColumn).toHaveBeenCalledTimes(1);
        });

        test("should reject duplicate column names", () => {
            mockDataStore.columnIndex = { age: { field: "age" } };
            const { result } = renderHook(() => useSlickGridReact());

            act(() => {
                result.current.handleAddColumn({
                    name: "age",
                    datatype: "text",
                    position: 1,
                });
            });

            expect(mockDataStore.addColumn).not.toHaveBeenCalled();
            expect(result.current.feedbackAlert).toEqual(
                expect.objectContaining({
                    type: "error",
                    title: "Add Column Error",
                    message: "Column age already exists",
                }),
            );
        });
    });

    describe("bulk edit", () => {
        test("should fill all visible cells in an editable numeric column", () => {
            const { result } = renderHook(() => useSlickGridReact());

            act(() => {
                result.current.handleBulkEdit({
                    action: "fill-all",
                    columnName: "age",
                    value: "42",
                });
            });

            expect(Array.from(mockOrderedParamColumns[0].data as ArrayLike<number>)).toEqual([42, 42, 42]);
            expect(mockDataStore.dataChanged).toHaveBeenCalledWith(["age"]);
        });

        test("should fill only empty cells in a text column", () => {
            mockOrderedParamColumns[1].editable = true;
            mockOrderedParamColumns[1].values = ["", "Bob", "Charlie"];
            mockOrderedParamColumns[1].data = new Uint8Array([0, 1, 0]);
            const { result } = renderHook(() => useSlickGridReact());

            act(() => {
                result.current.handleBulkEdit({
                    action: "fill-empty",
                    columnName: "name",
                    value: "Filled",
                });
            });

            expect(mockOrderedParamColumns[1].values).toContain("Filled");
            expect(Array.from(mockOrderedParamColumns[1].data as ArrayLike<number>)).toEqual([3, 1, 3]);
            expect(mockDataStore.dataChanged).toHaveBeenCalledWith(["name"]);
        });
    });

    describe("column width retention", () => {
        test("should sync unsaved resized widths into config immediately", () => {
            const { result } = renderHook(() => useSlickGridReact());
            const { columnsResizedHandler, gridInstance } = setupGrid(result);

            gridInstance.slickGrid.getColumns = vi.fn(() => [
                { field: "age", width: 180 },
                { field: "name", width: 260 },
                { field: "id", width: 100 },
            ] as any);

            act(() => {
                columnsResizedHandler?.({}, {});
            });

            expect(mockConfig.column_widths).toEqual({
                age: 180,
                name: 260,
            });
        });
    });
});
