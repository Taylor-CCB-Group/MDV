import TableChartReactComponent from "@/react/components/TableChartReactComponent";
import { render, screen } from "@testing-library/react";
import { describe, test, expect, beforeEach, vi } from "vitest";
import { observable } from "mobx";
import { createSlickGridMock } from "./hooks/testUtils/createSlickGridMock";

// Mock SlickgridReact to avoid browser API issues
vi.mock("slickgrid-react", () => ({
    SlickgridReact: ({ gridId }: { gridId: string }) => (
        <div data-testid="slickgrid-component" data-grid-id={gridId}>
            SlickGrid
        </div>
    ),
}));

vi.mock("@/react/components/FindAndReplaceDialog", () => ({
    default: () => <div data-testid="mock-find-replace-dialog" />,
}));

vi.mock("@/react/components/AddTableColumnDialog", () => ({
    default: () => <div data-testid="mock-add-column-dialog" />,
}));

vi.mock("@/react/components/BulkEditColumnDialog", () => ({
    default: () => <div data-testid="mock-bulk-edit-dialog" />,
}));

vi.mock("@/charts/dialogs/ReusableAlertDialog", () => ({
    default: () => <div data-testid="mock-alert-dialog" />,
}));

// Mock the hooks
let mockSlickGridReactReturn: any;
let mockFindReplaceReturn: any;
let mockEditCellReturn: any;

vi.mock("@/react/hooks/useSlickGridReact", () => ({
    default: () => mockSlickGridReactReturn,
}));

vi.mock("@/react/hooks/useFindReplace", () => ({
    default: () => mockFindReplaceReturn,
}));

vi.mock("@/react/hooks/useEditCell", () => ({
    default: () => mockEditCellReturn,
}));

describe("TableChartReactComponent", () => {
    beforeEach(() => {
        mockSlickGridReactReturn = {
            chartId: "test-chart",
            config: observable({}),
            dataStore: {},
            gridRef: { current: createSlickGridMock() },
            selectionSourceRef: { current: null },
            isFindReplaceOpen: false,
            orderedParamColumns: [],
            searchColumn: null,
            orderedParamColumnsRef: { current: [] },
            sortedFilteredIndices: new Uint32Array([]),
            sortedFilteredIndicesRef: { current: new Uint32Array([]) },
            options: {},
            columnDefs: [],
            handleGridCreated: vi.fn(),
            isColumnEditable: false,
            onDialogClose: vi.fn(),
            feedbackAlert: null,
            setFeedbackAlert: vi.fn(),
            isAddColumnDialogOpen: false,
            cloneableColumns: [],
            addColumnDefaultPosition: 1,
            closeAddColumnDialog: vi.fn(),
            handleAddColumn: vi.fn(),
            isBulkEditDialogOpen: false,
            bulkEditColumn: null,
            closeBulkEditDialog: vi.fn(),
            handleBulkEdit: vi.fn(),
        };

        mockFindReplaceReturn = {
            matchCount: null,
            disableFindPrev: true,
            disableFindNext: true,
            handleFind: vi.fn(),
            handleFindNext: vi.fn(),
            handleFindPrev: vi.fn(),
            handleReplace: vi.fn(),
            handleReplaceAll: vi.fn(),
            onReset: vi.fn(),
        };

        mockEditCellReturn = {
            handleBeforeEditCell: vi.fn(),
            handleCellChange: vi.fn(),
        };
    });

    test("should render main container", () => {
        render(<TableChartReactComponent />);

        const container = screen.getByTestId("slickgrid-react-container");
        expect(container).toBeDefined();
    });

    test("should render SlickGrid component", () => {
        render(<TableChartReactComponent />);

        const slickGrid = screen.getByTestId("slickgrid-component");
        expect(slickGrid).toBeDefined();
    });

    test("should render FindAndReplaceDialog wrapper", () => {
        render(<TableChartReactComponent />);

        const dialogWrapper = screen.getByTestId("find-replace-dialog-wrapper");
        expect(dialogWrapper).toBeDefined();
    });

    test("should not show feedback alert by default", () => {
        render(<TableChartReactComponent />);

        const feedbackDialog = screen.queryByTestId("feedback-alert-dialog");
        expect(feedbackDialog).toBeNull();
    });
});