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
            isSelectingRef: { current: false },
            isFindReplaceOpen: false,
            orderedParamColumns: [],
            searchColumn: null,
            orderedParamColumnsRef: { current: [] },
            sortedIndices: new Uint32Array([]),
            sortedIndicesRef: { current: new Uint32Array([]) },
            options: {},
            columnDefs: [],
            handleGridCreated: vi.fn(),
            isColumnEditable: false,
            onDialogClose: vi.fn(),
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