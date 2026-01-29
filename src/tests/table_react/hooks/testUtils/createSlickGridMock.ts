import type {
    SlickgridReactInstance,
    BasePubSub,
    SlickGrid,
    SlickDataView,
    ExtensionService,
    FilterService,
    GridService,
    GridEventService,
    GridStateService,
    HeaderGroupingService,
    ResizerService,
    SortService,
    TreeDataService,
} from "slickgrid-react";
import { vi } from "vitest";

/**
 * Creates a mock of a SlickEvent.
 * SlickEvents have subscribe/unsubscribe methods for event handling.
 */
function createMockSlickEvent<T = any>(): SlickGrid["onSelectedRowsChanged"] {
    return {
        subscribe: vi.fn(),
        unsubscribe: vi.fn(),
    } as unknown as SlickGrid["onSelectedRowsChanged"];
}

/**
 * Creates a type-safe mock of SlickgridReactInstance for testing.
 * This mock satisfies TypeScript's type checking without implementing behavior,
 * keeping tests focused on pure logic rather than integration details.
 */
export function createSlickGridMock(): SlickgridReactInstance {
    const mockPubSub: BasePubSub = {
        subscribe: vi.fn(() => ({
            unsubscribe: vi.fn(),
        })),
        publish: vi.fn(),
    };

    const mockSlickGrid = {
        setData: vi.fn(),
        render: vi.fn(),
        invalidate: vi.fn(),
        gotoCell: vi.fn(),
        setSortColumn: vi.fn(),
        setSortColumns: vi.fn(),
        getSortColumns: vi.fn(() => []),
        getSelectedRows: vi.fn(() => []),
        setSelectedRows: vi.fn(),
        scrollRowIntoView: vi.fn(),
        getColumns: vi.fn(() => []),
        getPubSubService: vi.fn(() => mockPubSub),
        onSelectedRowsChanged: createMockSlickEvent(),
        onSort: createMockSlickEvent(),
        onColumnsResized: createMockSlickEvent(),
        onColumnsReordered: createMockSlickEvent(),
    } as unknown as SlickGrid;

    return {
        slickGrid: mockSlickGrid,
        element: document.createElement("div"),
        dataView: {} as SlickDataView,
        dispose: vi.fn(),
        extensionService: {} as ExtensionService,
        filterService: {} as FilterService,
        gridService: {} as GridService,
        gridEventService: {} as GridEventService,
        gridStateService: {} as GridStateService,
        headerGroupingService: {} as HeaderGroupingService,
        resizerService: {} as ResizerService,
        sortService: {} as SortService,
        treeDataService: {} as TreeDataService,
    } satisfies SlickgridReactInstance;
}

