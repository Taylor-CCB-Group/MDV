import ViewManager from "@/charts/ViewManager";
import { beforeEach, describe, expect, test, vi } from "vitest";

vi.mock("@/utilities/Screenshot", () => ({
    default: vi.fn(() => Promise.resolve("data:image/png;base64,test")),
}));

vi.mock("@/dataloaders/DataLoaderUtil", () => ({
    getPostData: vi.fn(),
    getProjectRoot: vi.fn(() => ""),
}));

describe("ViewManager.saveView", () => {
    let chartManager: any;

    beforeEach(() => {
        chartManager = {
            getState: vi.fn(() => ({
                view: {
                    dataSources: {},
                    initialCharts: {},
                },
                currentView: "view-a",
                all_views: ["view-a", "view-b"],
                updatedColumns: {},
                metadata: {},
                chartErrors: [],
            })),
            _callListeners: vi.fn(),
            contentDiv: document.createElement("div"),
        };
        Reflect.defineProperty(window, "mdv", {
            configurable: true,
            value: { chartManager },
        });
    });

    test("includes updatedViews in the state_saved payload when provided", async () => {
        const viewManager = new ViewManager("view-a", ["view-a", "view-b"]);
        const updatedViews = {
            "view-b": {
                dataSources: {},
                initialCharts: {},
            },
        };

        await viewManager.saveView(undefined, updatedViews as any);

        expect(chartManager._callListeners).toHaveBeenCalledWith(
            "state_saved",
            expect.objectContaining({
                updatedViews,
            }),
        );
    });
});

describe("ViewManager.changeView", () => {
    test("discards unsaved added columns before loading the next view", async () => {
        const discardPendingAddedColumns = vi.fn();
        const chartManager = {
            viewData: {
                dataSources: {
                    ds1: {},
                },
            },
            dsIndex: {
                ds1: { name: "ds1" },
            },
            dataSources: [
                {
                    name: "ds1",
                    dataStore: {
                        discardPendingAddedColumns,
                    },
                },
            ],
            gridStack: {
                destroy: vi.fn(),
            },
            removeAllCharts: vi.fn(),
            contentDiv: document.createElement("div"),
            viewLoader: vi.fn(() =>
                Promise.resolve({
                    dataSources: { ds1: {} },
                    initialCharts: { ds1: [] },
                }),
            ),
            filterRemovedColumnsFromViewData: vi.fn((data) => data),
            _init: vi.fn(() => Promise.resolve()),
            getState: vi.fn(() => ({
                view: {
                    dataSources: {},
                    initialCharts: {},
                },
                currentView: "view-b",
                all_views: ["view-a", "view-b"],
                updatedColumns: {},
                metadata: {},
                chartErrors: [],
            })),
        };
        Reflect.defineProperty(window, "mdv", {
            configurable: true,
            value: { chartManager },
        });

        const viewManager = new ViewManager("view-a", ["view-a", "view-b"]);

        viewManager.changeView("view-b");
        await Promise.resolve();
        await Promise.resolve();

        expect(discardPendingAddedColumns).toHaveBeenCalledTimes(1);
        expect(chartManager._init).toHaveBeenCalledTimes(1);
    });
});
