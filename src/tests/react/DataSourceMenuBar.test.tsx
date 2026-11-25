import { act, cleanup, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type { DataSource } from "@/charts/charts";
import DataStore from "@/datastore/DataStore";

const { chartManager, renameDialogMock } = vi.hoisted(() => ({
    chartManager: {
        config: { permission: "edit" },
        layoutMenus: { people: { show: vi.fn() } },
        createInfoAlert: vi.fn(),
    },
    renameDialogMock: vi.fn(),
}));

vi.mock("@/react/hooks", () => ({ useChartManager: () => chartManager }));
vi.mock("@/charts/dialogs/AddChartDialogReact", () => ({ default: vi.fn() }));
vi.mock("@/charts/dialogs/AnnotationDialogReact", () => ({ AnnotationDialogReact: vi.fn() }));
vi.mock("@/charts/dialogs/ColorPaletteWrapper", () => ({ default: vi.fn() }));
vi.mock("@/charts/dialogs/RenameDatasourceDialog", () => ({ default: renameDialogMock }));

import DataSourceMenuBar from "@/react/components/DataSourceMenuBar";

function createDataSource(displayName?: string): DataSource {
    const dataStore: DataStore = Object.assign(Object.create(DataStore.prototype), {
        name: "people",
        config: { name: "people", display_name: displayName },
        size: 6,
        filterSize: 6,
        listeners: {},
        removeAllFilters: vi.fn(),
    });
    return {
        name: "people",
        dataStore,
        charts: [],
        columns: [],
        size: 6,
        contentDiv: document.createElement("div"),
        menuBar: document.createElement("div"),
    };
}

function setFilterSize(dataSource: DataSource, filterSize: number) {
    act(() => {
        dataSource.dataStore.filterSize = filterSize;
        dataSource.dataStore._callListeners("filtered", null);
    });
}

beforeEach(() => {
    chartManager.config.permission = "edit";
});

afterEach(() => {
    cleanup();
    vi.clearAllMocks();
});

describe("DataSourceMenuBar", () => {
    test("shows the datasource name and its row count", () => {
        const { container } = render(<DataSourceMenuBar dataSource={createDataSource()} />);

        expect(screen.getByText("people")).toBeTruthy();
        expect(screen.getByLabelText("6 of 6 rows").textContent).toBe("6");
        expect(container.querySelector<HTMLElement>(".mdv-ds-count-fill")?.style.width).toBe("100%");
    });

    test("formats large row counts", () => {
        const dataSource = createDataSource();
        dataSource.dataStore.size = 12000;
        dataSource.dataStore.filterSize = 12000;
        render(<DataSourceMenuBar dataSource={dataSource} />);

        expect(screen.getByLabelText(`${(12000).toLocaleString()} of ${(12000).toLocaleString()} rows`)).toBeTruthy();
    });

    test("follows the filtered row count", () => {
        const dataSource = createDataSource();
        const { container } = render(<DataSourceMenuBar dataSource={dataSource} />);

        setFilterSize(dataSource, 3);

        expect(screen.getByLabelText("3 of 6 rows").textContent).toBe("3");
        expect(container.querySelector<HTMLElement>(".mdv-ds-count-fill")?.style.width).toBe("50%");
    });

    test("stops listening to the datastore when unmounted", () => {
        const dataSource = createDataSource();
        const { unmount } = render(<DataSourceMenuBar dataSource={dataSource} />);
        expect(Object.keys(dataSource.dataStore.listeners)).toHaveLength(1);

        unmount();

        expect(Object.keys(dataSource.dataStore.listeners)).toHaveLength(0);
    });

    test("resets the datastore's filters", () => {
        const dataSource = createDataSource();
        render(<DataSourceMenuBar dataSource={dataSource} />);

        fireEvent.click(screen.getByLabelText("Reset All Filters"));

        expect(dataSource.dataStore.removeAllFilters).toHaveBeenCalledOnce();
    });

    test("opens the layout menu of its own datasource", () => {
        render(<DataSourceMenuBar dataSource={createDataSource()} />);

        fireEvent.click(screen.getByLabelText("Change layout"));

        expect(chartManager.layoutMenus.people.show).toHaveBeenCalledOnce();
    });

    test("shows the full name as a tooltip, as it is truncated when space is short", () => {
        render(<DataSourceMenuBar dataSource={createDataSource()} />);

        expect(screen.getByLabelText("people").textContent).toBe("people");
    });

    test("shows the datasource's display name in place of its identifier", () => {
        render(<DataSourceMenuBar dataSource={createDataSource("Patients")} />);

        expect(screen.getByLabelText("Patients").textContent).toBe("Patients");
        expect(screen.queryByText("people")).toBeNull();
    });

    test("opens the rename dialog for its datastore from the pencil button, not the name", () => {
        const dataSource = createDataSource();
        render(<DataSourceMenuBar dataSource={dataSource} />);

        expect(screen.getByText("people").closest("button")).toBeNull();
        fireEvent.click(screen.getByRole("button", { name: "Rename datasource" }));

        expect(renameDialogMock).toHaveBeenCalledWith(dataSource.dataStore);
    });

    test("has no rename button when the project is view only", () => {
        chartManager.config.permission = "view";
        render(<DataSourceMenuBar dataSource={createDataSource("Patients")} />);

        expect(screen.queryByRole("button", { name: "Rename datasource" })).toBeNull();
        expect(screen.getByText("Patients")).toBeTruthy();
    });

    test("shows whether its pane is full screen", () => {
        const dataSource = createDataSource();
        const pane = document.createElement("div");
        pane.append(dataSource.menuBar);
        render(<DataSourceMenuBar dataSource={dataSource} />);
        // jsdom does not implement the Fullscreen API
        const showFullscreen = (element: Element | null) => {
            Object.defineProperty(document, "fullscreenElement", { value: element, configurable: true });
            act(() => {
                document.dispatchEvent(new Event("fullscreenchange"));
            });
        };

        showFullscreen(pane);
        expect(screen.getByLabelText("Exit Full Screen")).toBeTruthy();

        // another element going full screen, such as a single chart, is not this pane
        showFullscreen(document.createElement("div"));
        expect(screen.getByLabelText("Full Screen")).toBeTruthy();

        showFullscreen(null);
    });
});
