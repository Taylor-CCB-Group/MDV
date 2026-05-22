import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import DataStore from "@/datastore/DataStore";

const { baseCloseMock, createMdvPortalMock } = vi.hoisted(() => ({
    baseCloseMock: vi.fn(),
    createMdvPortalMock: vi.fn(),
}));

vi.mock("@/react/react_utils", () => ({
    createMdvPortal: createMdvPortalMock,
}));

vi.mock("@/charts/dialogs/ColorPaletteComponent", () => ({
    default: function MockColorPaletteComponent() {
        return null;
    },
}));

vi.mock("@/utilities/Dialog.js", () => {
    class MockBaseDialog {
        config: Record<string, unknown>;
        dialog: HTMLDivElement;
        outer: HTMLDivElement;
        constructor(config: Record<string, unknown>, content: unknown) {
            this.config = config;
            this.outer = document.createElement("div");
            this.dialog = document.createElement("div");
            this.outer.append(this.dialog);
            document.body.append(this.outer);
            this.init(content);
        }
        init(_content: unknown) {}
        close() {
            baseCloseMock();
            this.outer.remove();
        }
    }

    return {
        BaseDialog: MockBaseDialog,
    };
});

import ColorPalette from "@/charts/dialogs/ColorPaletteWrapper";

afterEach(() => {
    createMdvPortalMock.mockReset();
    baseCloseMock.mockReset();
    document.body.innerHTML = "";
});

beforeEach(() => {
    createMdvPortalMock.mockReturnValue({
        unmount: vi.fn(),
    });
});

function createDataSource(
    name: string,
    overrides: Partial<{
        linkColumns: Array<{ columns: string[]; dataSource: string }>;
        syncColumnColors: Array<{
            columns: Array<{ col: string; link_to: string }>;
            dataSource: string;
        }>;
    }> = {},
) {
    const dataStore: DataStore = Object.assign(Object.create(DataStore.prototype), {
        dataChanged: vi.fn(),
        getColumnColors: vi.fn(() => ["#111111", "#222222"]),
        getColumnList: vi.fn(() => [{ field: "cell_type", name: "Cell type" }]),
        getColumnValues: vi.fn(() => ["T cell", "B cell"]),
        linkColumns: overrides.linkColumns ?? [],
        name,
        setColumnColors: vi.fn(),
        syncColumnColors: overrides.syncColumnColors ?? [],
    });
    return dataStore;
}

function createChartManagerDataSource(dataStore: ReturnType<typeof createDataSource>) {
    return {
        charts: [],
        columns: [],
        contentDiv: document.createElement("div"),
        dataStore,
        menuBar: document.createElement("div"),
        name: dataStore.name,
        size: 0,
    };
}

function createChartManager(
    dataSources: ReturnType<typeof createChartManagerDataSource>[],
    mainDataSource: ReturnType<typeof createDataSource>,
) {
    return {
        _sync_colors: vi.fn(),
        dataSources,
        getDataSource: vi.fn(() => mainDataSource),
    };
}

function setupDialog({
    mainDataSource = createDataSource("main-ds"),
    extraSources = [],
}: {
    mainDataSource?: ReturnType<typeof createDataSource>;
    extraSources?: ReturnType<typeof createChartManagerDataSource>[];
} = {}) {
    const mainSource = createChartManagerDataSource(mainDataSource);
    const chartManager = createChartManager([mainSource, ...extraSources], mainDataSource);
    const dialog = new ColorPalette(chartManager, mainSource);
    return { chartManager, dialog, mainDataSource, mainSource };
}

describe("ColorPalette wrapper", () => {
    test("mounts the React component with the selected datasource", () => {
        const { chartManager, dialog } = setupDialog();

        expect(chartManager.getDataSource).toHaveBeenCalledWith("main-ds");
        expect(createMdvPortalMock).toHaveBeenCalledTimes(1);
        expect(dialog.dialog.style.display).toBe("flex");

        const [, container] = createMdvPortalMock.mock.calls[0];
        expect(container).toBeInstanceOf(HTMLDivElement);
    });

    test("applyColors updates the selected datasource and linked color sync targets", () => {
        const mainDataSource = createDataSource("main-ds");
        const syncedTarget = createDataSource("synced-target", {
            syncColumnColors: [
                {
                    dataSource: "main-ds",
                    columns: [{ col: "linked_cell_type", link_to: "cell_type" }],
                },
                {
                    dataSource: "main-ds",
                    columns: [{ col: "secondary_cell_type", link_to: "cell_type" }],
                },
            ],
        });
        const linkedTarget = createDataSource("linked-target", {
            linkColumns: [
                {
                    dataSource: "main-ds",
                    columns: ["cell_type"],
                },
                {
                    dataSource: "main-ds",
                    columns: ["cell_type", "other_column"],
                },
            ],
        });
        const syncedSource = createChartManagerDataSource(syncedTarget);
        const linkedSource = createChartManagerDataSource(linkedTarget);
        const { chartManager, dialog } = setupDialog({
            mainDataSource,
            extraSources: [syncedSource, linkedSource],
        });
        dialog.applyColors("cell_type", ["#AA0000", "#00AA00"]);

        expect(mainDataSource.setColumnColors).toHaveBeenCalledWith("cell_type", ["#AA0000", "#00AA00"]);
        expect(mainDataSource.dataChanged).toHaveBeenCalledWith(["cell_type"], false);
        expect(chartManager._sync_colors).toHaveBeenCalledWith(
            [
                { col: "linked_cell_type", link_to: "cell_type" },
                { col: "secondary_cell_type", link_to: "cell_type" },
            ],
            mainDataSource,
            syncedTarget,
        );
        expect(syncedTarget.dataChanged).toHaveBeenCalledWith(["linked_cell_type", "secondary_cell_type"], false);
        expect(linkedTarget.setColumnColors).toHaveBeenCalledWith("cell_type", ["#AA0000", "#00AA00"]);
        expect(linkedTarget.dataChanged).toHaveBeenCalledWith(["cell_type"], false);
        expect(linkedTarget.setColumnColors).toHaveBeenCalledTimes(1);
    });

    test("close unmounts the React root and closes the base dialog", () => {
        const unmount = vi.fn();
        createMdvPortalMock.mockReturnValueOnce({ unmount });
        const { dialog } = setupDialog();
        dialog.close();

        expect(unmount).toHaveBeenCalledTimes(1);
        expect(baseCloseMock).toHaveBeenCalledTimes(1);
    });
});
