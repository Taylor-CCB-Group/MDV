import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";

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
    return {
        dataChanged: vi.fn(),
        getColumnColors: vi.fn(() => ["#111111", "#222222"]),
        getColumnList: vi.fn(() => [{ field: "cell_type", name: "Cell type" }]),
        getColumnValues: vi.fn(() => ["T cell", "B cell"]),
        linkColumns: overrides.linkColumns ?? [],
        name,
        setColumnColors: vi.fn(),
        syncColumnColors: overrides.syncColumnColors ?? [],
    };
}

describe("ColorPalette wrapper", () => {
    test("mounts the React component with the selected datasource", () => {
        const mainDataSource = createDataSource("main-ds");
        const chartManager = {
            _sync_colors: vi.fn(),
            dataSources: [{ dataStore: mainDataSource }],
            getDataSource: vi.fn(() => mainDataSource),
        };
        const ds = { dataStore: { name: "main-ds" }, name: "My Dataset" };

        const dialog = new ColorPalette(chartManager, ds);

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
            ],
        });
        const linkedTarget = createDataSource("linked-target", {
            linkColumns: [
                {
                    dataSource: "main-ds",
                    columns: ["cell_type"],
                },
            ],
        });

        const chartManager = {
            _sync_colors: vi.fn(),
            dataSources: [{ dataStore: mainDataSource }, { dataStore: syncedTarget }, { dataStore: linkedTarget }],
            getDataSource: vi.fn(() => mainDataSource),
        };
        const ds = { dataStore: { name: "main-ds" }, name: "My Dataset" };

        const dialog = new ColorPalette(chartManager, ds);
        dialog.applyColors("cell_type", ["#AA0000", "#00AA00"]);

        expect(mainDataSource.setColumnColors).toHaveBeenCalledWith("cell_type", ["#AA0000", "#00AA00"]);
        expect(mainDataSource.dataChanged).toHaveBeenCalledWith(["cell_type"], false, false);
        expect(chartManager._sync_colors).toHaveBeenCalledWith(
            [{ col: "linked_cell_type", link_to: "cell_type" }],
            mainDataSource,
            syncedTarget,
        );
        expect(syncedTarget.dataChanged).toHaveBeenCalledWith(["linked_cell_type"], false, false);
        expect(linkedTarget.setColumnColors).toHaveBeenCalledWith("cell_type", ["#AA0000", "#00AA00"]);
        expect(linkedTarget.dataChanged).toHaveBeenCalledWith(["cell_type"], false, false);
    });

    test("close unmounts the React root and closes the base dialog", () => {
        const unmount = vi.fn();
        createMdvPortalMock.mockReturnValueOnce({ unmount });

        const mainDataSource = createDataSource("main-ds");
        const chartManager = {
            _sync_colors: vi.fn(),
            dataSources: [{ dataStore: mainDataSource }],
            getDataSource: vi.fn(() => mainDataSource),
        };
        const ds = { dataStore: { name: "main-ds" }, name: "My Dataset" };

        const dialog = new ColorPalette(chartManager, ds);
        dialog.close();

        expect(unmount).toHaveBeenCalledTimes(1);
        expect(baseCloseMock).toHaveBeenCalledTimes(1);
    });
});
