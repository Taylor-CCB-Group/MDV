import { loadColumn } from "@/dataloaders/DataLoaderUtil";
import TagModel from "@/table/TagModel";
import { DataModel } from "@/table/DataModel";
import { beforeEach, describe, expect, test, vi } from "vitest";

vi.mock("@/table/DataModel", () => ({
    DataModel: vi.fn(),
}));
vi.mock("@/dataloaders/DataLoaderUtil", () => ({
    loadColumn: vi.fn(),
}));

function createMockDataModel(selection: number[]) {
    return {
        updateModel: vi.fn(),
        addListener: vi.fn(),
        data: new Int32Array(selection),
        dataStore: undefined,
    } as any;
}

function createMockDataStore(
    columnIndex: Record<string, any>,
    size: number,
    options?: { isFiltered?: () => boolean },
) {
    const listeners: Record<string, (type: string, data?: any) => void> = {};
    const store = {
        size,
        name: "test-store",
        columnIndex,
        filterArray: new Uint8Array(size),
        filterSize: size,
        highightedData: [] as number[],
        ...(options?.isFiltered ? { isFiltered: options.isFiltered } : {}),
        addListener: vi.fn((id, listener) => {
            listeners[id] = listener;
        }),
        removeListener: vi.fn((id) => {
            delete listeners[id];
        }),
        addColumn: vi.fn((spec, data) => {
            columnIndex[spec.field] = {
                ...spec,
                data: new Uint16Array(data),
                buffer: data,
            };
        }),
        getHighlightedData: vi.fn(() => store.highightedData),
        emit(type: string, data?: any) {
            Object.values(listeners).forEach((listener) => listener(type, data));
        },
        dataChanged: vi.fn(),
    };
    return store;
}

describe("TagModel", () => {
    beforeEach(() => {
        vi.clearAllMocks();
    });

    test("creates an empty __tags column with delimiter metadata and empty sentinels", async () => {
        const columnIndex: Record<string, any> = {};
        const dataStore = createMockDataStore(columnIndex, 4);
        const dataModel = createMockDataModel([0, 1, 2, 3]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });

        const tagModel = await TagModel.create(dataStore as never);

        expect(tagModel).toBeInstanceOf(TagModel);
        expect(dataStore.addColumn).toHaveBeenCalledTimes(1);
        expect(dataStore.addColumn).toHaveBeenCalledWith(
            expect.objectContaining({
                name: "__tags",
                field: "__tags",
                datatype: "multitext",
                delimiter: ";",
                stringLength: 1,
            }),
            expect.any(SharedArrayBuffer),
        );
        expect(Array.from(columnIndex.__tags.data)).toEqual([
            65535, 65535, 65535, 65535,
        ]);
    });

    test("getTags ignores undefined legacy values", async () => {
        const mockColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a", "b", undefined] as unknown as string[],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([0, 1, 2]),
            buffer: new SharedArrayBuffer(6),
        };
        const dataStore = createMockDataStore({ __tags: mockColumn }, 3);
        const dataModel = createMockDataModel([0, 1, 2]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(mockColumn as never);

        const tagModel = await TagModel.create(dataStore as never);

        expect(tagModel.getTags()).toEqual(new Set(["a", "b"]));
    });

    test("reads packed gates-style multitext columns as item-level tags", async () => {
        const gateColumn = {
            name: "__gates__",
            field: "__gates__",
            datatype: "multitext" as const,
            values: ["N/A", "a", "b"],
            delimiter: ",",
            stringLength: 3,
            data: new Uint16Array([
                0, 65535, 65535,
                2, 65535, 65535,
                1, 2, 65535,
            ]),
            buffer: new SharedArrayBuffer(18),
        };
        const dataStore = createMockDataStore({ __gates__: gateColumn }, 3);
        const dataModel = createMockDataModel([0, 1, 2]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(gateColumn as never);

        const tagModel = await TagModel.create(dataStore as never, "__gates__");

        expect(tagModel.getTags()).toEqual(new Set(["a", "b"]));
        expect(tagModel.getTagsInSelection()).toEqual(new Set(["a", "b"]));
    });

    test("setTag updates packed multitext rows without keeping the gates placeholder", async () => {
        const gateColumn = {
            name: "__gates__",
            field: "__gates__",
            datatype: "multitext" as const,
            values: ["N/A", "a", "b"],
            delimiter: ",",
            stringLength: 3,
            data: new Uint16Array([
                0, 65535, 65535,
                2, 65535, 65535,
                1, 2, 65535,
            ]),
            buffer: new SharedArrayBuffer(18),
        };
        const dataStore = createMockDataStore({ __gates__: gateColumn }, 3);
        const dataModel = createMockDataModel([0, 1]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(gateColumn as never);

        const tagModel = await TagModel.create(dataStore as never, "__gates__");
        tagModel.setTag("a", true);

        expect(Array.from(gateColumn.data)).toEqual([
            1, 65535, 65535,
            1, 2, 65535,
            1, 2, 65535,
        ]);
        expect(dataStore.dataChanged).toHaveBeenCalledWith(["__gates__"]);
    });

    test("can operate on highlighted rows instead of the filtered model selection", async () => {
        const tagColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a", "b"],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([0, 1, 65535]),
            buffer: new SharedArrayBuffer(6),
        };
        const dataStore = createMockDataStore({ __tags: tagColumn }, 3);
        dataStore.highightedData = [1, 2];
        const dataModel = createMockDataModel([0, 1, 2]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(tagColumn as never);

        const tagModel = await TagModel.create(
            dataStore as never,
            "__tags",
            "highlighted",
        );

        expect(tagModel.getSelectionLength()).toBe(2);
        expect(tagModel.getTagsInSelection()).toEqual(new Set(["b"]));

        tagModel.setTag("a", true);

        expect(tagColumn.values).toEqual(["a", "b", "a; b"]);
        expect(Array.from(tagColumn.data)).toEqual([0, 2, 0]);
        expect(dataStore.dataChanged).toHaveBeenCalledWith(["__tags"]);
    });

    test("updates highlighted-scope selection when highlighted rows change", async () => {
        const tagColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a", "b"],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([0, 1, 65535]),
            buffer: new SharedArrayBuffer(6),
        };
        const dataStore = createMockDataStore({ __tags: tagColumn }, 3);
        dataStore.highightedData = [1];
        const dataModel = createMockDataModel([0, 1, 2]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(tagColumn as never);

        const tagModel = await TagModel.create(
            dataStore as never,
            "__tags",
            "highlighted",
        );
        const listener = vi.fn();
        tagModel.addListener(listener);

        expect(tagModel.getTagsInSelection()).toEqual(new Set(["b"]));

        dataStore.highightedData = [0, 2];
        dataStore.emit("data_highlighted", { indexes: [0, 2] });

        expect(listener).toHaveBeenCalled();
        expect(tagModel.getTagsInSelection()).toEqual(new Set(["a"]));
    });

    test("refreshes when the tag column changes in place", async () => {
        const tagColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a"],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([0, 65535]),
            buffer: new SharedArrayBuffer(4),
        };
        const dataStore = createMockDataStore({ __tags: tagColumn }, 2);
        const dataModel = createMockDataModel([0, 1]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(tagColumn as never);

        const tagModel = await TagModel.create(dataStore as never, "__tags");
        const listener = vi.fn();
        tagModel.addListener(listener);

        tagColumn.values.push("a; b", "b");
        tagColumn.data[1] = 1;
        dataStore.emit("data_changed", { columns: ["__tags"] });

        expect(dataModel.updateModel).toHaveBeenCalled();
        expect(listener).toHaveBeenCalled();
        expect(tagModel.getTags()).toEqual(new Set(["a", "b"]));
        expect(tagModel.getTagsInSelection()).toEqual(new Set(["a", "b"]));
    });

    test("matches multitext rows by individual tag item within the active scope", async () => {
        const tagColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a", "b", "a; b", "c"],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([0, 2, 1, 3]),
            buffer: new SharedArrayBuffer(8),
        };
        const dataStore = createMockDataStore({ __tags: tagColumn }, 4);
        const dataModel = createMockDataModel([0, 1, 2, 3]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(tagColumn as never);

        const filteredTagModel = await TagModel.create(dataStore as never, "__tags");
        const highlightedTagModel = await TagModel.create(
            dataStore as never,
            "__tags",
            "highlighted",
        );
        dataStore.highightedData = [1, 3];

        expect(filteredTagModel.getMatchingRowIndices("b")).toEqual([1, 2]);
        expect(highlightedTagModel.getMatchingRowIndices("b")).toEqual([1]);
    });

    test("getAnnotationViewState mirrors tags, selection tags, and length", async () => {
        const tagColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a", "b"],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([0, 1]),
            buffer: new SharedArrayBuffer(4),
        };
        const dataStore = createMockDataStore({ __tags: tagColumn }, 2);
        const dataModel = createMockDataModel([0]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(tagColumn as never);

        const tagModel = await TagModel.create(dataStore as never, "__tags");
        const view = tagModel.getAnnotationViewState();

        expect(view.selectionCount).toBe(1);
        expect(view.tagList).toEqual(new Set(["a", "b"]));
        expect(view.tagsInSelection).toEqual(new Set(["a"]));
    });

    test("getAnnotationWarningFlags: no rows shows no-selection only", async () => {
        const tagColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a"],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([65535, 65535]),
            buffer: new SharedArrayBuffer(4),
        };
        const dataStore = createMockDataStore({ __tags: tagColumn }, 2);
        const dataModel = createMockDataModel([]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(tagColumn as never);

        const tagModel = await TagModel.create(dataStore as never, "__tags");
        const flags = tagModel.getAnnotationWarningFlags();

        expect(flags.showNoSelectionWarning).toBe(true);
        expect(flags.showFilteredScopeCoversWholeTableWarning).toBe(false);
    });

    test("getAnnotationWarningFlags: filtered scope with no active filter warns whole table", async () => {
        const tagColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a"],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([0]),
            buffer: new SharedArrayBuffer(2),
        };
        const dataStore = createMockDataStore({ __tags: tagColumn }, 1, {
            isFiltered: () => false,
        });
        const dataModel = createMockDataModel([0]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(tagColumn as never);

        const tagModel = await TagModel.create(dataStore as never, "__tags");
        const flags = tagModel.getAnnotationWarningFlags();

        expect(flags.showNoSelectionWarning).toBe(false);
        expect(flags.showFilteredScopeCoversWholeTableWarning).toBe(true);
    });

    test("getAnnotationWarningFlags: highlighted scope never warns whole table", async () => {
        const tagColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a"],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([0]),
            buffer: new SharedArrayBuffer(2),
        };
        const dataStore = createMockDataStore({ __tags: tagColumn }, 1, {
            isFiltered: () => false,
        });
        const dataModel = createMockDataModel([0]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(tagColumn as never);

        const tagModel = await TagModel.create(
            dataStore as never,
            "__tags",
            "highlighted",
        );
        dataStore.highightedData = [0];

        const flags = tagModel.getAnnotationWarningFlags();

        expect(flags.showNoSelectionWarning).toBe(false);
        expect(flags.showFilteredScopeCoversWholeTableWarning).toBe(false);
    });

    test("rejects read-only multitext columns", async () => {
        const tagColumn = {
            name: "__tags",
            field: "__tags",
            datatype: "multitext" as const,
            editable: false,
            values: ["a"],
            delimiter: ";",
            stringLength: 1,
            data: new Uint16Array([0]),
            buffer: new SharedArrayBuffer(2),
        };
        const dataStore = createMockDataStore({ __tags: tagColumn }, 1);
        const dataModel = createMockDataModel([0]);
        dataModel.dataStore = dataStore;
        vi.mocked(DataModel).mockImplementation(function MockDataModel() {
            return dataModel as never;
        });
        vi.mocked(loadColumn).mockResolvedValue(tagColumn as never);

        await expect(TagModel.create(dataStore as never)).rejects.toThrow(
            "column '__tags' is not editable",
        );
    });
});
