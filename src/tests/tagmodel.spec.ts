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

function createMockDataStore(columnIndex: Record<string, any>, size: number) {
    return {
        size,
        name: "test-store",
        columnIndex,
        filterArray: new Uint8Array(size),
        filterSize: size,
        addListener: vi.fn(),
        addColumn: vi.fn((spec, data) => {
            columnIndex[spec.field] = {
                ...spec,
                data: new Uint16Array(data),
                buffer: data,
            };
        }),
        dataChanged: vi.fn(),
    };
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
});
