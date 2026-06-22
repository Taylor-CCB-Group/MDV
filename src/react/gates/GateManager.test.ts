import { describe, expect, test, vi } from "vitest";
import type DataStore from "@/datastore/DataStore";
import { GateManager } from "./GateManager";
import type { Gate } from "./types";

vi.mock("@/dataloaders/DataLoaderUtil", () => ({
    loadColumn: vi.fn(),
}));

const GATES_COLUMN_NAME = "__gates__";
const GATE_CAPACITY = 24;

type FakeColumn = {
    field: string;
    name: string;
    datatype: string;
    data?: Float32Array | Uint16Array;
    values?: string[];
    stringLength?: number;
    delimiter?: string;
    editable?: boolean;
};

function createFakeDataStore(size = 3) {
    const columnIndex: Record<string, FakeColumn> = {
        x: {
            field: "x",
            name: "x",
            datatype: "double",
            data: new Float32Array([0.5, 2, 0.25]),
        },
        y: {
            field: "y",
            name: "y",
            datatype: "double",
            data: new Float32Array([0.5, 2, 0.25]),
        },
    };

    const dataStore = {
        size,
        name: "cells",
        config: {},
        gates: undefined,
        columnIndex,
        dirtyMetadata: new Set<string>(),
        dirtyColumns: {
            added: {},
            removed: {},
            data_changed: {},
            colors_changed: {},
        },
        addColumn: vi.fn((column: FakeColumn, data: SharedArrayBuffer | Float32Array | Uint16Array) => {
            dataStore.columnIndex[column.field] = {
                ...column,
                data: data instanceof SharedArrayBuffer ? new Uint16Array(data) : data,
            };
        }),
        dataChanged: vi.fn(),
    };

    return dataStore as unknown as DataStore;
}

function createGate(): Gate {
    return {
        id: "gate-1",
        name: "Gate A",
        columns: ["x", "y"],
        createdAt: 1,
        labelPosition: [0.5, 0.5],
        geometry: {
            type: "FeatureCollection",
            features: [
                {
                    type: "Feature",
                    properties: {},
                    geometry: {
                        type: "Polygon",
                        coordinates: [
                            [
                                [0, 0],
                                [1, 0],
                                [1, 1],
                                [0, 1],
                                [0, 0],
                            ],
                        ],
                    },
                },
            ],
        },
    };
}

describe("GateManager", () => {
    test("does not create a dense gates column during construction when there are no gates", () => {
        const dataStore = createFakeDataStore(10_000_000);

        new GateManager(dataStore);

        expect(dataStore.addColumn).not.toHaveBeenCalled();
        expect(dataStore.columnIndex[GATES_COLUMN_NAME]).toBeUndefined();
    });

    test("creates the gates column on demand when adding the first gate", async () => {
        const dataStore = createFakeDataStore();
        const gateManager = new GateManager(dataStore);

        await gateManager.addGate(createGate());

        expect(dataStore.addColumn).toHaveBeenCalledOnce();
        const gateColumn = dataStore.columnIndex[GATES_COLUMN_NAME];
        expect(gateColumn?.datatype).toBe("multitext");
        expect(gateColumn?.data).toBeInstanceOf(Uint16Array);
        expect(gateColumn?.data).toHaveLength(dataStore.size * GATE_CAPACITY);
        expect(gateColumn?.values).toEqual(["N/A", "Gate A"]);
    });
});
