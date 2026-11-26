import { describe, expect, test } from "vitest";
import type { DataColumn, DataType } from "@/charts/charts";
import { processArrayBuffer } from "@/dataloaders/DataLoaders";
import {
    COLUMN_KINDS,
    ENCODINGS,
    getColumnByteLength,
    getColumnKind,
    getValuesPerRow,
    LEGACY_DATATYPES,
    STORED_TYPES,
    viewColumnStorage,
} from "@/datastore/columnTypes";
import { isDatatypeCategorical, isDatatypeNumeric } from "@/lib/utils";

const DATATYPES = Object.keys(LEGACY_DATATYPES) as DataType[];

function column(datatype: DataType, stringLength?: number) {
    return { field: datatype, name: datatype, datatype, stringLength } as DataColumn<DataType>;
}

describe("column-type contract", () => {
    // The Python package reads this copy; `pnpm vitest -u` on this file regenerates it.
    test("python/mdvtools/column_types.json is generated from LEGACY_DATATYPES", async () => {
        const contract = {
            _generated_from:
                "src/datastore/columnTypes.ts; regenerate with `pnpm vitest -u src/tests/datastore/columnTypes.spec.ts`",
            storedTypes: STORED_TYPES,
            encodings: ENCODINGS,
            kinds: COLUMN_KINDS,
            datatypes: LEGACY_DATATYPES,
        };
        await expect(`${JSON.stringify(contract, null, 4)}\n`).toMatchFileSnapshot(
            "../../../python/mdvtools/column_types.json",
        );
    });

    test("views each datatype the way the DataStore holds it", () => {
        const buffer = new ArrayBuffer(8);
        expect(viewColumnStorage(column("integer"), buffer)).toBeInstanceOf(Float32Array);
        expect(viewColumnStorage(column("double"), buffer)).toBeInstanceOf(Float32Array);
        expect(viewColumnStorage(column("int32"), buffer)).toBeInstanceOf(Int32Array);
        expect(viewColumnStorage(column("text"), buffer)).toBeInstanceOf(Uint8Array);
        expect(viewColumnStorage(column("text16"), buffer)).toBeInstanceOf(Uint16Array);
        expect(viewColumnStorage(column("unique"), buffer)).toBeInstanceOf(Uint8Array);
        expect(viewColumnStorage(column("multitext"), buffer)).toBeInstanceOf(Uint16Array);
    });

    // a /get_data response packs columns back to back, so guessing one column's width misreads the rest
    test("refuses a datatype it does not define, naming the column", () => {
        const unknown = column("boolean" as DataType);
        const message = 'Column "boolean" has datatype "boolean", which the column-type contract does not define';
        expect(() => getValuesPerRow(unknown)).toThrow(message);
        expect(() => getColumnByteLength(unknown, 2)).toThrow(message);
        expect(() => viewColumnStorage(unknown, new ArrayBuffer(8))).toThrow(message);
        expect(() => processArrayBuffer(new ArrayBuffer(8), [unknown], 2)).toThrow(message);
    });

    test("counts values per row and bytes per column", () => {
        expect(getValuesPerRow(column("double"))).toBe(1);
        expect(getValuesPerRow(column("unique", 5))).toBe(5);
        expect(getValuesPerRow(column("multitext", 3))).toBe(3);
        // multitext saved by the client has no stringLength: one code per row
        expect(getValuesPerRow(column("multitext"))).toBe(1);
        // a unique column of empty strings is really zero bytes wide
        expect(getValuesPerRow(column("unique", 0))).toBe(0);
        expect(getColumnByteLength(column("unique", 0), 10)).toBe(0);

        expect(getColumnByteLength(column("integer"), 10)).toBe(40);
        expect(getColumnByteLength(column("text"), 10)).toBe(10);
        expect(getColumnByteLength(column("text16"), 10)).toBe(20);
        expect(getColumnByteLength(column("unique", 5), 10)).toBe(50);
        expect(getColumnByteLength(column("multitext", 3), 10)).toBe(60);
    });

    test("kinds agree with the existing numeric and categorical predicates", () => {
        for (const datatype of DATATYPES) {
            const kind = getColumnKind(datatype);
            expect(isDatatypeNumeric(datatype)).toBe(kind === "numeric");
            expect(isDatatypeCategorical(datatype)).toBe(kind === "categorical" || kind === "multi_valued");
        }
    });
});

describe("processArrayBuffer", () => {
    // columns in ChartManager's order: 4-byte types first, then 2-byte, then 1-byte
    test("splits a /get_data response into one buffer per column", () => {
        const size = 2;
        const columns = [
            column("double"),
            column("int32"),
            column("multitext", 2),
            column("text16"),
            column("text"),
            column("unique", 3),
        ];
        const parts = [
            new Float32Array([1.5, -2]),
            new Int32Array([7, -8]),
            new Uint16Array([0, 65535, 1, 0]),
            new Uint16Array([300, 2]),
            new Uint8Array([1, 0]),
            new TextEncoder().encode("abcde\0"),
        ];
        const total = parts.reduce((n, p) => n + p.byteLength, 0);
        const body = new Uint8Array(total);
        let offset = 0;
        for (const p of parts) {
            body.set(new Uint8Array(p.buffer, p.byteOffset, p.byteLength), offset);
            offset += p.byteLength;
        }

        const result = processArrayBuffer(body.buffer, columns, size);

        expect(result.map((r) => r.field)).toEqual(columns.map((c) => c.field));
        result.forEach((r, i) => {
            const view = viewColumnStorage(columns[i], r.data as SharedArrayBuffer);
            expect(Array.from(view)).toEqual(Array.from(parts[i]));
        });
    });
});
