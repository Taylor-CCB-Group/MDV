import { describe, expect, test } from "vitest";
import { type SortColumn, sortOrder } from "@/datastore/sortColumns";

const EMPTY = 65535;

function sorted(size: number, columns: SortColumn[]) {
    const order = Uint32Array.from({ length: size }, (_, i) => i);
    sortOrder(order, columns);
    return Array.from(order);
}

function uniqueBuffer(strings: string[], stringLength: number) {
    const bytes = new Uint8Array(strings.length * stringLength);
    const encoder = new TextEncoder();
    strings.forEach((s, i) => {
        bytes.set(encoder.encode(s), i * stringLength);
    });
    return bytes.buffer;
}

describe("sortOrder", () => {
    // `integer` is held as float32, like `double`; reading those bytes as int32 reverses negatives
    test("sorts integer columns, including negative values", () => {
        const buffer = new Float32Array([3, -1, -2, Number.NaN, 0]).buffer;
        expect(sorted(5, [{ datatype: "integer", buffer }])).toEqual([2, 1, 4, 0, 3]);
    });

    test("sorts missing values last in both directions, after Infinity", () => {
        const buffer = new Float32Array([Number.POSITIVE_INFINITY, Number.NaN, 1, Number.NaN]).buffer;
        expect(sorted(4, [{ datatype: "double", buffer }])).toEqual([2, 0, 1, 3]);
        expect(sorted(4, [{ datatype: "double", buffer, desc: true }])).toEqual([0, 2, 1, 3]);
    });

    test("sorts int32 columns", () => {
        const buffer = new Int32Array([5, -7, 2]).buffer;
        expect(sorted(3, [{ datatype: "int32", buffer }])).toEqual([1, 2, 0]);
    });

    test("sorts double columns descending", () => {
        const buffer = new Float32Array([0.5, 2.5, -1.5]).buffer;
        expect(sorted(3, [{ datatype: "double", buffer, desc: true }])).toEqual([1, 0, 2]);
    });

    test("sorts text and text16 columns by value, not code", () => {
        const values = ["pear", "apple", "fig"];
        const text = new Uint8Array([0, 1, 2]).buffer;
        const text16 = new Uint16Array([0, 1, 2]).buffer;
        expect(sorted(3, [{ datatype: "text", buffer: text, values }])).toEqual([1, 2, 0]);
        expect(sorted(3, [{ datatype: "text16", buffer: text16, values, desc: true }])).toEqual([0, 2, 1]);
    });

    test("sorts unique columns, ignoring padding", () => {
        const buffer = uniqueBuffer(["bob", "al", "cat"], 3);
        expect(sorted(3, [{ datatype: "unique", buffer, stringLength: 3 }])).toEqual([1, 0, 2]);
    });

    // DataStore gives a unique column of empty strings stringLength 0 and a zero-byte buffer
    test("treats a zero-width unique column as all equal, leaving ties to later columns", () => {
        const columns: SortColumn[] = [
            { datatype: "unique", buffer: new ArrayBuffer(0), stringLength: 0 },
            { datatype: "double", buffer: new Float32Array([3, 1, 2]).buffer },
        ];
        expect(sorted(3, columns)).toEqual([1, 2, 0]);
    });

    test("sorts packed multitext columns item by item, shorter first", () => {
        const values = ["b", "a"];
        // rows: [a], [b], [a, b], []
        const buffer = new Uint16Array([1, EMPTY, 0, EMPTY, 1, 0, EMPTY, EMPTY]).buffer;
        const column: SortColumn = { datatype: "multitext", buffer, values, stringLength: 2 };
        expect(sorted(4, [column])).toEqual([3, 0, 2, 1]);
        expect(sorted(4, [{ ...column, desc: true }])).toEqual([1, 2, 0, 3]);
    });

    test("sorts multitext columns saved without a stringLength", () => {
        const values = ["b", "a"];
        const buffer = new Uint16Array([0, 1, EMPTY]).buffer;
        expect(sorted(3, [{ datatype: "multitext", buffer, values }])).toEqual([2, 1, 0]);
    });

    test("breaks ties with later columns", () => {
        const first = new Uint8Array([1, 0, 1, 0]).buffer;
        const second = new Float32Array([4, 3, 2, 1]).buffer;
        const columns: SortColumn[] = [
            { datatype: "text", buffer: first, values: ["x", "y"] },
            { datatype: "double", buffer: second },
        ];
        expect(sorted(4, columns)).toEqual([3, 1, 2, 0]);
    });
});
