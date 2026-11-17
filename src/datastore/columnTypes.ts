/**
 * The column-type contract (ADR-0024).
 *
 * A column's type has three parts:
 * - **stored type**: the zarr v3 data type its values are kept in;
 * - **encoding**: plain values, categorical codes into the column's `values`, or several codes per row;
 * - **kind**: how charts and filters treat it.
 *
 * A column's `datatype` is the legacy field. {@link LEGACY_DATATYPES} says what each legacy datatype
 * means in those terms, as h5 projects and the DataStore hold it. This is the one place that says
 * so: the Python package reads a copy generated from it (`python/mdvtools/column_types.json`, kept
 * in step by `src/tests/datastore/columnTypes.spec.ts`; update it with `pnpm vitest -u` on that file).
 */
import type { DataType } from "@/charts/charts";
import { getMultitextCapacity } from "@/lib/multitext";

/** zarr v3 data type names. `null_terminated_bytes` is a zarr-python extension without a stable spec. */
export const STORED_TYPES = [
    "bool",
    "int8",
    "int16",
    "int32",
    "int64",
    "uint8",
    "uint16",
    "uint32",
    "uint64",
    "float32",
    "float64",
    "string",
    "null_terminated_bytes",
] as const;
export type StoredType = (typeof STORED_TYPES)[number];

export const ENCODINGS = ["plain", "categorical", "multi_valued"] as const;
export type Encoding = (typeof ENCODINGS)[number];

export const COLUMN_KINDS = ["numeric", "categorical", "unique", "multi_valued"] as const;
export type ColumnKind = (typeof COLUMN_KINDS)[number];

/** How many stored values each row has: one, or the column's `stringLength` (at least one). */
export type RowWidth = "one" | "stringLength";

export type LegacyDatatype = {
    storedType: StoredType;
    encoding: Encoding;
    kind: ColumnKind;
    rowWidth: RowWidth;
};

/**
 * What each legacy `datatype` means, as stored in `datafile.h5` and held in the DataStore.
 * Note that `integer` is float32, like `double`: integers above 2^24 lose precision.
 */
export const LEGACY_DATATYPES = {
    integer: { storedType: "float32", encoding: "plain", kind: "numeric", rowWidth: "one" },
    double: { storedType: "float32", encoding: "plain", kind: "numeric", rowWidth: "one" },
    int32: { storedType: "int32", encoding: "plain", kind: "numeric", rowWidth: "one" },
    text: { storedType: "uint8", encoding: "categorical", kind: "categorical", rowWidth: "one" },
    text16: { storedType: "uint16", encoding: "categorical", kind: "categorical", rowWidth: "one" },
    unique: { storedType: "null_terminated_bytes", encoding: "plain", kind: "unique", rowWidth: "stringLength" },
    multitext: { storedType: "uint16", encoding: "multi_valued", kind: "multi_valued", rowWidth: "stringLength" },
} as const satisfies Record<DataType, LegacyDatatype>;

/** The stored types the DataStore holds legacy columns in, and the bytes per value of each. */
const HELD_TYPE_BYTES = {
    float32: 4,
    int32: 4,
    uint16: 2,
    uint8: 1,
    null_terminated_bytes: 1,
} as const satisfies Partial<Record<StoredType, number>>;

export type StorageArray = Float32Array | Int32Array | Uint16Array | Uint8Array;

/** `stringLength` is optional: multitext columns saved by the client have none. */
type ColumnShape = { datatype: DataType; stringLength?: number };

/** `undefined` for a datatype the contract does not know, such as one from a malformed config. */
export function getColumnKind(datatype: DataType): ColumnKind | undefined {
    return LEGACY_DATATYPES[datatype]?.kind;
}

/** Stored values per row: 1, or `stringLength` for `unique` (bytes) and `multitext` (codes). */
export function getValuesPerRow(column: ColumnShape): number {
    return LEGACY_DATATYPES[column.datatype].rowWidth === "one" ? 1 : getMultitextCapacity(column);
}

/** Bytes a column of `size` rows occupies in the DataStore, and in `/get_data` responses. */
export function getColumnByteLength(column: ColumnShape, size: number): number {
    const bytesPerValue = HELD_TYPE_BYTES[LEGACY_DATATYPES[column.datatype].storedType];
    return size * getValuesPerRow(column) * bytesPerValue;
}

/**
 * A typed-array view of a column's stored values, as the DataStore holds them.
 * `length` counts values, not rows (see {@link getValuesPerRow}).
 */
export function viewColumnStorage(
    datatype: DataType,
    buffer: ArrayBufferLike,
    byteOffset = 0,
    length?: number,
): StorageArray {
    switch (LEGACY_DATATYPES[datatype].storedType) {
        case "float32":
            return new Float32Array(buffer, byteOffset, length);
        case "int32":
            return new Int32Array(buffer, byteOffset, length);
        case "uint16":
            return new Uint16Array(buffer, byteOffset, length);
        case "uint8":
        case "null_terminated_bytes":
            return new Uint8Array(buffer, byteOffset, length);
    }
}
