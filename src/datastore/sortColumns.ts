import type { DataType } from "@/charts/charts";
import { getMultitextCapacity } from "@/lib/multitext";

/** What `SortableDimension` sends the sort worker for each column it sorts by. */
export type SortColumn = {
    datatype: DataType;
    buffer: ArrayBufferLike;
    values?: string[];
    stringLength?: number;
    desc?: boolean;
};

type RowComparator = (a: number, b: number) => number;

const MULTITEXT_EMPTY = 65535;

/** Sorts `order` (row indices) in place, by each column in turn. */
export function sortOrder(order: Uint32Array, columns: SortColumn[]) {
    const comparators = columns.map(getRowComparator);
    order.sort((a, b) => {
        for (const compare of comparators) {
            const r = compare(a, b);
            if (r !== 0) return r;
        }
        return 0;
    });
}

function getRowComparator(column: SortColumn): RowComparator {
    const ascending = getAscendingComparator(column);
    return column.desc ? (a, b) => -ascending(a, b) : ascending;
}

// Typed arrays here must match how `DataStore.setColumnData` holds each datatype:
// `integer` and `double` are both float32.
function getAscendingComparator(column: SortColumn): RowComparator {
    switch (column.datatype) {
        case "integer":
        case "double":
            return getNumberComparator(new Float32Array(column.buffer));
        case "int32":
            return getNumberComparator(new Int32Array(column.buffer));
        case "text":
            return getTextComparator(new Uint8Array(column.buffer), column.values ?? []);
        case "text16":
            return getTextComparator(new Uint16Array(column.buffer), column.values ?? []);
        case "unique":
            return getUniqueComparator(new Uint8Array(column.buffer), column.stringLength ?? 1);
        case "multitext":
            return getMultiTextComparator(
                new Uint16Array(column.buffer),
                column.values ?? [],
                getMultitextCapacity(column),
            );
    }
}

/** Missing (NaN) values sort as if they were the largest value. */
function getNumberComparator(data: Float32Array | Int32Array): RowComparator {
    return (a, b) => {
        const va = Number.isNaN(data[a]) ? Number.MAX_VALUE : data[a];
        const vb = Number.isNaN(data[b]) ? Number.MAX_VALUE : data[b];
        return va > vb ? 1 : va < vb ? -1 : 0;
    };
}

/** Rank of each value code when the values are sorted alphabetically. */
function getValueRanks(values: string[]) {
    const ranks = new Int32Array(values.length);
    values
        .map((value, code) => ({ value, code }))
        .sort((x, y) => x.value.localeCompare(y.value))
        .forEach(({ code }, rank) => {
            ranks[code] = rank;
        });
    return ranks;
}

function getTextComparator(data: Uint8Array | Uint16Array, values: string[]): RowComparator {
    const ranks = getValueRanks(values);
    return (a, b) => ranks[data[a]] - ranks[data[b]];
}

function getUniqueComparator(data: Uint8Array, stringLength: number): RowComparator {
    const decoder = new TextDecoder();
    const rows = Math.floor(data.length / stringLength);
    const names: string[] = new Array(rows);
    for (let row = 0; row < rows; row++) {
        const bytes = data.subarray(row * stringLength, (row + 1) * stringLength);
        // fixed-width strings are padded with zero bytes
        names[row] = decoder.decode(bytes).replace(/\0+$/, "");
    }
    return (a, b) => names[a].localeCompare(names[b]);
}

/** Compares item by item; a row that runs out of items first sorts first. */
function getMultiTextComparator(data: Uint16Array, values: string[], capacity: number): RowComparator {
    const ranks = getValueRanks(values);
    const rank = (code: number) => (code === MULTITEXT_EMPTY ? -1 : ranks[code]);
    return (a, b) => {
        for (let i = 0; i < capacity; i++) {
            const ca = data[a * capacity + i];
            const cb = data[b * capacity + i];
            if (ca === MULTITEXT_EMPTY && cb === MULTITEXT_EMPTY) return 0;
            const r = rank(ca) - rank(cb);
            if (r !== 0) return r;
        }
        return 0;
    };
}
