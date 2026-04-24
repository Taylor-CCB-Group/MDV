import {
    getOwnerVisibleRows,
    isRowFilteredByOtherOwner,
} from "@/lib/filterOwnership";
import { describe, expect, test } from "vitest";

function toRows(rows: Uint32Array) {
    return Array.from(rows);
}

describe("filter ownership", () => {
    test("keeps aggregate rows when there are no local filters", () => {
        const aggregateRows = new Uint32Array([0, 1, 2]);
        const globalFilter = new Uint8Array([0, 0, 0]);

        const ownerVisibleRows = getOwnerVisibleRows({
            aggregateRows,
            globalFilter,
        });

        expect(ownerVisibleRows).toBe(aggregateRows);
        expect(toRows(ownerVisibleRows)).toEqual([0, 1, 2]);
    });

    test("treats global-only filters as external", () => {
        const globalFilter = new Uint8Array([0, 1, 0]);

        const ownerVisibleRows = getOwnerVisibleRows({
            aggregateRows: new Uint32Array([0, 2]),
            globalFilter,
        });

        expect(toRows(ownerVisibleRows)).toEqual([0, 2]);
        expect(isRowFilteredByOtherOwner(1, globalFilter)).toBe(true);
    });

    test("keeps locally filtered rows visible to the owning chart", () => {
        const globalFilter = new Uint8Array([0, 1, 0]);
        const localFilter = new Uint8Array([0, 1, 0]);

        const ownerVisibleRows = getOwnerVisibleRows({
            aggregateRows: new Uint32Array([0, 2]),
            globalFilter,
            localFilter,
            localFilterIsActive: true,
        });

        expect(toRows(ownerVisibleRows)).toEqual([0, 1, 2]);
        expect(isRowFilteredByOtherOwner(1, globalFilter, localFilter)).toBe(false);
    });

    test("hides rows filtered by both the owner and another filter", () => {
        const globalFilter = new Uint8Array([0, 2, 1, 1]);
        const localFilter = new Uint8Array([0, 1, 1, 0]);

        const ownerVisibleRows = getOwnerVisibleRows({
            aggregateRows: new Uint32Array([0]),
            globalFilter,
            localFilter,
            localFilterIsActive: true,
        });

        expect(toRows(ownerVisibleRows)).toEqual([0, 2]);
        expect(isRowFilteredByOtherOwner(1, globalFilter, localFilter)).toBe(true);
        expect(isRowFilteredByOtherOwner(3, globalFilter, localFilter)).toBe(true);
    });

    test("treats local filter value 2 as chart-local exclusion", () => {
        const globalFilter = new Uint8Array([0, 2, 0]);
        const localFilter = new Uint8Array([0, 2, 0]);

        const ownerVisibleRows = getOwnerVisibleRows({
            aggregateRows: new Uint32Array([0, 2]),
            globalFilter,
            localFilter,
            localFilterIsActive: true,
        });

        expect(toRows(ownerVisibleRows)).toEqual([0, 2]);
        expect(isRowFilteredByOtherOwner(1, globalFilter, localFilter)).toBe(false);
    });

    test("applies chart-scope predicates before ownership checks", () => {
        const globalFilter = new Uint8Array([0, 1, 0, 1]);
        const localFilter = new Uint8Array([0, 1, 0, 1]);
        const scopePredicate = (rowIndex: number) => rowIndex !== 1 && rowIndex !== 2;

        const ownerVisibleRows = getOwnerVisibleRows({
            aggregateRows: new Uint32Array([0, 2]),
            globalFilter,
            localFilter,
            localFilterIsActive: true,
            scopePredicate,
        });

        expect(toRows(ownerVisibleRows)).toEqual([0, 3]);
        expect(isRowFilteredByOtherOwner(1, globalFilter, localFilter, scopePredicate)).toBe(false);
    });
});
