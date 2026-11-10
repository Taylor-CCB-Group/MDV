export type FilterArrayLike = ArrayLike<number>;
export type RowScopePredicate = (rowIndex: number) => boolean;

type OwnerVisibleRowsInput = {
    aggregateRows: Uint32Array;
    globalFilter: FilterArrayLike;
    localFilter?: FilterArrayLike | null;
    localFilterIsActive?: boolean;
    scopePredicate?: RowScopePredicate | null;
};

export function isRowInOwnedFilterScope(
    rowIndex: number,
    localFilter?: FilterArrayLike | null,
    scopePredicate?: RowScopePredicate | null,
) {
    if (localFilter?.[rowIndex] === 2) {
        return false;
    }
    return scopePredicate ? scopePredicate(rowIndex) : true;
}

export function isRowVisibleToFilterOwner(
    rowIndex: number,
    globalFilter: FilterArrayLike,
    localFilter?: FilterArrayLike | null,
    scopePredicate?: RowScopePredicate | null,
) {
    if (!isRowInOwnedFilterScope(rowIndex, localFilter, scopePredicate)) {
        return false;
    }

    const globalValue = globalFilter[rowIndex] ?? 0;
    if (globalValue === 0) {
        return true;
    }

    return localFilter != null && globalValue === localFilter[rowIndex];
}

export function isRowFilteredByOtherOwner(
    rowIndex: number,
    globalFilter: FilterArrayLike,
    localFilter?: FilterArrayLike | null,
    scopePredicate?: RowScopePredicate | null,
) {
    if (!isRowInOwnedFilterScope(rowIndex, localFilter, scopePredicate)) {
        return false;
    }

    const globalValue = globalFilter[rowIndex] ?? 0;
    if (globalValue === 0) {
        return false;
    }

    return localFilter == null || globalValue !== localFilter[rowIndex];
}

/**
 * Rows the owning chart should draw, scanned from the current global filter (the
 * worker-computed `aggregateRows` can lag behind it). Runs over every row on each filter
 * change, so the per-row test is inlined with the cheap checks first.
 */
export function getOwnerVisibleRows({
    globalFilter,
    localFilter,
    localFilterIsActive = false,
    scopePredicate,
}: OwnerVisibleRowsInput) {
    const local = localFilterIsActive ? (localFilter ?? null) : null;
    const size = globalFilter.length;
    const rows = new Uint32Array(size);
    let count = 0;
    for (let rowIndex = 0; rowIndex < size; rowIndex++) {
        const globalValue = globalFilter[rowIndex] ?? 0;
        if (local) {
            // see isRowVisibleToFilterOwner
            const localValue = local[rowIndex];
            if (localValue === 2) continue;
            if (globalValue !== 0 && globalValue !== localValue) continue;
        } else if (globalValue !== 0) {
            continue;
        }
        if (scopePredicate && !scopePredicate(rowIndex)) continue;
        rows[count++] = rowIndex;
    }
    return rows.slice(0, count);
}
