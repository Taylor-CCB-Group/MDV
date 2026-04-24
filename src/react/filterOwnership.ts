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

export function getOwnerVisibleRows({
    aggregateRows,
    globalFilter,
    localFilter,
    localFilterIsActive = false,
    scopePredicate,
}: OwnerVisibleRowsInput) {
    if (!localFilterIsActive) {
        return aggregateRows;
    }

    const rows: number[] = [];
    for (let rowIndex = 0; rowIndex < globalFilter.length; rowIndex += 1) {
        if (
            isRowVisibleToFilterOwner(
                rowIndex,
                globalFilter,
                localFilter,
                scopePredicate,
            )
        ) {
            rows.push(rowIndex);
        }
    }
    return Uint32Array.from(rows);
}
