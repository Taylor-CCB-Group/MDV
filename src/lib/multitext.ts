export const MULTITEXT_EMPTY_INDEX = 65535;

export type MultiTextColumnLike = {
    delimiter?: string;
    field?: string;
    name?: string;
    stringLength?: number;
    values: Array<string | undefined>;
    data?: ArrayLike<number>;
};

type MultiTextItemOptions = {
    emptyItem?: string;
};

function trimNonEmpty(items: Iterable<string | null | undefined>) {
    return Array.from(items)
        .flatMap((item) => (typeof item === "string" ? [item.trim()] : []))
        .filter((item) => item.length > 0);
}

function filterSemanticItems(items: string[], options?: MultiTextItemOptions) {
    if (!options?.emptyItem) {
        return items;
    }
    return items.filter((item) => item !== options.emptyItem);
}

export function getMultitextDelimiter(
    column: Pick<MultiTextColumnLike, "delimiter" | "values">,
) {
    if (column.delimiter) {
        return column.delimiter;
    }
    return column.values.some((value) => value?.includes(";")) ? ";" : ",";
}

export function getMultitextJoinDelimiter(delimiter: string) {
    return delimiter.endsWith(" ") ? delimiter : `${delimiter} `;
}

export function splitMultitextItems(value: string | undefined, delimiter = ",") {
    if (!value) {
        return [];
    }
    return trimNonEmpty(value.split(delimiter));
}

export function joinMultitextItems(items: Iterable<string>, delimiter = ",") {
    return trimNonEmpty(items).join(getMultitextJoinDelimiter(delimiter));
}

export function normalizeMultitextItems(items: Iterable<string>) {
    return Array.from(new Set(trimNonEmpty(items))).sort((left, right) =>
        left.localeCompare(right),
    );
}

export function getMultitextCapacity(
    column: Pick<MultiTextColumnLike, "stringLength">,
) {
    return Math.max(column.stringLength || 1, 1);
}

export function isPackedMultitextColumn(
    column: Pick<MultiTextColumnLike, "stringLength">,
) {
    return getMultitextCapacity(column) > 1;
}

export function inferMultitextEmptyItem(column: Pick<MultiTextColumnLike, "values" | "stringLength">) {
    if (!isPackedMultitextColumn(column)) {
        return undefined;
    }
    return column.values[0] === "N/A" ? "N/A" : undefined;
}

export function getMultitextValueItems(
    column: Pick<MultiTextColumnLike, "delimiter" | "stringLength" | "values">,
    valueIndex: number,
    options?: MultiTextItemOptions,
) {
    const value = column.values[valueIndex];
    if (value == null) {
        return [];
    }
    const items = isPackedMultitextColumn(column)
        ? [value]
        : splitMultitextItems(value, getMultitextDelimiter(column));
    return filterSemanticItems(items, options);
}

export function getMultitextItemValues(
    column: Pick<MultiTextColumnLike, "delimiter" | "stringLength" | "values">,
    options?: MultiTextItemOptions,
) {
    if (isPackedMultitextColumn(column)) {
        return filterSemanticItems(
            trimNonEmpty(column.values),
            options,
        );
    }
    const delimiter = getMultitextDelimiter(column);
    const values = column.values.flatMap((value) =>
        splitMultitextItems(value, delimiter),
    );
    return Array.from(new Set(filterSemanticItems(values, options)));
}

export function getRowMultitextValueIndices(
    column: Pick<MultiTextColumnLike, "data" | "stringLength">,
    rowIndex: number,
) {
    if (!column.data) {
        return [];
    }
    const capacity = getMultitextCapacity(column);
    const start = rowIndex * capacity;
    const indices: number[] = [];
    for (let offset = 0; offset < capacity; offset += 1) {
        const valueIndex = column.data[start + offset];
        if (
            valueIndex == null ||
            valueIndex === MULTITEXT_EMPTY_INDEX
        ) {
            break;
        }
        indices.push(valueIndex);
    }
    return indices;
}

export function getRowMultitextItems(
    column: Pick<MultiTextColumnLike, "data" | "delimiter" | "stringLength" | "values">,
    rowIndex: number,
    options?: MultiTextItemOptions,
) {
    const indices = getRowMultitextValueIndices(column, rowIndex);
    const items = indices.flatMap((valueIndex) =>
        getMultitextValueItems(column, valueIndex, options),
    );
    return filterSemanticItems(items, options);
}
