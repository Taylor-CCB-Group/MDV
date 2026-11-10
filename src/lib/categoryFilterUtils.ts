import type { CategoricalDataType } from "@/charts/charts";
import {
    getMultitextCapacity,
    getMultitextValueItems,
    inferMultitextEmptyItem,
    MULTITEXT_EMPTY_INDEX,
} from "@/lib/multitext";
import { isArray } from "@/lib/utils";

type CategoryFilterColumn = {
    data?: ArrayLike<number>;
    datatype: CategoricalDataType;
    values: string[];
    delimiter?: string;
    stringLength?: number;
};

export function getCategoryFilterIndices(
    data: number[] | Uint32Array,
    column?: CategoryFilterColumn,
    category?: string | string[] | null,
) {
    if (!column || !category || !column.data) {
        return [];
    }
    const matches = compileCategoryFilter(column, category);
    return data.filter((rowIndex) => matches(rowIndex));
}

/**
 * Compiles a category filter into a per-row test. The selected values are copied and resolved
 * to a lookup over the column's value indices once, so each row costs a table lookup (a few,
 * for multitext).
 *
 * Callers typically pass `category` straight from chart config, which is MobX-observable:
 * re-deriving it per row put proxy reads and string work in every iteration of whole-dataset
 * scans. Changes to `category` after compiling are not seen; compile again.
 *
 * Column edits are seen: the lookup is rebuilt when the column's data or values array is
 * replaced (a new multitext layout always reallocates data) or values are appended (new gates
 * or tags, which can give a multitext row a new combined value). That check costs a few ms per
 * million rows; inlining the lookup into the same closure measured slower.
 */
export function compileCategoryFilter(
    column: CategoryFilterColumn,
    category: string | string[] | null,
): (rowIndex: number) => boolean {
    const selectedValues = (isArray(category) ? [...category] : [category]).filter(
        (value): value is string => typeof value === "string" && value.length > 0,
    );
    if (selectedValues.length === 0) {
        return () => false;
    }
    const selected = new Set(selectedValues);

    let matches = compileRowMatcher(column, selectedValues, selected);
    let data = column.data;
    let values = column.values;
    let valueCount = values.length;
    return (rowIndex) => {
        if (column.data !== data || column.values !== values || column.values.length !== valueCount) {
            matches = compileRowMatcher(column, selectedValues, selected);
            data = column.data;
            values = column.values;
            valueCount = values.length;
        }
        return matches(rowIndex);
    };
}

/** The lookup for the column as it is now; see {@link compileCategoryFilter}. */
function compileRowMatcher(
    column: CategoryFilterColumn,
    selectedValues: string[],
    selected: Set<string>,
): (rowIndex: number) => boolean {
    const data = column.data;
    if (!data) {
        return () => false;
    }
    const values = column.values;
    const matches = new Uint8Array(values.length);
    for (const value of selectedValues) {
        const valueIndex = values.indexOf(value);
        if (valueIndex !== -1) matches[valueIndex] = 1;
    }

    if (column.datatype !== "multitext") {
        return (rowIndex) => matches[data[rowIndex]] === 1;
    }

    // a multitext row also matches if any item of one of its values is selected
    const emptyItem = inferMultitextEmptyItem(column);
    for (let valueIndex = 0; valueIndex < values.length; valueIndex++) {
        if (getMultitextValueItems(column, valueIndex, { emptyItem }).some((item) => selected.has(item))) {
            matches[valueIndex] = 1;
        }
    }
    const capacity = getMultitextCapacity(column);
    return (rowIndex) => {
        const start = rowIndex * capacity;
        for (let offset = 0; offset < capacity; offset++) {
            const valueIndex = data[start + offset];
            if (valueIndex == null || valueIndex === MULTITEXT_EMPTY_INDEX) return false;
            if (matches[valueIndex] === 1) return true;
        }
        return false;
    };
}

export type { CategoryFilterColumn };
