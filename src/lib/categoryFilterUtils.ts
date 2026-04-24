import type { CategoricalDataType } from "@/charts/charts";
import {
    getRowMultitextItems,
    getRowMultitextValueIndices,
    inferMultitextEmptyItem,
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
    if (!column || !category) {
        return [];
    }

    const selectedValues = (isArray(category) ? category : [category]).filter(
        (value): value is string => typeof value === "string" && value.length > 0,
    );

    if (selectedValues.length === 0 || !column.data) {
        return [];
    }

    return data.filter((rowIndex) =>
        rowMatchesCategoryFilter(rowIndex, column, selectedValues),
    );
}

export function rowMatchesCategoryFilter(
    rowIndex: number,
    column: CategoryFilterColumn,
    category: string | string[] | null,
) {
    if (!category || !column.data) {
        return false;
    }

    const selectedValues = (isArray(category) ? category : [category]).filter(
        (value): value is string => typeof value === "string" && value.length > 0,
    );

    if (selectedValues.length === 0) {
        return false;
    }

    if (column.datatype === "multitext") {
        const emptyItem = inferMultitextEmptyItem(column);
        const selectedValueIndices = new Set(
            selectedValues
                .map((value) => column.values.indexOf(value))
                .filter((valueIndex) => valueIndex !== -1),
        );
        const rowValueIndices = getRowMultitextValueIndices(
            column,
            rowIndex,
        );
        if (
            rowValueIndices.some((valueIndex) =>
                selectedValueIndices.has(valueIndex),
            )
        ) {
            return true;
        }

        const rowItems = getRowMultitextItems(column, rowIndex, {
            emptyItem,
        });
        return selectedValues.some((value) => rowItems.includes(value));
    }

    const categoryValueIndices = selectedValues
        .map((value) => column.values.indexOf(value))
        .filter((valueIndex) => valueIndex !== -1);

    if (categoryValueIndices.length === 0) {
        return false;
    }

    return categoryValueIndices.includes(column.data[rowIndex]);
}

export type { CategoryFilterColumn };
