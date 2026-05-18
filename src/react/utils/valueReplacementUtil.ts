import type { DataType, LoadedDataColumn } from "@/charts/charts";
import {
    getMultitextCapacity,
    getMultitextDelimiter,
    getMultitextJoinDelimiter,
    splitMultitextItems,
} from "@/lib/multitext";

/**
 * Replaces all occurrences of find value with replace value in the given cell value.
 *
 * @param findValue - Find value
 * @param replaceValue - Replace value
 * @param cellValue - Cell value in the table
 * @returns Replaced cell value
 */
export const replaceValueInString = (findValue: string, replaceValue: string, cellValue: string) => {
    if (!findValue || !cellValue) {
        return cellValue;
    }
    
    // Escape special characters
    const escapedFind = findValue.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");

    // Using case-insensitive replacement
    const regex = new RegExp(escapedFind, "gi");
    return cellValue.replace(regex, replaceValue);
};

const MULTITEXT_EMPTY_VALUE = 65535;
type DecodedCategoricalValue = string | null;

/**
 * Narrows a loaded column to the categorical column shape used by table editing.
 *
 * MDV categorical columns do not store labels directly in `data`; they store numeric
 * indices into `values`. `multitext` uses the same dictionary, with `65535` as an
 * empty slot marker.
 *
 * @param column Loaded data column being edited.
 * @returns The column's categorical data buffer, values dictionary, and datatype.
 * @throws If the column is not a supported categorical datatype or has no values dictionary.
 */
function assertCategoricalColumn(column: LoadedDataColumn<DataType>) {
    const { data, values, datatype } = column;

    if (datatype !== "text" && datatype !== "text16" && datatype !== "multitext") {
        throw new Error(`Unsupported categorical datatype: ${datatype} for column: ${column.field}`);
    }

    if (!values) {
        throw new Error(`No values array found for column: ${column.field}`);
    }

    return { data, values, datatype };
}

/**
 * Converts an encoded categorical data buffer into display labels.
 *
 * The write path edits labels first because that mirrors what users see in the
 * table. After edits are applied, {@link rewriteCategoricalDataFromLabels} encodes
 * the labels back into row indices and keeps `column.values` in sync.
 *
 * @param column Loaded text, text16, or multitext column.
 * @returns One decoded label per encoded data slot. `null` represents an empty multitext slot.
 * @throws If any non-empty encoded index does not exist in `column.values`.
 */
function decodeCategoricalData(column: LoadedDataColumn<DataType>): DecodedCategoricalValue[] {
    const { data, values, datatype } = assertCategoricalColumn(column);
    const decoded: DecodedCategoricalValue[] = [];

    // Categorical row buffers store numbers. Decode them into labels so the write
    // path can reason in table-visible values, then encode back once at the end.
    for (let i = 0; i < data.length; i++) {
        const valueIndex = data[i];
        if (datatype === "multitext" && valueIndex === MULTITEXT_EMPTY_VALUE) {
            decoded.push(null);
            continue;
        }

        if (valueIndex < 0 || valueIndex >= values.length) {
            throw new Error(`Invalid value index ${valueIndex} for column: ${column.field}. Index out of bounds.`);
        }

        decoded.push(values[valueIndex] ?? "");
    }

    return decoded;
}

/**
 * Rebuilds a categorical column from decoded labels.
 *
 * This is the inverse of {@link decodeCategoricalData}. It removes unused labels
 * from `column.values`, appends new labels, remaps the encoded `column.data`
 * indices, and preserves existing category colors for labels that remain.
 *
 * @param column Loaded text, text16, or multitext column to mutate.
 * @param decodedValues Labels to encode back into `column.data`; `null` writes an empty multitext slot.
 * @param maxValues Maximum allowed dictionary size for the datatype.
 * @throws If the updated values dictionary would exceed the datatype limit.
 */
function rewriteCategoricalDataFromLabels(
    column: LoadedDataColumn<DataType>,
    decodedValues: DecodedCategoricalValue[],
    maxValues: number,
) {
    const { data, values } = assertCategoricalColumn(column);

    const usedLabels = new Set<string>();
    for (const label of decodedValues) {
        if (label !== null) {
            usedLabels.add(label);
        }
    }

    if (usedLabels.size > maxValues) {
        throw new Error(`Column exceeded ${maxValues} values while setting: ${column.field}`);
    }

    // Preserve existing dictionary order for labels that are still used. This keeps
    // existing category colors aligned. New labels are appended in row encounter order.
    const nextValues: string[] = [];
    const nextColors: typeof column.colors = column.colors ? [] : undefined;
    for (let oldIndex = 0; oldIndex < values.length; oldIndex++) {
        const label = values[oldIndex];
        if (!usedLabels.has(label)) {
            continue;
        }
        nextValues.push(label);
        if (nextColors && column.colors?.[oldIndex] !== undefined) {
            nextColors.push(column.colors[oldIndex]);
        }
    }

    const nextValueIndexes = new Map(nextValues.map((label, index) => [label, index]));
    for (const label of decodedValues) {
        if (label === null || nextValueIndexes.has(label)) {
            continue;
        }
        nextValueIndexes.set(label, nextValues.length);
        nextValues.push(label);
    }

    decodedValues.forEach((label, i) => {
        if (label === null) {
            data[i] = MULTITEXT_EMPTY_VALUE;
            return;
        }
        const nextValueIndex = nextValueIndexes.get(label);
        if (nextValueIndex === undefined) {
            throw new Error(`No encoded value found for "${label}" in column: ${column.field}`);
        }
        data[i] = nextValueIndex;
    });

    values.splice(0, values.length, ...nextValues);
    if (column.colors && nextColors) {
        column.colors.splice(0, column.colors.length, ...nextColors);
    }
}

/**
 * Gets the current cell value for the given datatype.
 *
 * @param column - Column object with metadata
 * @param dataIndex - Index of the value in the data array
 * @returns Cell value of the column in the table
 * @throws Error if column data is invalid or missing
 */
export const getCellValueAsString = (column: LoadedDataColumn<DataType>, dataIndex: number): string => {
    if (!column.data) {
        throw new Error(`No data found in column: ${column.field}`);
    }

    const { data, values, stringLength, datatype } = column;

    // text / text16
    if (datatype === "text" || datatype === "text16") {
        if (!values) {
            throw new Error(`No values array found for column: ${column.field}`);
        }

        const valueIndex = data[dataIndex];

        if (valueIndex === undefined || valueIndex === null) return "";
        if (valueIndex < 0 || valueIndex >= values.length) {
            throw new Error(`Invalid value index ${valueIndex} for column: ${column.field}. Index out of bounds.`);
        }

        return values[valueIndex] ?? "";
    }

    // numeric
    if (datatype === "double" || datatype === "int32" || datatype === "integer") {
        const value = data[dataIndex];

        if (value === undefined || value === null || Number.isNaN(value)) return "";
        return String(value);
    }

    // multitext
    if (datatype === "multitext") {
        if (!values) {
            throw new Error(`Missing values or stringLength for multitext column: ${column.field}`);
        }

        const capacity = getMultitextCapacity(column);
        const baseIndex = dataIndex * capacity;
        const cellData = data.slice(baseIndex, baseIndex + capacity);

        return Array.from(cellData)
            .filter((x) => x !== 65535) // 65535 = empty
            .map((index) => {
                if (index < 0 || index >= values.length) {
                    throw new Error(
                        `Invalid multitext index ${index} for column: ${column.field}. Index out of bounds.`,
                    );
                }
                return values[index];
            })
            .join(getMultitextJoinDelimiter(getMultitextDelimiter(column)));
    }

    // unique
    if (datatype === "unique") {
        if (!stringLength) {
            throw new Error(`Missing stringLength for column: ${column.field}`);
        }

        const baseIndex = dataIndex * stringLength;
        if (baseIndex + stringLength > data.length) {
            throw new Error(
                `Invalid index ${baseIndex + stringLength} for column: ${column.field}. Index out of bounds.`,
            );
        }

        const decoder = new TextDecoder();
        const decoded = decoder.decode(data.slice(baseIndex, baseIndex + stringLength));
        // strip null padding
        return decoded.replace(/\0+$/g, "");
    }

    throw new Error(`Unsupported datatype: ${datatype} for column: ${column.field}`);
};

/**
 * Sets the cell value for the given datatype.
 *
 * @param column - Column object with metadata
 * @param dataIndex - The index of the value in data array
 * @param newValue - New value to replace the existing cell value
 * @throws Error if the operation fails
 */
export const setCellValueFromString = (
    column: LoadedDataColumn<DataType>,
    dataIndex: number,
    newValue: string,
): void => {
    if (!column.data) {
        throw new Error(`No data found in column: ${column.field}`);
    }

    const { data, values, stringLength, datatype } = column;

    // text / text16
    if (datatype === "text" || datatype === "text16") {
        if (!values) {
            throw new Error(`No values array found for column: ${column.field}`);
        }

        const maxValues = datatype === "text" ? 256 : 65536;
        const decodedValues = decodeCategoricalData(column);
        decodedValues[dataIndex] = newValue;
        rewriteCategoricalDataFromLabels(column, decodedValues, maxValues);
        return; 
    }

    // numeric
    if (datatype === "double" || datatype === "int32" || datatype === "integer") {
        const numValue = Number.parseFloat(newValue);
        if (!Number.isFinite(numValue)) {
            throw new Error(`Invalid number value "${newValue}" for ${datatype} column: ${column.field}`);
        }

        data[dataIndex] = numValue;
        return; 
    }

    // multitext
    if (datatype === "multitext") {
        if (!values) {
            throw new Error(`Missing values or stringLength for multitext column: ${column.field}`);
        }

        const capacity = getMultitextCapacity(column);
        const baseIndex = dataIndex * capacity;

        if (baseIndex + capacity > data.length) {
            throw new Error(`Index out of bounds for multitext column: "${column.field}"`);
        }

        const maxValues = 65536;
        const delimiter = getMultitextDelimiter(column);

        const parts = splitMultitextItems(newValue, delimiter);
        const decodedValues = decodeCategoricalData(column);
        for (let i = 0; i < capacity; i++) {
            decodedValues[baseIndex + i] = i < parts.length ? parts[i] : null;
        }

        rewriteCategoricalDataFromLabels(column, decodedValues, maxValues);

        return;
    }

    // unique
    if (datatype === "unique") {
        if (!stringLength) {
            throw new Error(`Missing stringLength for column: ${column.field}`);
        }

        const encoder = new TextEncoder();
        const encoded = encoder.encode(newValue);

        if (encoded.length > stringLength) {
            throw new Error(
                `Value "${newValue}" too long for column "${column.field}". Maximum length: ${stringLength}, actual length: ${encoded.length}`,
            );
        }

        const baseIndex = dataIndex * stringLength;
        if (baseIndex + stringLength > data.length) {
            throw new Error(`Index out of bounds when writing unique value to column: "${column.field}"`);
        }

        // Clear and set new value
        for (let i = 0; i < stringLength; i++) {
            data[baseIndex + i] = i < encoded.length ? encoded[i] : 0;
        }

        return;
    }

    throw new Error(`Unsupported datatype: ${datatype} for column: ${column.field}`);
};

/**
 * Replaces occurrences of findValue with replaceValue in the cell value (Excel-like behavior).
 * Instead of replacing the whole cell, it replaces the findValue substring within the cell's display value.
 *
 * @param searchColumn - Column name in which replace is performed
 * @param column - column object with metadata
 * @param findValue - Find value
 * @param replaceValue - Replace value
 * @param dataIndex - Index of the value in the data array
 * @throws Error if the operation fails
 * @returns True if a replacement was made, false if no match was found (but no error occurred)
 */
export const replaceMatches = (
    searchColumn: string,
    column: LoadedDataColumn<DataType> | undefined,
    findValue: string,
    replaceValue: string,
    dataIndex: number,
): boolean => {
    if (!column) {
        throw new Error(`No column found for replace operation in column: ${searchColumn}`);
    }

    if (!column.data) {
        throw new Error(`No data found in the column: ${searchColumn}`);
    }

    if (!findValue) {
        throw new Error(`Find value is empty for column: ${searchColumn}`);
    }

    let currentValue: string;
    try {
        currentValue = getCellValueAsString(column, dataIndex);
    } catch (error) {
        console.error(`Failed to get current cell value in column "${searchColumn}": ${error instanceof Error ? error.message : String(error)}`);
        throw new Error(
            `Failed to get current cell value in column "${searchColumn}": ${error instanceof Error ? error.message : String(error)}`,
        );
    }

    // No match found
    if (!currentValue.toLowerCase().includes(findValue.toLowerCase())) {
        return false;
    }

    //! case-insensitive
    const updatedCellValue = replaceValueInString(findValue, replaceValue, currentValue);

    if (updatedCellValue === currentValue) {
        return false;
    }

    // Set the new value
    try {
        setCellValueFromString(column, dataIndex, updatedCellValue);
        return true;
    } catch (error) {
        throw new Error(
            `Failed to set replaced value in column "${searchColumn}": ${error instanceof Error ? error.message : String(error)}`,
        );
    }
};
