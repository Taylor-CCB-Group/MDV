import type { DataType, LoadedDataColumn } from "@/charts/charts";

//? What if find and cell value are empty strings "" and we want to update it
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

/**
 * Gets or adds the index of the new value in values array
 * If the new value doesn't exist in the array, it adds and returns the new index
 *
 * @param replaceValue - Replace value to find or add in values array
 * @param values - Array of unique strings
 * @param maxValues - Maximum number of values allowed (256 for text, 65536 for text16 or multitext)
 * @returns index of the new value in values array
 * @throws Error if length of values array exceeds the maxValues
 *
 */
const getValueIndex = (replaceValue: string, values: string[], maxValues: number) => {
    let valueIndex = values.indexOf(replaceValue);

    if (valueIndex === -1) {
        // Check if length exceeds the max values of the array
        if (values.length >= maxValues) {
            throw new Error(`Column exceeded ${maxValues} values while adding: ${replaceValue}`);
        }
        // Create a new value if it doesn't exist in values array
        values.push(replaceValue);
        valueIndex = values.length - 1;
    }

    return valueIndex;
};

/**
 * Gets the current cell value for the given datatype.
 *
 * @param column - Column object with metadata
 * @param dataIndex - Index of the value in the data array
 * @returns Cell value of the column in the table
 * @throws Error if column data is invalid or missing
 */
const getCurrentCellValue = (column: LoadedDataColumn<DataType>, dataIndex: number): string => {
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
        if (!values || !stringLength) {
            throw new Error(`Missing values or stringLength for multitext column: ${column.field}`);
        }

        const baseIndex = dataIndex * stringLength;
        const cellData = data.slice(baseIndex, baseIndex + stringLength);

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
            .join(", ");
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
        const valueIndex = getValueIndex(newValue, values, maxValues);
        data[dataIndex] = valueIndex;
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

    //! review this
    // multitext
    if (datatype === "multitext") {
        if (!values || !stringLength) {
            throw new Error(`Missing values or stringLength for multitext column: ${column.field}`);
        }

        const baseIndex = dataIndex * stringLength;

        const maxValues = 65536;

        // Split by comma and trim each value
        const parts = newValue
            .split(",")
            .map((p) => p.trim())
            .filter((p) => p.length > 0);

        // Clear existing values
        for (let i = 0; i < stringLength; i++) {
            data[baseIndex + i] = 65535; // Empty marker
        }

        // Set new values
        for (let i = 0; i < Math.min(parts.length, stringLength); i++) {
            const valueIndex = getValueIndex(parts[i], values, maxValues);
            data[baseIndex + i] = valueIndex;
        }

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
        currentValue = getCurrentCellValue(column, dataIndex);
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
