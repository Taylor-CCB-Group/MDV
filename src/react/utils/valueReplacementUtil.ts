import type { DataType, LoadedDataColumn } from "@/charts/charts";

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
 * Replaces the current value with the new replace value
 * 
 * Text datatype: Get or add replace value index and update the index in data array
 * Numeric datatype: Replace the value in the data array with the new value
 * Multitext datatype: Get or add new value index, find current value and replace the index
 * Unique datatype: Encode the replace value to bytes, update the existing value with new value (in bytes)
 * 
 * @param searchColumn - Column name of column in which replace is performed
 * @param column - Column Object with metadata
 * @param findValue - Find value
 * @param replaceValue - Replace value
 * @param dataIndex - Row index of data array in which replacement should be performed
 * @returns true if successful, false if not
 */
export const replaceMatches = (
    searchColumn: string,
    column: LoadedDataColumn<DataType> | undefined,
    findValue: string,
    replaceValue: string,
    dataIndex: number,
) => {
    if (!column) {
        console.error("No column found for replace value: ", replaceValue);
        return false;
    }

    if (!column.data) {
        console.error("No data found in the column: ", searchColumn);
        return false;
    }

    // Text datatype
    if (column.datatype === "text" || column.datatype === "text16") {
        if (!column.values) {
            console.error("No values found in the column: ", searchColumn);
            return false;
        }
        // Based on datasource.md
        const maxValues = column.datatype === "text" ? 256 : 65536;
        
        const valueIndex = getValueIndex(replaceValue, column.values, maxValues);
        column.data[dataIndex] = valueIndex;
        return true;
    }

    // Numeric datatype
    if (column.datatype === "double" || column.datatype === "int32" || column.datatype === "integer") {
        const replaceNumber = Number.parseFloat(replaceValue);

        if (!Number.isFinite(replaceNumber)) {
            console.error("Non-numeric value: ", searchColumn);
            return false;
        }

        //? Can we assign a decimal to the datatype integer or int32?
        //? What if a user inputs a decimal for those columns?

        column.data[dataIndex] = replaceNumber;
        return true;
    }

    // Multitext datatype
    if (column.datatype === "multitext") {
        const { values, data, stringLength } = column;
        if (!values || !stringLength) {
            console.error("No values found in the column: ", searchColumn);
            return false;
        }

        const baseIndex = dataIndex * stringLength;

        const findLower = findValue.toLowerCase();
        const replaceIndex = getValueIndex(replaceValue, column.values, 65536);

        let replaced = false;
        for (let rowIndex = 0; rowIndex < stringLength; rowIndex++) {
            const index = data[baseIndex + rowIndex];

            // No value check
            if (index === 65535) continue;
        
            // Out of bounds check
            if (index < 0 || index >= values.length) {
                console.error(`Index out of bounds for column: ${searchColumn}, skipping.`)
                continue;
            }

            const currentValue = values[index];

            //! Need to take a look at this
            // Only update the searched value
            if (currentValue.toLowerCase() === findLower) {
                data[baseIndex + rowIndex] = replaceIndex;
                replaced = true;
            }
        }

        return replaced;
    }

    // Unique datatype
    if (column.datatype === "unique") {
        const { data, stringLength } = column;

        if (!stringLength) {
            console.error("Missing string length for the column: ", searchColumn);
            return false;
        }

        const encoder = new TextEncoder();
        const replaceEncoded = encoder.encode(replaceValue);

        if (replaceEncoded.length > stringLength) {
            console.error("Value too long for the column: ", searchColumn);
            return false;
        }

        const baseIndex = dataIndex * stringLength;

        // Index out of bound check
        if (baseIndex + stringLength > data.length) {
            console.error("Data index out of bounds for the column: ", searchColumn);
            return false;
        }

        for (let i = 0; i < stringLength; i++) {
            // Update the existing value
            if (i < replaceEncoded.length) {
                data[baseIndex + i] = replaceEncoded[i];
            } else {
                // Add zeros at the end
                data[baseIndex + i] = 0;
            }
        }

        console.log("updated data: ", data);

        return true;
    }

    console.error("Replace not supported for datatype: ", column.datatype);
    return false;
};
