import type { DataType, LoadedDataColumn } from "@/charts/charts";

const getValueIndex = (replaceValue: string, values: string[], maxValues: number) => {
    let valueIndex = values.indexOf(replaceValue);

    if (valueIndex === -1) {
        if (values.length >= maxValues) {
            throw new Error(`Column exceeded ${maxValues} values while adding: ${replaceValue}`);
        }
        values.push(replaceValue);
        valueIndex = values.length - 1;
    }

    return valueIndex;
};

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

    if (column.datatype === "double" || column.datatype === "int32" || column.datatype === "integer") {
        const replaceNumber = Number.parseFloat(replaceValue);
        console.log("replace value: ", replaceValue, typeof replaceValue);
        console.log("replace number: ", replaceNumber, typeof replaceNumber);
        if (!Number.isFinite(replaceNumber)) {
            console.error("Non-numeric value: ", searchColumn);
            return false;
        }

        //? Can we assign a decimal to the datatype integer or int32?
        //? What if a user inputs a decimal for those columns?

        column.data[dataIndex] = replaceNumber;
        return true;
    }

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

            if (index === 65535) continue;

            const currentValue = values[index];
            if (currentValue.toLowerCase() === findLower) {
                data[baseIndex + rowIndex] = replaceIndex;
                replaced = true;
            }
        }

        return replaced;
    }

    if (column.datatype === "unique") {
        const { data, stringLength } = column;

        if (!stringLength) {
            console.error("Missing string length for the column: ", searchColumn);
            return false;
        }

        const encoder = new TextEncoder();
        const replaceEncoded = encoder.encode(replaceValue);
        console.log("encoded replace value: ", replaceEncoded);
        if (replaceEncoded.length > stringLength) {
            console.error("Value too long for the column: ", searchColumn);
            return false;
        }

        const baseIndex = dataIndex * stringLength;

        if (baseIndex + stringLength > data.length) {
            console.error("Data index out of bounds for the column: ", searchColumn);
            return false;
        }

        console.log("baseIndex: ", baseIndex);
        console.log("string length: ", stringLength);

        for (let i = 0; i < stringLength; i++) {
            if (i < replaceEncoded.length) {
                console.log("in if value: ", data[baseIndex + i], replaceEncoded[i]);
                data[baseIndex + i] = replaceEncoded[i];
            } else {
                console.log("in else value: ", data[baseIndex + i]);
                data[baseIndex + i] = 0;
            }
        }

        console.log("updated data: ", data);

        return true;
    }

    console.error("Replace not supported for datatype: ", column.datatype);
    return false;
};
