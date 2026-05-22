import type { DataType, LoadedDataColumn } from "@/charts/charts";
import {
    replaceMatches,
    replaceValueInString,
    setCellValueFromString,
} from "@/react/utils/valueReplacementUtil";
import { describe, expect, test } from "vitest";

describe("valueReplacementUtils", () => {
    // replaceValueInString
    describe("replaceValueInString", () => {
        test("empty args", () => {
            const cellValue = replaceValueInString("", "", "");
            expect(cellValue).toBe("");
        });

        test("valid args", () => {
            const updatedCellValue = replaceValueInString("p", "t", "Apple");
            expect(updatedCellValue).toBe("Attle");
        });

        test("empty replace value", () => {
            const updatedCellValue = replaceValueInString("A", "", "Apple");
            expect(updatedCellValue).toBe("pple");
        });
    });

    // setCellValueFromString
    describe("setCellValueFromString", () => {
        let column: LoadedDataColumn<DataType>;
        let dataIndex: number;
        let replaceValue: string;

        // text column
        describe("text column", () => {
            test("update value successfully", () => {
                column = {
                    field: "Col1",
                    data: [0],
                    values: ["Cat"],
                    datatype: "text",
                } as any;

                dataIndex = column.data.length - 1;
                replaceValue = "Ball";

                setCellValueFromString(column, dataIndex, replaceValue);
                expect(column.data[dataIndex]).toBe(0);
                expect(column.values).toEqual([replaceValue]);
                expect(column.values[column.data[dataIndex]]).toBe(replaceValue);
            });

            test("removes unused old labels and remaps text row data", () => {
                column = {
                    field: "Col1",
                    data: new Uint8Array([0, 1, 0]),
                    values: ["Value A", "Other"],
                    colors: ["#111111", "#222222"],
                    datatype: "text",
                } as any;

                setCellValueFromString(column, 0, "Value B");
                setCellValueFromString(column, 2, "Value B");

                expect(column.values).toEqual(["Other", "Value B"]);
                expect(Array.from(column.data)).toEqual([1, 0, 1]);
                expect(column.colors).toEqual(["#222222"]);
            });

            test("throw error for no values array", () => {
                column = {
                    field: "Col1",
                    data: ["A"],
                    values: undefined,
                    datatype: "text",
                } as any;

                expect(() => setCellValueFromString(column, 0, "")).toThrowError(
                    `No values array found for column: ${column.field}`,
                );
            });
        });

        // numeric column
        describe("numeric column", () => {
            test("update value successfully", () => {
                column = {
                    field: "Col1",
                    data: [1, 2, 3],
                    datatype: "double",
                } as any;

                dataIndex = 0;
                replaceValue = "4";

                setCellValueFromString(column, dataIndex, replaceValue);
                expect(column.data[dataIndex]).toBe(Number.parseFloat(replaceValue));
            });

            test("throw error for invalid number value", () => {
                column = {
                    field: "Col1",
                    data: [1, 2, 3],
                    datatype: "double",
                } as any;

                expect(() => setCellValueFromString(column, 0, "")).toThrowError(
                    `Invalid number value "" for ${column.datatype} column: ${column.field}`,
                );
            });
        });

        // multitext column
        describe("multitext column", () => {
            test("update value successfully", () => {
                column = {
                    field: "Col1",
                    data: [0, 1, 2],
                    values: ["val1", "val2", "val3"],
                    stringLength: 3,
                    datatype: "multitext",
                } as any;

                dataIndex = 0;
                replaceValue = "val4, val5, val6";

                setCellValueFromString(column, dataIndex, replaceValue);

                expect(Array.from(column.data)).toEqual([0, 1, 2]);
                expect(column.values).toEqual(["val4", "val5", "val6"]);
            });

            test("uses the column delimiter when parsing replacements", () => {
                column = {
                    field: "Col1",
                    data: [0, 1, 2],
                    values: ["val1", "val2", "val3"],
                    stringLength: 3,
                    delimiter: ";",
                    datatype: "multitext",
                } as any;

                dataIndex = 0;
                replaceValue = "val4; val5; val6";

                setCellValueFromString(column, dataIndex, replaceValue);

                expect(column.values).toEqual(["val4", "val5", "val6"]);
            });

            test("reserves the multitext empty-slot sentinel as unavailable for labels", () => {
                const values = Array.from({ length: 65535 }, (_, index) => `value-${index}`);
                const data = new Uint16Array(65536);
                for (let index = 0; index < values.length; index++) {
                    data[index] = index;
                }
                data[65535] = 65535;
                column = {
                    field: "Col1",
                    data,
                    values,
                    stringLength: 1,
                    datatype: "multitext",
                } as any;

                expect(() => setCellValueFromString(column, 65535, "new value")).toThrowError(
                    `Column exceeded 65535 values while setting: ${column.field}`,
                );
            });

            test("throw error for empty values array", () => {
                column = {
                    field: "Col1",
                    data: [1, 2, 3],
                    values: undefined,
                    datatype: "multitext",
                } as any;

                expect(() => setCellValueFromString(column, 0, "")).toThrowError(
                    `Missing values or stringLength for multitext column: ${column.field}`,
                );
            });
        });

        // unique column
        describe("unique column", () => {
            test("update value successfully", () => {
                column = {
                    field: "Col1",
                    data: [65, 66],
                    datatype: "unique",
                    stringLength: 1,
                } as any;

                dataIndex = 0;
                replaceValue = "C";
                const encodedLength = 67;

                setCellValueFromString(column, dataIndex, replaceValue);
                expect(column.data[dataIndex]).toBe(encodedLength);
            });

            test("throw error for invalid string length", () => {
                column = {
                    field: "Col1",
                    data: [1, 2, 3],
                    datatype: "unique",
                    stringLength: 0,
                } as any;

                expect(() => setCellValueFromString(column, 0, "")).toThrowError(
                    `Missing stringLength for column: ${column.field}`,
                );
            });

            test("throw error for maximum string length", () => {
                column = {
                    field: "Col1",
                    data: [65, 66],
                    datatype: "unique",
                    stringLength: 1,
                } as any;

                dataIndex = 0;
                replaceValue = "abcd";
                const encodedLength = 4;

                expect(() => setCellValueFromString(column, dataIndex, replaceValue)).toThrowError(
                    `Value "${replaceValue}" too long for column "${column.field}". Maximum length: ${column.stringLength}, actual length: ${encodedLength}`,
                );
            });

            test("throw error for index value out of bounds", () => {
                column = {
                    field: "Col1",
                    data: [65, 66],
                    datatype: "unique",
                    stringLength: 1,
                } as any;

                dataIndex = 2;
                replaceValue = "a";

                expect(() => setCellValueFromString(column, dataIndex, replaceValue)).toThrowError(
                    `Index out of bounds when writing unique value to column: "${column.field}"`,
                );
            });
        });

        test("throw error for no column data array", () => {
            column = {
                field: "Col1",
                data: undefined,
            } as any;

            expect(() => setCellValueFromString(column, 0, "")).toThrowError(
                `No data found in column: ${column.field}`,
            );
        });

        test("throw error for unsupported datatype of column", () => {
            column = {
                field: "Col1",
                datatype: "unsupported",
                data: [1],
            } as any;

            expect(() => setCellValueFromString(column, 0, "")).toThrowError(
                `Unsupported datatype: ${column.datatype} for column: ${column.field}`,
            );
        });
    });

    // replaceMatches
    describe("replaceMatches", () => {
        const searchColumn = "Col1";
        let column: LoadedDataColumn<DataType>;
        let findValue: string;
        let replaceValue: string;
        let dataIndex: number;

        test("replace matches successfully", () => {
            column = {
                field: searchColumn,
                data: [1, 2, 3],
                datatype: "double",
            } as any;
            findValue = "1";
            replaceValue = "4";
            dataIndex = 0;
            const replaced = replaceMatches(searchColumn, column, findValue, replaceValue, dataIndex);
            expect(replaced).toBeTruthy();
        });

        test("current value does not include find value", () => {
            column = {
                field: searchColumn,
                data: [1, 2, 3],
                datatype: "double",
            } as any;
            findValue = "0";
            replaceValue = "4";
            dataIndex = 0;
            const replaced = replaceMatches(searchColumn, column, findValue, replaceValue, dataIndex);
            expect(replaced).toBeFalsy();
        });

        test("current value equal to updated value", () => {
            column = {
                field: searchColumn,
                data: [1, 2, 3],
                datatype: "double",
            } as any;
            findValue = "1";
            replaceValue = "1";
            dataIndex = 0;
            const replaced = replaceMatches(searchColumn, column, findValue, replaceValue, dataIndex);
            expect(replaced).toBeFalsy();
        });

        test("throw error if column object is invalid", () => {
            expect(() => replaceMatches(searchColumn, undefined, "", "", 0)).toThrowError(
                `No column found for replace operation in column: ${searchColumn}`,
            );
        });

        test("throw error if data array is invalid", () => {
            column = {
                data: undefined,
            } as any;
            expect(() => replaceMatches(searchColumn, column, "", "", 0)).toThrowError(
                `No data found in the column: ${searchColumn}`,
            );
        });

        test("throw error if find value is invalid", () => {
            column = {
                data: [1],
            } as any;
            expect(() => replaceMatches(searchColumn, column, "", "", 0)).toThrowError(
                `Find value is empty for column: ${searchColumn}`,
            );
        });
    });
});
