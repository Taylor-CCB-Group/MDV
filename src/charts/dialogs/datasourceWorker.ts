import { parse, type ParseResult } from "papaparse";

type ColumnType = "float" | "integer" | "boolean" | "date" | "string" | "unknown";
type DataValue = string | number | boolean | null;
type FileType = 'csv' | 'tsv';

interface FileMessage {
    file: File;
    fileType?: FileType;
}

interface PreviewData {
    columnNames: string[];
    columnTypes: ColumnType[];
    secondRowValues: DataValue[];
    previewRowCount: number;
    columnCount: number;
    rowCount?: number;
}

interface ErrorResponse {
    error: string;
}

type WorkerResponse = PreviewData | ErrorResponse;
type ParsedRow = Record<string | number, DataValue>;

self.onmessage = (event: MessageEvent<File | FileMessage>) => {
    const file: File = 'file' in event.data ? event.data.file : event.data;
    const fileType: FileType = ('fileType' in event.data ? event.data.fileType : 'csv') as FileType;
    const delimiter: string = fileType === 'tsv' ? '\t' : ',';

    let totalRowCount = 0;
    let previewData: PreviewData = {
        columnNames: [],
        columnTypes: [],
        secondRowValues: [],
        previewRowCount: 0,
        columnCount: 0,
    };

    const countRows = (): Promise<number> => {
        return new Promise((resolve, reject) => {
            parse(file, {
                worker: true,
                delimiter,
                step: () => {
                    totalRowCount++;
                },
                complete: () => resolve(totalRowCount),
                error: (error: Error) => reject(error),
            });
        });
    };

    const getPreview = (): Promise<void> => {
        return new Promise((resolve, reject) => {
            parse(file, {
                header: true,
                delimiter,
                dynamicTyping: true,
                worker: true,
                preview: 2,
                complete: (results: ParseResult<ParsedRow>) => {
                    const columnNames = results.meta.fields || [];
                    const columnTypes: ColumnType[] = [];
                    const secondRowValues = results.data[1] || {};
                    const previewRowCount = results.data.length;
                    const columnCount = columnNames.length;

                    columnNames.forEach((name: string | number) => {
                        const columnValues = results.data.map(row => row[name]);
                        const uniqueTypes = new Set(
                            columnValues.map(value => typeof value)
                        );

                        if (uniqueTypes.has("number")) {
                            const isFloat = columnValues.some(
                                value => typeof value === 'number' && 
                                value !== Math.floor(value)
                            );
                            columnTypes.push(isFloat ? "float" : "integer");
                        } else if (uniqueTypes.has("boolean")) {
                            columnTypes.push("boolean");
                        } else if (uniqueTypes.has("string")) {
                            const isDate = columnValues.every(
                                value => typeof value === 'string' && 
                                !Number.isNaN(Date.parse(value))
                            );
                            columnTypes.push(isDate ? "date" : "string");
                        } else {
                            columnTypes.push("unknown");
                        }
                    });

                    const secondRowValuesArray = columnNames.map(
                        name => secondRowValues[name]
                    );

                    previewData = {
                        columnNames,
                        columnTypes,
                        secondRowValues: secondRowValuesArray,
                        previewRowCount,
                        columnCount,
                    };

                    resolve();
                },
                error: (error: Error) => reject(error),
            });
        });
    };

    getPreview()
        .then(() => countRows())
        .then((totalRowCount) => {
            previewData.rowCount = totalRowCount;
            self.postMessage(previewData as WorkerResponse);
        })
        .catch((error: Error) => {
            self.postMessage({ error: error.message } as WorkerResponse);
        });
};