import { parse, ParseResult } from 'papaparse';

interface PreviewData {
  columnNames: string[];
  columnTypes: string[];
  secondRowValues: (string | number | boolean | null)[];
  previewRowCount: number;
  columnCount: number;
  rowCount?: number;
}

self.onmessage = (event: MessageEvent) => {
  const file: File = event.data;

  // Initialize variables for row counting
  let totalRowCount = 0;
  let previewData: PreviewData = {
    columnNames: [],
    columnTypes: [],
    secondRowValues: [],
    previewRowCount: 0,
    columnCount: 0
  };

  // Function to count the total number of rows
  const countRows = (): Promise<number> => {
    return new Promise((resolve, reject) => {
      parse(file, {
        worker: true,
        step: () => {
          totalRowCount++;
        },
        complete: () => {
          resolve(totalRowCount);
        },
        error: (error) => {
          reject(error);
        }
      });
    });
  };

  // Function to get the preview rows and column information
  const getPreview = (): Promise<void> => {
    return new Promise((resolve, reject) => {
      parse(file, {
        header: true,
        dynamicTyping: true,
        worker: true,
        preview: 2, // Parse the first two rows
        complete: (results: ParseResult<any>) => {
          const columnNames = results.meta.fields || [];
          const columnTypes: string[] = [];
          const secondRowValues = results.data[1] || {};
          const previewRowCount = results.data.length;
          const columnCount = columnNames.length;

          // Infer column types based on the first few rows
          columnNames.forEach((name) => {
            const columnValues = results.data.map((row: any) => row[name]);
            const uniqueTypes = new Set(columnValues.map((value: any) => typeof value));

            if (uniqueTypes.has('number')) {
              const isFloat = columnValues.some((value: any) => Number(value) !== Math.floor(Number(value)));
              columnTypes.push(isFloat ? 'float' : 'integer');
            } else if (uniqueTypes.has('boolean')) {
              columnTypes.push('boolean');
            } else if (uniqueTypes.has('string')) {
              const isDate = columnValues.every((value: any) => !isNaN(Date.parse(value)));
              columnTypes.push(isDate ? 'date' : 'string');
            } else {
              columnTypes.push('unknown');
            }
          });

          const secondRowValuesArray = columnNames.map((name) => secondRowValues[name]);

          previewData = { columnNames, columnTypes, secondRowValues: secondRowValuesArray, previewRowCount, columnCount };

          resolve();
        },
        error: (error) => {
          reject(error);
        }
      });
    });
  };

  // Perform row counting and get the preview in sequence
  getPreview()
    .then(() => countRows())
    .then((totalRowCount) => {
      // Update the rowCount with totalRowCount
      previewData.rowCount = totalRowCount;
      self.postMessage(previewData);
    })
    .catch((error) => {
      self.postMessage({ error: error.message });
    });
};
