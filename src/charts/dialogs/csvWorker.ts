import { parse } from 'papaparse';

self.onmessage = (event: MessageEvent) => {
  const file = event.data;

  parse(file, {
    header: true,
    dynamicTyping: true,
    worker: true,
    preview: 2, // Parse the first two rows
    complete: (results) => {
      const columnNames = results.meta.fields || [];
      const columnTypes: string[] = [];
      const secondRowValues = results.data[1] || {};
      const rowCount = results.data.length;
      const columnCount = columnNames.length;

      // Infer column types based on the first few rows
      columnNames.forEach((name) => {
        const columnValues = results.data.map((row) => row[name]);
        const uniqueTypes = new Set(columnValues.map((value) => typeof value));

        if (uniqueTypes.has('number')) {
          const isFloat = columnValues.some((value) => Number(value) !== Math.floor(Number(value)));
          columnTypes.push(isFloat ? 'float' : 'integer');
        } else if (uniqueTypes.has('boolean')) {
          columnTypes.push('boolean');
        } else if (uniqueTypes.has('string')) {
          const isDate = columnValues.every((value) => !isNaN(Date.parse(value)));
          columnTypes.push(isDate ? 'date' : 'string');
        } else {
          columnTypes.push('unknown');
        }
      });

      const secondRowValuesArray = columnNames.map((name) => secondRowValues[name]);

      self.postMessage({ columnNames, columnTypes, secondRowValues: secondRowValuesArray, rowCount, columnCount });
    },
    error: (error) => {
      self.postMessage({ error: error.message });
    }
  });
};
