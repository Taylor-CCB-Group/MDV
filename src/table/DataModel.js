import { getExportCsvStream } from "@/datastore/dataExportUtils";
import {
    getCellValueAsString,
    setCellValueFromString,
} from "@/react/utils/valueReplacementUtil";

const DEFAULT_MULTITEXT_CAPACITY = 24;
const DEFAULT_UNIQUE_STRING_LENGTH = 64;

/**
 * Helper to clone data type agnostically 
 */
function cloneTypedArrayBuffer(data) {
    const buffer = new SharedArrayBuffer(data.byteLength);
    const ctor = data.constructor;
    const copy = new ctor(buffer);
    copy.set(data, 0);
    return buffer;
}

/**
 * Helper to parse positive integer or use fallback
 */
function parsePositiveIntOrDefault(value, fallback, label) {
    if (value == null || value === "") {
        return fallback;
    }
    const parsed = Number.parseInt(value, 10);
    if (!Number.isInteger(parsed) || parsed <= 0) {
        throw new Error(`${label} must be a positive integer`);
    }
    return parsed;
}

class DataModel {
    /**
     * * @typedef {import("@/datastore/DataStore").default} DataStore
     * @param {DataStore} dataStore
     * @param {Object} config
     * @param {boolean?} config.autoupdate if true (or by default), listen for "filtered" in the dataStore and update the model
     */
    constructor(dataStore, config = {}) {
        /** @type {DataStore} */
        this.dataStore = dataStore;
        // why??? maybe it should make an initial call to updateModel?
        this.data = new Int32Array(0);
        const len = 0;
        for (let n = 0; n < len; n++) {
            this.data[n] = n;
        }
        this.size = len;
        this.listeners = {};
        if (config.autoupdate == null || config.autoupdate === true) {
            dataStore.addListener("sss", (type) => {
                if (type === "filtered") {
                    this.updateModel();
                }
            });
        }
    }

    getLength() {
        return this.size;
    }

    addListener(id, func) {
        this.listeners[id] = func;
    }

    setColumns(columns) {
        this.columns = columns;
    }

    sort(column, order) {
        if (column === "__index__") {
            this.data.sort();
            if (order === "desc") {
                this.data.reverse();
            }
            return;
        }
        const col = this.dataStore.getRawColumn(column);
        const c = this.dataStore.columnIndex[column];

        if (c.datatype === "unique") {
            const names = {};
            const tc = new TextDecoder();
            const sl = c.stringLength;
            const mu = order === "desc" ? 1 : -1;
            for (let i = 0; i < this.data.length; i++) {
                const index = this.data[i];
                names[index] = tc.decode(
                    col.slice(index * sl, index * sl + sl),
                );
            }
            this.data.sort((a, b) => {
                return names[a].localeCompare(names[b]) * mu;
            });
        } else {
            if (order === "desc") {
                this.data.sort((a, b) => {
                    let va = col[a];
                    let vb = col[b];
                    va = Number.isNaN(va) ? Number.MAX_VALUE : va;
                    vb = Number.isNaN(vb) ? Number.MAX_VALUE : vb;
                    return va - vb;
                });
            } else {
                this.data.sort((a, b) => {
                    let va = col[a];
                    let vb = col[b];
                    va = Number.isNaN(va) ? Number.MIN_VALUE : va;
                    vb = Number.isNaN(vb) ? Number.MIN_VALUE : vb;
                    return vb - va;
                });
            }
        }
    }

    /**
     * @param {string|Object} column - `string` when we clone a column 
     * and `object`: `{name, datatype, stringLength, delimiter}` when we create a new empty column
     * 
     * @param {string} [cloneCol] - optional legacy clone source
     */
    createColumn(column, cloneCol) {
        const columnSpec =
            typeof column === "string"
                ? {
                    name: column,
                    cloneColumn: cloneCol,
                    datatype: cloneCol ? undefined : "text",
                }
                : column;

        if (!columnSpec?.name) {
            throw new Error("Column name is required");
        }

        if (columnSpec.cloneColumn) {
            this._createClonedColumn(columnSpec.name, columnSpec.cloneColumn);
            return;
        }

        this._createEmptyColumn(columnSpec);
    }

    /**
     * 
     * @param {string} columnName - name of the new column
     * @param {string} cloneColumnName - name of the column to be cloned from
     * 
     * Fetch the column spec of the actual column and assign them to the new column specs
     * Create a deep copy of the data of the actual column to avoid mutations
     * Add the new column to the datastore
     */
    _createClonedColumn(columnName, cloneColumnName) {
        const clonedColumn = this.dataStore.columnIndex[cloneColumnName];
        if (!clonedColumn?.data) {
            throw new Error(`Column ${cloneColumnName} not found`);
        }

        const newColumn = {
            name: columnName,
            field: columnName,
            datatype: clonedColumn.datatype,
            editable: true,
        };

        // Assign the specs of the cloned column to new column
        if (clonedColumn.values) {
            newColumn.values = clonedColumn.values.slice();
        }
        if (clonedColumn.stringLength) {
            newColumn.stringLength = clonedColumn.stringLength;
        }
        if (clonedColumn.delimiter) {
            newColumn.delimiter = clonedColumn.delimiter;
        }
        if (clonedColumn.is_url) {
            newColumn.is_url = clonedColumn.is_url;
        }

        this.dataStore.addColumn(
            newColumn,
            // clone the data of existing column to avoid mutations
            cloneTypedArrayBuffer(clonedColumn.data),
            true,
        );
    }

    /**
     * 
     * @param {object} columnSpec - Column specs of new column to be created
     * 
     * Assign column specs based on the datatype and add the column to datastore
     */
    _createEmptyColumn(columnSpec) {
        const datatype = columnSpec.datatype || "text";
        const col = {
            name: columnSpec.name,
            field: columnSpec.name,
            datatype,
            editable: true,
        };

        let data;
        // Prefill the column specs with default values based on the datatype
        if (datatype === "text" || datatype === "text16") {
            col.values = [""];
            data = new Array(this.dataStore.size).fill("");
        } else if (
            datatype === "double" ||
            datatype === "integer" ||
            datatype === "int32"
        ) {
            //! int32 cannot truly represent NaN, so empty int32 cells currently coerce during storage and need a real missing-value strategy later.
            data = new Array(this.dataStore.size).fill(Number.NaN);
        } else if (datatype === "multitext") {
            const stringLength = parsePositiveIntOrDefault(
                columnSpec.stringLength,
                DEFAULT_MULTITEXT_CAPACITY,
                "multitext capacity",
            );
            col.values = [];
            col.stringLength = stringLength;
            col.delimiter = columnSpec.delimiter || ",";
            const buffer = new SharedArrayBuffer(this.dataStore.size * stringLength * 2);
            const arr = new Uint16Array(buffer);
            arr.fill(65535);
            data = buffer;
        } else if (datatype === "unique") {
            const stringLength = parsePositiveIntOrDefault(
                columnSpec.stringLength,
                DEFAULT_UNIQUE_STRING_LENGTH,
                "unique stringLength",
            );
            col.stringLength = stringLength;
            data = new SharedArrayBuffer(this.dataStore.size * stringLength);
        } else {
            throw new Error(`Unsupported datatype: ${datatype}`);
        }
        this.dataStore.addColumn(col, data, true);
    }

    removeColumn(col) {
        this.dataStore.removeColumn(col, true, true);
    }

    fillColumn(columnName, value, rowIndices, emptyOnly = false) {
        const col = this.dataStore.columnIndex[columnName];
        if (!col) {
            throw new Error(`Column ${columnName} not found`);
        }
        if (!col.editable) {
            throw new Error(`Column ${columnName} not editable`);
        }

        for (const rowIndex of rowIndices) {
            if (emptyOnly) {
                const currentValue = getCellValueAsString(col, rowIndex);
                if (currentValue !== "") {
                    continue;
                }
            }
            setCellValueFromString(col, rowIndex, value);
        }

        this.dataStore.dataChanged([columnName]);
    }

    async getDataAsBlob(delimiter = "\t", newline = "\n") {
        // return this.dataStore.getDataAsBlob(
        //     this.columns,
        //     this.data,
        //     delimiter,
        //     newline,
        // );
        const stream = await getExportCsvStream(this, delimiter, newline);
        const response = new Response(stream);
        return await response.blob();
    }

    getValueSuggestion(val, column) {
        const values = this.dataStore.getColumnValues(column);
        for (const v of values) {
            if (v === val) {
                return v;
            }
            if (v.startsWith(val)) {
                return v;
            }
        }
        return null;
    }

    _getValueIndex(value, col) {
        let index = col.values.indexOf(value);
        if (index === -1) {
            col.values.push(value);
            if (col.values.length > 256)
                throw new Error(
                    `text column '${col.name}' exceeded 256 values when adding '${value}'`,
                );
            index = col.values.length - 1;
        }
        return index;
    }

    updateValue(value, row, column, notify = true) {
        const col = this.dataStore.columnIndex[column];
        const valPos = this._getValueIndex(value, col);
        col.data[row] = valPos;
        if (notify) {
            this.dataStore.dataChanged([column]);
        }
    }

    updateRange(col, start, end, dir) {
        const arr = this.dataStore.getRawColumn(col);
        const val = dir === "+" ? arr[this.data[start]] : arr[this.data[end]];
        let i = dir === "+" ? start : end;
        const inc = dir === "+" ? +1 : -1;
        const lim = dir === "+" ? end + 1 : start - 1;
        for (; i !== lim; i += inc) {
            arr[this.data[i]] = val;
        }

        this.dataStore.dataChanged([col]);
    }

    replaceValues(value, replace, column, notify = true) {
        const col = this.dataStore.columnIndex[column];

        const valPos = this._getValueIndex(value, col);
        if (replace === "_all_") {
            for (let i = 0; i < this.data.length; i++) {
                col.data[this.data[i]] = valPos;
            }
        } else if (replace === "_blank_") {
            const bIndex = col.values.indexOf(value);
            //no blanks to replace
            if (bIndex === -1) {
                return;
            }
            for (let i = 0; i < this.data.length; i++) {
                const index = this.data[i];
                if (col.data[index] === bIndex) {
                    col.data[this.data[i]] = valPos;
                }
            }
        } else if (replace === "_delete_") {
            this.dataStore.cleanColumnData(column);
            return;
        } else {
            const repPos = this._getValueIndex(replace, col);
            for (let i = 0; i < this.data.length; i++) {
                const index = this.data[i];
                if (col.data[index] === valPos) {
                    col.data[this.data[i]] = repPos;
                }
            }
        }
        if (notify) {
            this.dataStore.dataChanged([column]);
        }
    }

    updateModel() {
        const arr = this.dataStore.filterArray;
        this.data = new Int32Array(this.dataStore.filterSize);
        const data = this.data;
        const len = this.dataStore.size;
        let index = 0;
        for (let n = 0; n < len; n++) {
            if (arr[n] === 0) {
                data[index++] = n;
            }
        }
        this.size = this.dataStore.filterSize;
        for (const id in this.listeners) {
            this.listeners[id]();
        }
    }

    getId(index) {
        return this.data[index];
    }

    getItemField(index, column) {
        return this.dataStore.getRowAsObject(this.data[index], [column])[
            column
        ];
    }

    getItem(index) {
        return this.dataStore.getRowAsObject(this.data[index], this.columns);
    }
}

export { DataModel };
