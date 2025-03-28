import { getExportCsvStream } from "@/datastore/dataExportUtils";

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
     * @param {string} name
     * @param {string} cloneCol
     */
    createColumn(name, cloneCol) {
        const col = {
            name: name,
            field: name,
            datatype: "text",
            values: [""],
        };
        const buff = new SharedArrayBuffer(this.dataStore.size);
        if (cloneCol) {
            const arr = this.dataStore.getRawColumn(cloneCol);
            const newArr = new Uint8Array(buff);
            newArr.set(arr);
            col.values = this.dataStore.getColumnValues(cloneCol).slice(0);
        }
        this.dataStore.addColumn(col, buff, true);
    }

    removeColumn(col) {
        this.dataStore.removeColumn(col, true, true);
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
