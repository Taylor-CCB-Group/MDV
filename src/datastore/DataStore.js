import Dimension from "./Dimension.js";
//register dimensions
import "./CategoryDimension.js";
import "./RangeDimension.js";
import "./CatColDimension.js";
import "./CatRangeDimension.js";
import "./DensityDimension.js";
import "./ValueSetDimension";
import "./SortableDimension.js";
import "./DeepToolsDimension.js";
import { scaleLinear, scaleSymlog } from "d3-scale";
import { getColorLegend, getColorBar } from "../utilities/Color.js";
import { quantileSorted } from "d3-array";
import { makeObservable, observable, action } from "mobx";
import { isColumnNumeric, isColumnText } from "../utilities/Utilities";
import { isDatatypeNumeric } from "@/lib/utils";


class DataStore {
    /**
     * Creates a data store object based on the given configuration. Does not
     * initially load any data, but will hold references to column data once loaded.
     *
     * [Data Source](../../docs/extradocs/datasource.md)
     * @typedef {import("@/charts/charts").DataSource} DataSource
     * @param {number} size - the number of rows(items) of the data structure
     * @param {DataSource} [config] - setup information for the datastore.
     * @param {Object[]} [config.columns] - an array of column objects, specifying
     * the metadata for data structure, see {@link DataStore#addColumn}
     * @param {Object[]} [config.columnGroups] - an array of objects each has
     * name, and a list of columns in that group
     * @param {Object} [config.links] - an object describing how this DataStore
     * links with other DataStore
     * @param {Object} [config.images] - an object describing thumbnails which
     * are associated with each item/row in the DataStore
     * @param {Object} [config.large_images] - an object describing large images which
     * are associated with each item/row in the DataStore
     * @param {Object} [config.offsets] - an object specifying which columns can
     * have values that can be transformed and rotated and any transformations/
     * rotations applied to them
     */
    constructor(size, config = {}, dataLoader = null) {
        // by keeping a reference to the original config metadata, we can later compare
        // it to new config from server to determine if it has changed etc.
        this.config = config;
        /** @type {number} */
        this.size = size;
        /** @type {number} */
        this.filterSize = size;
        /**
         * @typedef {import("@/charts/charts.js").DataType} DataType
         * @type {import("@/charts/charts.js").DataColumn<DataType>[]} */
        this.columns = [];
        /**
         * @type {{[k: string]: import("@/charts/charts.js").DataColumn<DataType> | undefined}} */
        this.columnIndex = {};
        this.listeners = {};
        this.indexes = {};
        this.filterBuffer = new SharedArrayBuffer(size);
        this.filterArray = new Uint8Array(this.filterBuffer);
        //keep track of dimensions and the columns they represent
        this.dimensions = [];
        this.name = config.name;
        this.avivator = config.avivator;

        this.textDecoder = new TextDecoder();
        this.columnGroups = {};
        this.subtypeToGroup = {};
        this.interactions = config.interactions;
        this.deeptools = config.deeptools;

        //info about subgroups and their datasources
        this.subgroups = {};
        this.subgroupDataSources = new Set();
        this.accessOtherDataStore = [];
        this.syncColumnColors = [];
        this.linkColumns = [];
        this.regions = config.regions;

        //for react / mobx... makeAutoObservable causes problems with webworker.postMessage:
        // `DOMException: Failed to execute 'postMessage' on 'Worker': [object Array] could not be cloned.`
        // because it adds a mobx proxy to the array, which can't be cloned...
        // so for now, rather than using makeAutoObservable, we'll just manually add the mobx stuff we need:
        makeObservable(this, {
            _filteredIndicesPromise: observable,
            _callListeners: action,
            //columns: observable, //will probably have clone issues...
        });

        // for re-usable filteredIndices
        this.addListener(
            "invalidateFilteredIndicesCache",
            action((type, _data) => {
                //if (this._filteredIndicesPromise) this._filteredIndicesPromise.cancel(); // relevant? any test-cases to consider?
                // Make the promise null only if the data is filtered, changed or added
                if (type === "filtered" || type === "data_changed" || type === "data_added") {
                    this._filteredIndicesPromise = null;
                }
            }),
        );

        if (config.row_data_loader) {
            if (!dataLoader?.rowDataLoader) {
                console.warn(
                    `datasource ${this.name} requires a row data loader but none has been supplied`,
                );
            } else {
                this.rowDataLoader = dataLoader.rowDataLoader;
                this.rowData = new Map();
            }
        }
        if (config.binary_data_loader) {
            if (!dataLoader?.binaryDataLoader) {
                console.warn(
                    `datasource ${this.name} requires a binary data loader but none has been supplied`,
                );
            } else {
                this.binaryDataLoader = dataLoader.binaryDataLoader;
                this.binaryData = new Map();
            }
        }

        /** @type {string[]} column names */
        this.columnsWithData = [];
        this.dirtyColumns = {
            added: {},
            removed: {},
            data_changed: {},
            colors_changed: {},
        };
        this.dirtyMetadata = new Set();

        if (config.offsets) {
            config.offsets.values = config.offsets.values || {};
            this.offsets = config.offsets;
        }
        this.images = config.images;
        this.genome_browser = config.genome_browser;
        this.tree_diagram = config.tree_diagram;

        if (config.columns) {
            for (const c of config.columns) {
                this.addColumn(c, c.data);
            }
        }

        if (config.columnGroups) {
            for (const item of config.columnGroups) {
                this.addColumnGroup(item);
            }
        }
        if (config.large_images) {
            this.large_images = config.large_images;
        }
        if (config.links) {
            this.links = config.links;
            for (const ds in config.links) {
                const link = config.links[ds];
                if (link.rows_as_columns) {
                    const sg = link.rows_as_columns.subgroups;
                    this.subgroupDataSources.add(ds);
                    for (const t in sg) {
                        this.subgroups[t] = {
                            subgroup: sg[t],
                            dataSource: ds,
                            name_column: link.rows_as_columns.name_column,
                        };
                    }
                }
                if (link.access_data) {
                    this.accessOtherDataStore.push({
                        dataSource: ds,
                        index: link.index,
                    });
                }
                if (link.sync_column_colors) {
                    this.syncColumnColors.push({
                        dataSource: ds,
                        columns: link.sync_column_colors,
                    });
                }
                if (link.columns) {
                    this.linkColumns.push({
                        dataSource: ds,
                        columns: link.columns,
                    });
                }
            }
        }
    }

    /**
     *
     * @returns {boolean} - true if the DataStore has been filtered
     *
     */
    isFiltered() {
        return this.filterSize !== this.size;
    }

    /**
     * Adds a listener to the datastore that will be called when an event occurs,
     * passing the event type and any data. There are the following different types
     * of event:-
     * <ul>
     *   <li> filtered - called when a filter is applied. The dimension doing the
     *   the filtering is passed as data </li>
     *   <li> data_highlighted - called when certain indexes have been highlighted </li>
     *   <li> column_removed - called just before a column is removed </li>
     *   <li> data_changed- called when data has changed giving a list of columns where
     *   the data has changed </li>
     * </ul>
     * @param {string} id - a unique id identifying the listener
     * @param {(type: string, data: any) => void} listener - a function that accepts two paramaters, the type
     * of event and the data associated with it.
     */
    addListener(id, listener) {
        //PJT: XXX: this replaces any existing listener with the same id, probably not intended
        this.listeners[id] = listener;
    }

    /**
     * Removes the specified listener from the datastore
     * @param {string} id The id of the listener to remove
     */
    removeListener(id) {
        delete this.listeners[id];
    }

    _callListeners(type, data) {
        for (const id in this.listeners) {
            try {
                this.listeners[id](type, data);
            } catch (e) {
                console.error(
                    `Error calling listener '${id}' for event '${type}'`,
                    e,
                );
            }
        }
    }

    /**
     * Removes all filters from the datastore,
     * More efficient than removing each filter individually
     * Filters on dimnensions with a noclear property will not
     * be removed
     */
    removeAllFilters() {
        this.filterArray.fill(0);
        const noclear = [];
        for (const dim of this.dimensions) {
            if (dim.noClear) {
                noclear.push(dim);
                continue;
            }
            dim.filterArray.fill(0);
            dim.filterMethod = null;
            if (dim.bgfArray) {
                for (let i = 0; i < this.size; i++) {
                    if (dim.bgfArray[i] === 0) {
                        dim.filterArray[i] = 2;
                    }
                }
            }
        }
        this.filterSize = this.size;
        //need to re-add the noclear filters (if any)
        for (const dim of noclear) {
            const f = dim.filterArray;
            for (let n = 0; n < f.length; n++) {
                if (f[n] === 1) {
                    if (++this.filterArray[n] === 1) {
                        this.filterSize--;
                    }
                }
            }
        }
        this._callListeners("filtered", "all_removed");
    }

    /**
     * This method should be called if any data has been modified, specifying the
     * columns involved.
     * Any dimensions will re-filter if necessary i.e. if the modified columns are
     * involved in the filter and single filtered event will be broadcast.
     * All 'data_changed' listeners will be informed with the columns changed and whether
     * filtering has already occurred (in which case updating may have already occurred
     * if the object listening to data changes also listens to filtered events)
     * @param {string[]} columns - a list of column/fields whose data has been modified
     * @param {boolean} [is_dirty=true] - if false (default is true) then the column
     * will not be tagged as dirty. This my be the case if the change was pulled from
     * the backend for example
     */
    dataChanged(columns, is_dirty = true) {
        let hasFiltered = false;
        for (const d of this.dimensions) {
            //this method will not call any listeners- wait until all done
            if (d.reFilterOnDataChanged(columns)) {
                hasFiltered = true;
            }
        }
        if (is_dirty) {
            for (const c of columns) {
                this.setColumnIsDirty(c);
            }
        }
        //at least one of the dimensions has had to refilter?
        if (hasFiltered) {
            this._callListeners("filtered");
        }
        this._callListeners("data_changed", {
            columns: columns,
            hasFiltered: hasFiltered,
        });
    }

    /**
     * This method calls any listeners to 'highlight' any rows specified e.g.
     * rows in a table or points in a scatter plot
     * @param {number[]} indexes an array of indexes to items that should be highlighted
     * @param {object} source the object doing the highlighting
     */
    dataHighlighted(indexes, source) {
        indexes.sort((a, b) => a - b);
        this.highightedData = indexes;
        this.loadRowData(indexes[0]).then((data) =>
            this._callListeners("data_highlighted", { indexes, source, data }),
        );
    }

    async loadRowData(index) {
        //no data loader
        if (!this.rowDataLoader) {
            return null;
        }
        // is it in cache (stores null values as well)
        if (this.rowData.has(index)) {
            return this.rowData.get(index);
        }
        const rd = await this.rowDataLoader(this.name, index);
        this.rowData.set(index, rd);
        return rd;
    }

    async loadBinaryData(name) {
        //no data loader
        if (!this.binaryDataLoader) {
            return null;
        }
        // is it in cache (stores null values as well)
        if (this.binaryData.has(name)) {
            return this.binaryData.get(name);
        }
        const rd = await this.binaryDataLoader(this.name, name);
        this.binaryData.set(name, rd);
        return rd;
    }

    /**
     * @returns The indexes of items that have been highligted
     */
    getHighlightedData() {
        return this.highightedData;
    }

    /**
     * Broadcast a filter event to all listeners
     */
    triggerFilter() {
        this._callListeners("filtered");
    }

    /**
     * Tag that the colunm's data has changed and is not synched with the backend
     * @param {string} col - the column that has been modified
     */
    setColumnIsDirty(col) {
        //new column anyway so column already 'dirty'
        if (this.dirtyColumns.added[col]) {
            return;
        }
        this.dirtyColumns["data_changed"][col] = true;
    }

    /**
     * Specify that all data (and metadata) has beem synched with
     * the backend
     */
    setAllColumnsClean() {
        this.dirtyColumns.added = {};
        this.dirtyColumns.removed = {};
        this.dirtyColumns.data_changed = {};
        this.dirtyColumns.colors_changed = {};
        this.dirtyMetadata.clear();
    }

    /**
     * Returns the current filter, which is just an array corresponding
     * to the index of the row, which contains 0 if it is present or
     * greater than 0 if it has been filtrered out. Do not modify the array.
     * To check if an row is filtered use {@link DataStore#isRowFiltered}
     * @returns {Uint8Array}
     */
    getFilter() {
        return this.filterArray;
    }

    /**
     * Returns true if the row is in the filter and false if it
     * has been filtered out
     * @param {number} index The index of the row
     * @returns {boolean} whether the row has been filtered out
     */
    isRowFiltered(index) {
        return this.filterArray[index] === 0;
    }

    //delete
    _calculateCategories(column) {
        const vs = column.values;
        const d = column.data;
        const ci = {};
        for (const n of vs) {
            ci[n] = 0;
        }
        for (let i = 0; i < d.length; i++) {
            const v = vs[d[i]];
            ci[v]++;
        }
    }

    /*
     * todo infer generic type for column, not sure how to do this in jsdoc
     * @typedef {import("@/charts/charts.js").DataColumn<any>} Column
     *  */
    /**
     * Adds a column's metadata and optionally it's data to the DataStore
     * [Data Source](../../docs/extradocs/datasource.md)
     * @param {Omit<Column, 'getValue'>} column An object describing the column
     * @param {Array|SharedArrayBuffer} [data] The data for the column
     * @param {boolean} [dirty=false] true for columns that are not synched with the backend,
     * ~~for example 'virtual' columns added from a `rows_as_columns` link~~
     */
    addColumn(column, data = null, dirty = false) {
        /** @type {Column} */
        const c = {
            ...column, //may be useful to get any other properties we didn't explicitly copy before.
            getValue: (index) => {
                const col = this.columnIndex[column.field];
                if (!col.data) {
                    console.error(`Column ${column.field} has no data`);
                    return;
                }
                let v = col.data[index];
                if (col.datatype === "text" || col.datatype === "text16") {
                    v = col.values[v];
                } else if (
                    col.datatype === "double" ||
                    col.datatype === "integer" ||
                    col.datatype === "int32"
                ) {
                    if (Number.isNaN(v)) {
                        v = "missing";
                    }
                }
                //multitext displayed as comma delimited values
                else if (col.datatype === "multitext") {
                    const delim = ", ";
                    const d = col.data.slice(
                        index * col.stringLength,
                        index * col.stringLength + col.stringLength,
                    );
                    v = Array.from(d.filter((x) => x !== 65535))
                        .map((x) => col.values[x])
                        .join(delim);
                } else {
                    v = this.textDecoder.decode(
                        col.data.slice(
                            index * col.stringLength,
                            index * col.stringLength + col.stringLength,
                        ),
                    );
                    v = v.replaceAll("\0", "");
                }
                return v;
            },
            subgroup: column.subgroup,
        };
        if (!c.field) {
            c.field = column.name;
        }
        if (isColumnText(column) && !c.values) {
            c.values = [""];
        }
        this.columns.push(c);
        this.columnIndex[c.field] = c;
        if (data) {
            this.setColumnData(column.field, data);
        }
        if (dirty) {
            this.dirtyColumns.added[column.field] = true;
        }
        // can we infer/check the type of column here?
        return c;
    }

    /**
     * This method will return (case-insensitive) any values in the column
     * which contain the specified text (unique/text/multitext columns only)
     * @param {string} text - the query value
     * @param {string} column - the column to query
     * @param {number} number  - the maximum number of results to return
     * @returns {{value: string, index: number}[]} An array of objects with
     *  - value - the actual text match
     *  - index - for unique columns, the row index and for
     *   text/multitext its index in the column's values array
     */
    getSuggestedValues(text, column, number = 10) {
        const col = this.columnIndex[column];
        const tupper = text.toLowerCase();
        const tlower = text.toUpperCase();
        const matches = [];
        if (col.datatype === "unique") {
            const e = new TextEncoder();
            const bupper = e.encode(tupper);
            const blower = e.encode(tlower);
            const len = text.length;
            const d = col.data;
            const sl = col.stringLength;
            for (let i = 0; i < d.length - len; i++) {
                let match = true;
                for (let a = 0; a < len; a++) {
                    if (bupper[a] !== d[i + a] && blower[a] !== d[i + a]) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    const index = Math.floor(i / sl);
                    const v = this.textDecoder
                        .decode(col.data.slice(index * sl, index * sl + sl))
                        .replaceAll("\0", "");
                    matches.push({
                        value: v,
                        index: index,
                    });
                    if (matches.length > number) {
                        break;
                    }
                }
            }
        } else if (col.datatype === "text" || col.datatype === "multitext") {
            const tlength = text.length;
            for (let i = 0; i < col.values.length; i++) {
                const v = col.values[i];
                let match = true;
                for (let n = 0; n < v.length - tlength; n++) {
                    match = true;
                    for (let a = 0; a < tlength; a++) {
                        if (v[a + n] !== tupper[a] && v[a + n] !== tlower[a]) {
                            match = false;
                            break;
                        }

                        console.log("match");
                    }
                    if (match) {
                        break;
                    }

                    console.log(n);
                }
                if (match) {
                    matches.push({
                        value: v,
                        index: i,
                    });
                    if (matches.length > number) {
                        break;
                    }
                }
            }
        }
        return matches;
    }

    /**
     * For columns where the metadata is not listed in the `columns` of the DataSource/DataStore[*],
     * in particular in the case of columns that result from a 'rows_as_columns' link,
     * and perhaps in future for other features involving 'virtual' columns, such as generating synthetic data
     * or computing new statistics at runtime.
     *
     * This method will use a formatted "field" string to create a 'virtual' column.
     * We may revise this specification to have a more structured object representation
     * in place of the string format that will need to be parsed. The current system will
     * lead to edge-cases for example if column names contain the '|' character that may lead to unpleasant
     * user-experience at a later date. It also makes it less easy to reason about the code.
     *
     * Note that while this method will add relevant metadata to the DataStore, it will not fetch (or otherwise compute)
     * the data.
     *
     * [*] remarkable that I find this distinction confusing, this is perhaps an aspect of the design that could be improved
     *
     * @param {string} field - a formatted string that contains the metadata for the column.
     * The format is `subtype|name|index` where (for now):
     * - `subtype` a name of a "subgroup" (whatever that means), like "Gene scores"
     * - `name` the name to be given to the new virtual column, like "GENE (Gene scores)"
     * - `index` the index of the column in the source data
     */
    addColumnFromField(field) {
        const data = field.split("|");
        if (data.length !== 3) {
            throw new Error(
                `field ${field} cannot be found in datasource ${this.name}`,
            );
        }
        let g = this.subtypeToGroup[data[0]];
        let sg = null;
        if (!g) {
            g = this.subgroups[data[0]];
            sg = g.subgroup;
        } else {
            sg = g.subgroups[data[0]];
        }
        if (!this.columnIndex[field]) {
            return this.addColumn({
                name: data[1],
                field: field,
                datatype: "double",
                subgroup: sg.name,
                sgindex: data[2],
                sgtype: sg.type,
            });
        }
        return this.columnIndex[field];
    }

    /**
     * returns a list of column name and fields (which have data) sorted by name
     * @param {string[]|string} [filter] - can be either be a string -'number', 'string', 'all' or a column type.
     * Or an array of column types
     * @param {boolean} [addNone=false] if true then an extra object will be added to the list
     * with name 'None' and field '__none__'
     * @returns {Partial<Column>[]}  An array of objects containing name,field and datatype
     * properties. The columns are ordered by main columns followed by subgroups.
     * Each group is ordered alphabetically
     */
    getColumnList(filter = null, addNone = false) {
        const columns = [];
        const sgDataSources = {};
        const sgs = Object.keys(this.subgroups).map((x) => {
            sgDataSources[x] = [];
            return x;
        });
        const f_array = Array.isArray(filter);
        const has_sgs = sgs.length !== 0;
        for (const f in this.columnIndex) {
            const c = this.columnIndex[f];
            if (filter) {
                if (f_array) {
                    if (filter.indexOf(c.datatype) === -1) {
                        continue;
                    }
                } else if (filter === "number") {
                    if (!isColumnNumeric(c)) {
                        continue;
                    }
                } else if (filter === "string") {
                    if (!isColumnText(c)) {
                        continue;
                    }
                } else if (filter !== "all" && filter !== c.datatype) {
                    continue;
                }
            }
            //subgroup columns separate
            if (has_sgs) {
                let has = false;
                for (const s of sgs) {
                    if (c.field.startsWith(`${s}|`)) {
                        const sg = this.subgroups[s];
                        sgDataSources[s].push({
                            name: c.name,
                            field: c.field,
                            datatype: c.datatype,
                            subgroup: {
                                dataSource: sg.dataSource,
                                name_column: sg.name_column,
                            },
                        });
                        has = true;
                        break;
                    }
                }
                if (has) {
                    continue;
                }
            }
            columns.push({
                name: c.name,
                field: c.field,
                datatype: c.datatype,
            });
        }
        let cols = columns.sort((a, b) => a.name.localeCompare(b.name));
        if (has_sgs) {
            for (const ds in sgDataSources) {
                sgDataSources[ds].sort((a, b) => a.name.localeCompare(b.name));
                cols = cols.concat(sgDataSources[ds]);
            }
        }
        if (addNone) {
            cols.push({ name: "None", field: "__none__", datatype: "double" });
        }
        return cols;
    }

    /**
     * Creates and returns a dimension that it used to filter and group the data
     * @param {string} type - the dimension type , the built-in dimensions are
     * 'category_dimension' for text  columns and 'range_dimension' for number
     * columns
     * @returns {Dimension} A dimension that can be used for grouping/filtering
     */
    getDimension(type) {
        if (!Dimension.types[type]) {
            throw `Adding non existent Dimension: ${type}`;
        }
        const dim = new Dimension.types[type](this);
        this.dimensions.push(dim);

        return dim;
    }

    /**
     * Returns an object, representing the row/item containing key/value pairs
     * for all columns. As an object is created, this method is slow,
     * so it is advisable not to use it for many rows at once.
     * @param {number} index - The index of the row
     * @param {string[]?} columns - The ids of columns to use. If not provided, `columnsWithData` will be used.
     * @returns {Record<string, object>} An object containing key(field)/value pairs. An extra
     * variable 'index' containing the row index is also added to the object
     */
    getRowAsObject(index, columns) {
        if (!columns) {
            columns = this.columnsWithData;
        }
        const obj = {}; // pjt consider using Map if this is a bottleneck
        for (const c of columns) {
            obj[c] = this.columnIndex[c].getValue(index);
        }
        obj["__index__"] = index;
        return obj;
    }

    /**
     * Returns the index of the first filtered item - slow
     * @returns the index of the first filtered item
     */
    getFirstFilteredIndex() {
        for (let n = 0; n < this.size; n++) {
            if (this.filterArray[n] === 0) {
                return n;
            }
        }
    }

    _filteredIndicesPromise = null;
    _filteredIndexWorker = new Worker(
        new URL("./filteredIndexWorker.ts", import.meta.url),
    );
    /**
     * Get an array of indexes of all the items that are not filtered.
     * Subsequent calls will return a cached version of a reference to the promised array.
     * This cached reference is invalidated by any event on this DataStore
     * (and will trigger mobx/react when re-assigned).
     * @returns {Promise<Uint32Array>} A promise that resolves to an array of indexes
     */
    async getFilteredIndices() {
        if (this._filteredIndicesPromise) {
            // will be null if the filter has changed since the last call
            return this._filteredIndicesPromise;
        }
        action(
            () =>
                (this._filteredIndicesPromise = new Promise(
                    (resolve, reject) => {
                        const worker = this._filteredIndexWorker;
                        const byteLength =
                            this.filterSize * Uint32Array.BYTES_PER_ELEMENT;
                        const outputBuffer = new SharedArrayBuffer(byteLength);
                        // todo atomics... not really needed for this, but would be good to illustrate / learn.
                        worker.postMessage({
                            filterBuffer: this.filterBuffer,
                            outputBuffer,
                        });
                        worker.onmessage = (e) => {
                            resolve(new Uint32Array(outputBuffer));
                        };
                    },
                )),
        )();

        return this._filteredIndicesPromise;
    }

    /**
     * Return an array of containing all the filtered values for
     * the specified column - inefficient for large data sets
     * @param {string} column - the column's field.id
     * @returns {string[]|number[]}  An array of filtered values
     */
    getFilteredValues(column) {
        const arr = new Array(this.filterSize);
        let index = 0;
        for (let n = 0; n < this.size; n++) {
            if (this.filterArray[n] === 0) {
                arr[index] = this.getRowText(n, column);
                index++;
            }
        }
        return arr;
    }

    /**
     * Returns the value for the given row index and column
     * @param {number} index - the index of the row,
     * @param {string} column - the columns's field/id
     * @returns {string|number} - the value for the given index and field
     */
    getRowText(index, column) {
        //not very efficent
        return this.getRowAsObject(index, [column])[column];
    }

    /**
     * Sets the columns offsets
     * @param {object} data - information about offset/rotation (should only contain offsets or rotation)
     * @param {string} [data.filter] - if the offsets have a background filter - the value for this filter
     * @param {string} data.group - the group (category) to offset
     * @param {number} [data.rotation] - The amount to rotate (in degrees)
     * @param {number[]} [data.offset] - The x, y offset values
     * @param {boolean} update - update all dependants - default is false
     */
    setColumnOffset(data, update) {
        if (!this.offsets) {
            throw new Error(
                "Attempting to set offsets but none are specified in config",
            );
        }
        data.filter = data.filter || "all";
        //give default values to this group if it has none
        const o = this.offsets.values;
        let fg = o[data.filter];
        if (!fg) {
            fg = {};
            o[data.filter] = fg;
        }
        let inf = fg[data.group];
        if (!inf) {
            inf = { offsets: [0, 0], rotation: 0 };
            fg[data.group] = inf;
        }
        this.updateColumnOffsets(data, data.rotation, update);
        this.dirtyMetadata.add("offsets");
    }

    //gets the rotation origin of the group of points i.e. the center
    _getRotationOrigin(column, group, filter) {
        const o = this.offsets;
        const gc = this.columnIndex[column];
        const gv = gc.values.indexOf(group);
        const gd = gc.data;
        const x = this.columnIndex[o.param[0]];
        const y = this.columnIndex[o.param[1]];

        //get background filters if any
        let fv = null;
        let fc = null;
        if (o.background_filter) {
            fc = this.columnIndex[o.background_filter];
            fv = fc.values.indexOf(filter);
        }

        const mmx = [Number.MAX_VALUE, Number.MIN_VALUE];
        const mmy = [Number.MAX_VALUE, Number.MIN_VALUE];

        for (let n = 0; n < this.size; n++) {
            if (fc && fc.data[n] !== fv) {
                continue;
            }
            if (gd[n] === gv) {
                mmx[0] = Math.min(mmx[0], x.originalData[n]);
                mmx[1] = Math.max(mmx[1], x.originalData[n]);
                mmy[0] = Math.min(mmy[0], y.originalData[n]);
                mmy[1] = Math.max(mmy[1], y.originalData[n]);
            }
        }
        const mid = (x) => x[0] + (x[1] - x[0]) / 2;
        return [mid(mmx), mid(mmy)];
    }

    initiateOffsets() {
        const o = this.offsets;
        const x = this.columnIndex[o.param[0]];
        const y = this.columnIndex[o.param[1]];
        //store the original values in the column object
        x.originalData = new Float32Array(x.data);
        y.originalData = new Float32Array(y.data);
    }

    //gets the offset values
    //single - the only value to offset or rotate otherwise
    //the values are obtained from the offsets config
    //rotate true or false- to rotate if not then offset
    _getOffsetValues(single, rotation) {
        const o = this.offsets;
        let filterData = null;
        let filterValues = null;
        if (o.background_filter) {
            const c = this.columnIndex[o.background_filter];
            filterData = c.data;
            filterValues = c.values;
        }
        const gc = this.columnIndex[o.groups];
        const values = [];

        for (const fv in o.values) {
            if (single && fv !== single.filter) {
                continue;
            }

            for (const i in o.values[fv]) {
                if (single && i !== single.group) {
                    continue;
                }
                const v = o.values[fv][i];

                const val = {
                    index: gc.values.indexOf(i),
                };
                if (fv !== "all") {
                    val.filterData = filterData;
                    val.filterValue = filterValues.indexOf(fv);
                }

                if (rotation) {
                    if (!v.rotation_center) {
                        v.rotation_center = this._getRotationOrigin(
                            o.groups,
                            i,
                            fv,
                        );
                    }
                    val.rotation_center = [
                        v.rotation_center[0] + v.offsets[0],
                        v.rotation_center[1] + v.offsets[1],
                    ];
                }

                if (single) {
                    if (rotation) {
                        val.rotation = single.rotation;
                        v.rotation += single.rotation;
                    } else {
                        v.offsets[0] += single.offsets[0];
                        v.offsets[1] += single.offsets[1];
                        val.offsets = single.offsets;
                    }
                } else {
                    if (rotation) {
                        val.rotation = v.rotation;
                    } else {
                        val.offsets = v.offsets;
                    }
                }
                if (rotation) {
                    const radians = (Math.PI / 180) * val.rotation;
                    val.cos = Math.cos(radians);
                    val.sin = Math.sin(radians);
                }
                values.push(val);
            }
        }
        return values;
    }

    /**
     * resets the columns offsets to default values
     * @param {string|null} filter - the filter value - or null if no filter
     * @param {string} group - The group/category to reset
     * @param {boolean} update - whether to inform listeners data has changed
     */
    resetColumnOffsets(filter, group, update) {
        filter = filter || "all";
        const o = this.offsets;
        if (!o) {
            throw new Error(
                "Attempting to reset offsets but none are specified in config",
            );
        }
        delete this.offsets.values[filter][group];
        let fc = null;
        let fv = null;
        if (filter !== "all") {
            fc = this.columnIndex[this.offsets.background_filter];
            fv = fc.values.indexOf(filter);
        }
        const gr = this.getColumnValues(o.groups).indexOf(group);
        const grData = this.columnIndex[o.groups].data;
        const x = this.columnIndex[o.param[0]];
        const y = this.columnIndex[o.param[1]];

        for (let n = 0; n < this.size; n++) {
            if (fc && fc[n] !== fv) {
                continue;
            }
            if (grData[n] === gr) {
                y.data[n] = y.originalData[n];
                x.data[n] = x.originalData[n];
            }
        }
        if (update) {
            this.dataChanged([o.param[0], o.param[1]], false);
        }
        this.dirtyMetadata.add("offsets");
    }

    /**
     * Updates the columns offsets
     * @param {object} single - information about the group to rotate/offset
     * if `null` then all groups will be offset according to the values in offsets
     * @param {boolean|number} rotation -  will rotate rather than translate if "truthy"
     * @param {boolean} update - send message to update all dependants
     */
    updateColumnOffsets(single, rotation = false, update = false) {
        const o = this.offsets;
        const gc = this.columnIndex[o.groups];
        const groupData = gc.data;
        const x = this.columnIndex[o.param[0]];
        const y = this.columnIndex[o.param[1]];
        const values = this._getOffsetValues(single, rotation);

        for (let n = 0; n < this.size; n++) {
            for (const v of values) {
                if (v.filterData && v.filterData[n] !== v.filterValue) {
                    continue;
                }
                if (groupData[n] === v.index) {
                    if (rotation) {
                        const cx = v.rotation_center[0];
                        const cy = v.rotation_center[1];
                        const nx =
                            v.cos * (x.data[n] - cx) +
                            v.sin * (y.data[n] - cy) +
                            cx;
                        const ny =
                            v.cos * (y.data[n] - cy) -
                            v.sin * (x.data[n] - cx) +
                            cy;
                        x.data[n] = nx;
                        y.data[n] = ny;
                    } else {
                        x.data[n] = x.data[n] + v.offsets[0];
                        y.data[n] = y.data[n] + v.offsets[1];
                    }
                }
            }
        }
        if (update) {
            //update listeners but although column data has changed - don't synch with backend
            this.dataChanged([o.param[0], o.param[1]], false);
        }
    }

    /**
     * Gets an object whose keys contains the supplied column's values pointing
     * to the row's index. Only for text and unique columns. Will give unpredictable
     * results if a value is present more than once in the column.
     * The index is cached, so calling multiple times will not affect performance.
     * @param {string[]} column The column's field/id
     * @returns {object} An object whose keys are the columns values which point
     * to the row's index
     */
    getColumnIndex(column) {
        if (this.indexes[column]) {
            return this.indexes[column];
        }
        const col = this.columnIndex[column];
        const index = {};
        if (col.datatype === "unique") {
            for (let n = 0; n < this.size; n++) {
                let v = this.textDecoder
                    .decode(
                        col.data.slice(
                            n * col.stringLength,
                            n * col.stringLength + col.stringLength,
                        ),
                    )
                    .replaceAll("\0", "");
                //needed for some cell ids
                if (v.includes("#")) {
                    v = v.split("#")[1];
                }
                index[v] = n;
            }
        } else {
            for (let n = 0; n < this.size; n++) {
                index[col.values[col.data[n]]] = n;
            }
        }
        this.indexes[column] = index;
        return index;
    }

    /**
     * @deprecated use dataExportUtils `getExportCsvStream` instead
     * Gets a text file blob of the data which can be downloaded
     * @param {string[]} columns A list of column fields/ids to create the file with
     * @param {string|number[]} rows a value of 'filtered'  will only add filtered rows,
     * alternatively a list of row indexes can be given. Any other value will result
     * @param {string} [delimiter='\t'] The column delimiter to use
     * @param {string} [delimiter='\n'] The newline delimiter default '\n'
     * @returns {Blob} A data blob which can be downloaded
     */
    getDataAsBlob(columns, rows = "all", delimiter = "\t", newline = "\n") {
        let indexes = null;
        if (typeof rows !== "string") {
            indexes = rows;
        }

        const arr = [];
        const headers = ["index"].concat(columns);
        const len = indexes ? indexes.length : this.size;

        arr.push(headers.join(delimiter));
        const chunkSize = 10000; // You can adjust this depending on performance
        const dataChunks = [];

        for (let i = 0; i < len; i++) {
            let index = i;
            if (indexes) {
                index = indexes[i];
            } else {
                if (rows === "filtered" && this.filterArray[i] === 1) {
                    continue;
                }
            }
            const o = this.getRowAsObject(index, columns);
            const line = [i].concat(columns.map((x) => o[x]));
            arr.push(line.join(delimiter));

            // Flush the data in chunks (doesn't help much, see getExportCsvStream for something more scalable)
            if (arr.length >= chunkSize) {
                dataChunks.push(arr.join(newline));
                arr.length = 0; // Reset the array
            }
        }

        if (arr.length > 0) {
            dataChunks.push(arr.join(newline)); // Add the final chunk
        }

        // Join all the chunks into a single blob
        return new Blob([dataChunks.join(newline)], { type: "text/plain" });
    }

    /**
     * Adds data to the store for the specified column. If the data is a JavaScript array,
     * then it will be converted to the correct internal data structures and the only previously
     * supplied metadata required are field, name and datatype.
     * If the data is a shared array buffer then, the data should be in the correct
     * format see [columns](../../docs/extradocs/datasource.md)
     * @param {string}  column - The field/id of the column.
     * @param {SharedArrayBuffer|Array} data  either a javascript array or shared array buffer
     */
    setColumnData(column, data) {
        // if (Math.random() < 0.1) throw new Error("This is a test error");
        const c = this.columnIndex[column];
        if (!c) {
            throw `column ${column} is not present in data store`;
        }
        let buffer = null;
        if (Array.isArray(data)) {
            buffer = this._convertColumn(c, data);
        } else {
            buffer = data;
        }
        c.buffer = buffer;
        if (
            c.datatype === "integer" ||
            c.datatype === "double" ||
            c.datatype === "int32"
        ) {
            const dataArray = (c.data =
                c.datatype === "int32"
                    ? new Int32Array(buffer)
                    : new Float32Array(buffer));
            if (!c.minMax) {
                console.log(`Calculating min/max for ${c.name}`);
                let min = Number.MAX_VALUE;
                let max = Number.MIN_VALUE;
                for (let i = 0; i < dataArray.length; i++) {
                    const value = dataArray[i];
                    if (Number.isNaN(value)) {
                        continue;
                    }
                    min = value < min ? value : min;
                    max = value > max ? value : max;
                }
                c.minMax = [min, max];
            }
            if (!c.quantiles) {
                const a = c.data.filter((x) => !Number.isNaN(x)).sort();
                c.quantiles = {};
                for (const q of [0.05, 0.01, 0.001]) {
                    c.quantiles[q] = [
                        quantileSorted(a, q),
                        quantileSorted(a, 1 - q),
                    ];
                }
            }
        } else if (c.datatype === "multitext" || c.datatype === "text16") {
            c.data = new Uint16Array(buffer);
        } else {
            c.data = new Uint8Array(buffer);
        }
        this.columnsWithData.push(column);
    }

    //experimental
    appendColumnData(column, data, newSize) {
        const c = this.columnIndex[column];
        // handy utility type - ideally we wouldn't import this from luma.gl
        /** @type {import("@luma.gl/core").TypedArrayConstructor} */
        let arrType = Uint8Array;
        let size = newSize;
        if (c.datatype === "integer" || c.datatype === "double") {
            size = size * 4;
            arrType = Float32Array;
        } else if (c.datatype === "int32") {
            size = size * 4;
            arrType = Int32Array;
        } else if (c.datatype === "unique") {
            size = size * c.stringLength;
        }

        const newBuffer = new SharedArrayBuffer(size);
        const newArr = new arrType(newBuffer);
        newArr.set(c.data);
        newArr.set(new arrType(data), c.data.length);
        c.data = newArr;
        c.buffer = newBuffer;
    }

    //experimental
    addDataToStore(columnData, size) {
        const newSize = size + this.size;
        for (const c of columnData) {
            this.appendColumnData(c.column, c.data, newSize);
        }
        const newBuffer = new SharedArrayBuffer(newSize);
        const newArr = new Uint8Array(newBuffer);
        newArr.set(this.filterArray);

        this.filterBuffer = newBuffer;
        this.filterArray = newArr;
        this.size = newSize;
        this.filterSize += size;
        //update dimensions and redo any filters
        for (const d of this.dimensions) {
            d.updateSize();
        }
        this._callListeners("data_added", this.size);
    }

    //experimental
    cleanColumnData(column) {
        const index = {};
        const col = this.columnIndex[column];
        for (const v in col.values) {
            index[v] = 0;
        }
        for (let i = 0; i < this.size; i++) {
            index[col.data[i]]++;
        }
        const newVals = [];
        const oldToNew = {};
        for (let i = 0; i < col.values.length; i++) {
            if (index[i] !== 0) {
                newVals.push(col.values[i]);
                oldToNew[i] = newVals.length - 1;
            }
        }
        for (let i = 0; i < this.size; i++) {
            col.data[i] = oldToNew[col.data[i]];
        }
        col.values = newVals;
        this.setColumnIsDirty(column);
        this.dataChanged([column]);
    }

    //converts a JavaScript array to the intenal data structure (Typed Array)
    //adding any mentadata requires e.g values. stringLength
    _convertColumn(col, arr) {
        const len = arr.length;
        //pjt: note - looked like a bug here with `if (col.datatype==="text" || "text16")` - always truthy
        //also, should this handle "multitext" as well? No - that has a separate case below.
        //this was the source of bug with saving multitext columns for __tags
        if (col.datatype === "text" || col.datatype === "text16") {
            const t8 = col.datatype === "text";
            const buff = new SharedArrayBuffer(t8 ? this.size : this.size * 2);
            //work out number or rows with each value
            const v_to_n = {};
            for (let i = 0; i < len; i++) {
                const v = arr[i];
                if (v_to_n[v] === undefined) {
                    v_to_n[v] = 1;
                } else {
                    v_to_n[v]++;
                }
            }
            //sort values by number of rows
            const li = [];
            for (const v in v_to_n) {
                li.push([v, v_to_n[v]]);
            }
            col.values = [];
            //dictionary to convert value to index
            const v_to_i = {};
            li.sort((a, b) => b[1] - a[1]);
            for (let i = 0; i < li.length; i++) {
                col.values.push(li[i][0]);
                v_to_i[li[i][0]] = i;
            }

            const a = t8 ? new Uint8Array(buff) : new Uint16Array(buff);
            for (let i = 0; i < len; i++) {
                a[i] = v_to_i[arr[i]];
            }
            return buff;
        }
        if (
            col.datatype === "integer" ||
            col.datatype === "double" ||
            col.datatype === "int32"
        ) {
            const buff = new SharedArrayBuffer(this.size * 4);
            const a =
                col.datatype === "int32"
                    ? new Int32Array(buff)
                    : new Float32Array(buff);

            for (let i = 0; i < len; i++) {
                a[i] = arr[i];
            }
            return buff;
        }
        if (col.datatype === "multitext") {
            const delim = col.delimiter || ",";
            const vals = new Set();
            let max = 0;
            //first parse - get all possible values and max number
            //of values in a single field
            for (let i = 0; i < len; i++) {
                const v = arr[i];
                const vs = v.split(",");
                max = Math.max(max, vs.length);
                vs.forEach((x) => vals.add(x.trim()));
            }
            vals.delete("");
            const buff = new SharedArrayBuffer(this.size * max * 2);
            const data = new Uint16Array(buff);
            data.fill(65535);
            const values = new Array(vals.size);
            //more efficent than using indexOf in array
            const map = {};
            let index = 0;
            for (const v of vals) {
                map[v] = index;
                values[index] = v;
                index++;
            }

            for (let i = 0; i < len; i++) {
                const b = i * max;
                const v = arr[i];
                if (v === "") {
                    continue;
                }
                const vs = v.split(",");
                vs.sort();
                for (let n = 0; n < vs.length; n++) {
                    data[b + n] = map[vs[n].trim()];
                }
            }
            col.values = values;
            col.stringLength = max;
            return buff;
        }

        let max = 0;
        for (let i = 0; i < len; i++) {
            max = Math.max(max, arr[i].length);
        }
        col.stringLength = max;
        const enc = new TextEncoder();
        const buff = new SharedArrayBuffer(this.size * max);
        for (let i = 0; i < len; i++) {
            const a = enc.encode(arr[i].substring(0, max));
            const b = new Uint8Array(buff, i * max, max);
            b.set(a, 0);
        }

        return buff;
    }

    /**
     * Returns a function which gives the appropriate color for the value of
     * the specified column, when supplied with the index of a row/item in the datastore,
     * @param {string|Partial<Column>} column The column id(field) to use for the function . Can be null, if
     * useValue = true i.e.  the functions will be supplied with a value rather than an index
     * @param {object} [config] An optional config with extra parameters
     * @param {number} [config.bins=100] For columns with continuous data (doubles/integers),
     * bins are calculated across the data range so that only a limited number of
     * color values need to be calculated. The default is 100, although it can be
     * altered here.
     * @param {boolean} [config.asArray=false] By default the, function will return
     * a JavaScript compatible string specifying the color. If asArray is true then an array
     * of length 3 containing RGB values will be returned.
     * @param {object} [config.overideValues] an object containing values to use, instead of
     * the columns default values - can include
     * <ul>
     * <li> min - the minumum value </li>
     * <li> max - the maximum value </li>
     * <li> colors - the color scheme to use </li>
     * <li> colorLogScale- whether to use a log scale </li>
     * <li> fallbackOnZero- whether to return the fallback color for zeros, as well as NaN </li>
     * </ul
     * @param {boolean} [config.useValue=false]  The returned function will require the
     * columns value, not index
     * @param {string} [config.datatype] Only required if a column is not supplied
     * @returns {(i: number) => [number, number, number]} The function, which when given a row index (or value if this
     * was specified) will return a color.
     */
    getColorFunction(column, config = {}) {
        let c = this.columnIndex[column];
        //not actually related to column -add dummy column
        if (column == null) {
            //wtf is this? how is `column` used after this?
            c = column = {
                datatype: config.datatype,
            };
            //"the returned function will require the columns value, not index"
            //clear as mud...
            config.useValue = true;
        }

        const data = c.data;
        const ov = config.overideValues || {};
        const colors = this.getColumnColors(column, config);
        function isFallback(v) {
            return !Number.isFinite(v) || (ov.fallbackOnZero && v === 0);
        }
        //simply return the color associated with the value
        if (
            c.datatype === "text" ||
            c.datatype === "text16" ||
            c.datatype === "multitext"
        ) {
            return (x) => colors[data[x]];
        }
        if (
            c.datatype === "integer" ||
            c.datatype === "double" ||
            c.datatype === "int32"
        ) {
            const min = ov.min == null ? c.minMax[0] : ov.min;
            const max = ov.max == null ? c.minMax[1] : ov.max;
            const bins = config.bins || 100;
            const interval_size = (max - min) / bins;
            const fallbackColor = config.asArray ? [255, 255, 255] : "#ffffff";
            //the actual function - bins the value and returns the color for that bin
            function getColor(v) {
                if (isFallback(v)) return fallbackColor;
                // interval_size will be 0 for e.g. column containing only zeros.
                // don't want to use fallbackColor, just the zero-th color
                if (interval_size <= Number.EPSILON) return colors[0];
                let bin = Math.floor((v - min) / interval_size);
                if (bin < 0) {
                    bin = 0;
                } else if (bin >= colors.length) {
                    bin = colors.length - 1;
                }
                /// we could go belt-and-braces on error checking here, but should be robust now and this will have some overhead.
                // if (!Number.isFinite(bin)) {
                //     console.error('numerical issue computing color bin');
                //     bin = 0;
                // }
                return colors[bin];
            }
            if (config.useValue) {
                return getColor;
            }
            return x => getColor(data[x]);
        }
    }

    /**
     * For a given column will returns the given
     * color for that category
     * @param {string} column - the coloumn's field/id
     * @param {string} cat - the category
     * @returns {string|number[]} - the hex value of the color
     */
    getColorForCategory(column, cat) {
        const c = this.columnIndex[column];
        const i = c.values.indexOf(cat);
        return this.getColumnColors(column)[i];
    }

    /**
     * Makes a color bar/legend based on the give column
     * @param {string} column - the field/id of the column
     * @param {object} config - see [here]{@link DataStore#getColorFunction}
     * @returns {HTMLElement} - a color bar or color legend
     */
    getColorLegend(column, config = {}) {
        const colors = this.getColumnColors(column, config);
        const c = this.columnIndex[column];
        const name = config.name || c.name;
        if (
            c.datatype === "integer" ||
            c.datatype === "double" ||
            c.datatype === "int32"
        ) {
            const [min, max] = this.getMinMaxForColumn(column);
            let range = [min, max];
            if (config.overideValues) {
                const ov = config.overideValues;
                range = [
                    ov.min == null ? min : ov.min,
                    ov.max == null ? max : ov.max,
                ];
            }
            return getColorBar(colors, { range: range, label: name });
        }
        if (isColumnText(c)) {
            return getColorLegend(colors, c.values, { label: name });
        }
    }

    /**
     * This method returns an object whose keys are categories
     * and values are the colors of the categories
     * @param {string} column The column id(field) (only text/multitext columns)
     * @returns {object} an object of values to colors
     */
    getValueToColor(column) {
        const vc = {};
        const colors = this.getColumnColors(column);
        const values = this.getColumnValues(column);
        for (let i = 0; i < values.length; i++) {
            vc[values[i]] = colors[i];
        }
        return vc;
    }

    /**
     * Returns the n and 1-n qauntile values where n is the given percentile
     * If no percentile is passed or is 'all'  than the column's max/min
     * value will be returned
     * @param {string} column - the column's field/id
     * @param {string}  [per] - the percentile - either 0.001.0.01 or 0.05
     * @returns {number[]} - the n and 1-n percentiles
     */
    getColumnQuantile(column, per) {
        const col = this.columnIndex[column];
        if (per && per !== "none") {
            if (col.quantiles && col.quantiles !== "NA") {
                return [col.quantiles[per][0], col.quantiles[per][1]];
            }
        } else {
            return col.minMax;
        }
    }

    /**
     * Changes the columns current colors
     * @param {string} column - the column's field or id
     * @param {string[]} colors - the column's new colors
     * @param {boolean} [notify=false] - notify any listeners that
     * the data has changed
     */
    setColumnColors(column, colors, notify = false) {
        const c = this.columnIndex[column];
        if (colors.length !== c.values.length) {
            throw new Error(
                `${column} being set with incorrect number of columns `,
            );
        }
        c.colors = colors.slice(0);
        this.dirtyColumns.colors_changed[column] = true;
        if (notify) {
            this.dataChanged([column], false);
        }
    }

    /**
     * Returns the colors that the column uses or default ones if none have been set
     * @param {string|Partial<Column>} column - the column's field/id, or a column object
     * @param {object} config - see [here]{@link DataStore#getColorFunction}
     * @returns {(string|number[])[]} An array of colors. For text/multitext, the colors will correspond
     * to the column's values. For double/integers it will contain binned values from the min
     * to max value
     */
    getColumnColors(column, config = {}) {
        /** @type {Column} */
        const c = typeof column === "string" ? this.columnIndex[column] : column;
        const rArr = config.asArray;
        if (
            c.datatype === "double" ||
            c.datatype === "integer" ||
            c.datatype === "int32"
        ) {
            const ov = config.overideValues || {};
            const c_colors = ov.colors || c.colors || defaultIPalette;
            const min = ov.min === undefined ? c.minMax[0] : ov.min;
            const max = ov.max === undefined ? c.minMax[1] : ov.max;
            //calculate the color of each bin
            const ls = linspace(min, max, c_colors.length);
            const scale = scaleLinear().domain(ls).range(c_colors).clamp(true);
            const bins = ov.bins || 100;
            const interval_size = (max - min) / bins;
            let colors = new Array(bins + 1);
            for (let i = 0; i < bins + 1; i++) {
                colors[i] = scale(min + i * interval_size);
                //convert to rgb array
                if (rArr) {
                    colors[i] = rgbToRGB(colors[i]);
                }
            }
            let useLog = c.colorLogScale;
            if (ov.colorLogScale != null) {
                useLog = ov.colorLogScale;
            }
            if (useLog) {
                //calculate new colors based on a log scale
                const logScale = scaleSymlog()
                    .domain([min, max])
                    .range([0, bins])
                    .clamp(true);
                const newColors = new Array(bins + 1);
                for (let i = 0; i < bins + 1; i++) {
                    newColors[i] =
                        colors[Math.floor(logScale(min + i * interval_size))];
                }
                colors = newColors;
            }
            return colors;
        }
        if (
            c.datatype === "text" ||
            c.datatype === "multitext" ||
            c.datatype === "text16"
        ) {
            let colors = c.colors;
            if (!colors) {
                const vlen = c.values.length;
                const dlen = defaultPalette.length;
                if (vlen < dlen) {
                    colors = defaultPalette.slice(0, vlen);
                } else {
                    colors = [];
                    const times = Math.floor(vlen / dlen);
                    for (let n = 0; n < times; n++) {
                        colors = colors.concat(defaultPalette);
                    }
                    colors = colors.concat(
                        defaultPalette.slice(0, vlen % dlen),
                    );
                }
            }

            if (rArr) {
                colors = colors.map((x) => hexToRGB(x));
            }
            return colors;
        }
    }

    addColumnGroup(group) {
        if (!group.columns) {
            group.columns = [];
        }
        this.columnGroups[group.name] = group;
    }

    getColumnGroup(name) {
        const cg = this.columnGroups[name];
        if (!cg) {
            return null;
        }
        return cg.columns;
    }

    getRawColumn(column) {
        return this.columnIndex[column].data;
    }

    /**
     * Returns the min/max values for a given column
     * @param {string} column The column id(field). Should be a numeric column, otherwise an error will be thrown
     * @returns {[number, number]} An array - the first value being the min value and the second the max value
     */
    getMinMaxForColumn(column) {
        const c = this.columnIndex[column];
        if (!isDatatypeNumeric(c.datatype)) {
            throw new Error(`Trying to get minMax for non-numeric column '${column}'`);
        }
        // columns loaded via rows_as_columns_link may not automatically have this metadata, so we can compute here.
        if (!c.minMax) {
            if (!c.data) {
                throw new Error(`Attempting to compute minMax for column '${column}' which is not loaded...`);
                //considered returning some default vals here
                // return [0, 1];
            }
            // Return the calculated min and max values
            let min = Number.MAX_VALUE;
            let max = Number.MIN_VALUE;
            for (let i = 0; i < c.data.length; i++) {
                const value = c.data[i];
                if (Number.isNaN(value)) {
                    continue;
                }
                min = value < min ? value : min;
                max = value > max ? value : max;
            }
            c.minMax = [min, max]
            return c.minMax;
        }
        return c.minMax;
    }

    getColumnRange(column) {
        const [min, max] = this.getMinMaxForColumn(column);
        return max - min;
    }

    /**
     * @param {string} column if this is in the form of a link `FieldName` it will first attempt to
     * `addColumnFromField()`...
     */
    getColumnInfo(column) {
        if (column === undefined) {
            throw new Error("getColumnInfo needs a column argument!");
        }
        if (column.includes("|")) {
            this.addColumnFromField(column);
        }
        const c = this.columnIndex[column];
        if (!c) throw `no metadata for column '${column}' in ds '${this.name}'`;
        return {
            name: c.name,
            field: c.field,
            datatype: c.datatype,
            stringLength: c.stringLength,
            subgroup: c.subgroup,
            sgindex: c.sgindex,
            sgtype: c.sgtype,
        };
    }

    getLoadedColumns() {
        return this.columns.filter((x) => x.data != null).map((x) => x.field);
    }

    getAllColumns() {
        return this.columns.map((x) => x.field);
    }

    /**
     * @param {string | undefined} column - the column's field/id
     * @param {"name_value" | undefined} format - if "name_value", returns `{name: string, value: string}[]` for use in settings 'dropdown' widget
     * @returns {string[]} - the column's values
     */
    getColumnValues(column, format = null) {
        // if the column is not categorical, it will return undefined - that may be ok
        // if the column is a linked field, there is a chance that it won't be in the columnIndex
        // that is a problem because then we throw an error here.
        if (column === undefined) {
            console.warn('getColumnValues(undefined)');
            return [];
        }
        const v = this.columnIndex[column]?.values;
        if (!v) {
            console.error(`no values for column '${column}' in ds '${this.name}'`);
            return [];
        }
        // could throw here if v is undefined (ie bad column arg)
        if (format === "name_value") {
            const ls = Array.from(v, (x) => ({ name: x, value: x })).sort(
                (a, b) => a.name.localeCompare(b.name),
            );
            return ls;
        }
        return v;
    }

    getColumnName(col) {
        const c = this.columnIndex[col];
        return c ? c.name : null;
    }

    /**
     * Removes the column and all its data
     * @param {string} column - the columns id/field
     * @param {boolean} [dirty=false] if true, tags that the column should also be removed from the
     * backend
     * @param {boolean} [notify=false] if true notifies any listeners that the column has been removed
     */

    removeColumn(column, dirty = false, notify = false) {
        const c = this.columnIndex[column];
        c.data = null;
        c.buffer = null;
        this.columns = this.columns.filter((c) => c.field !== column);
        delete this.columnIndex[column];
        const i = this.columnsWithData.indexOf(column);
        if (i !== -1) {
            this.columnsWithData.splice(i, 1);
        }
        if (dirty) {
            //added and removed without saving
            if (this.dirtyColumns.added[column]) {
                delete this.dirtyColumns.added[column];
            } else {
                this.dirtyColumns.removed[column] = true;
            }
        }
        if (notify) {
            this._callListeners("column_removed", column);
        }
    }
}

function linspace(start, end, n) {
    const out = [];
    const delta = (end - start) / (n - 1);

    let i = 0;
    while (i < n - 1) {
        out.push(start + i * delta);
        i++;
    }

    out.push(end);
    return out;
}

function hexToRGB(hex) {
    hex = hex.replace("#", "");
    const bigint = Number.parseInt(hex, 16);
    const r = (bigint >> 16) & 255;
    const g = (bigint >> 8) & 255;
    const b = bigint & 255;
    return [r, g, b];
}

function RGBToHex(rgb) {
    return `#${rgb
        .map((x) => {
            const hex = x.toString(16);
            return hex.length === 1 ? `0${hex}` : hex;
        })
        .join("")}`;
}
function rgbToRGB(rgb) {
    if (!rgb) {
        return [255, 255, 255];
    }
    rgb = rgb.substring(4, rgb.length - 1).split(", ");
    return rgb.map((x) => Number.parseInt(x));
}

const defaultPalette = [
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    "#ffff33",
    "#a65628",
    "#f781bf",
    "#999999",
    "#1CE6FF",
    "#FF34FF",
    "#FF4A46",
    "#008941",
    "#676FA6",
    "#A30059",
    "#FFDBE5",
    "#7A4900",
    "#0000A6",
    "#63FFAC",
    "#B79762",
    "#004D43",
    "#8FB0FF",
    "#997D87",
    "#5A0007",
    "#809693",
    "#FEFFE6",
    "#1B4400",
    "#4FC601",
    "#3B5DFF",
    "#4A3B53",
    "#FF2F80",
    "#61615A",
    "#BA0900",
    "#6B7900",
    "#00C2A0",
    "#FFAA92",
    "#FF90C9",
    "#B903AA",
    "#D16100",
    "#DDEFFF",
    "#000035",
    "#7B4F4B",
    "#A1C299",
    "#300018",
    "#0AA6D8",
    "#013349",
    "#00846F",
    "#372101",
    "#FFB500",
    "#C2FFED",
    "#A079BF",
    "#CC0744",
    "#C0B9B2",
    "#C2FF99",
    "#001E09",
    "#00489C",
    "#6F0062",
    "#0CBD66",
    "#EEC3FF",
    "#456D75",
    "#B77B68",
    "#7A87A1",
    "#788D66",
    "#885578",
    "#FAD09F",
    "#FF8A9A",
    "#D157A0",
    "#BEC459",
    "#456648",
    "#0086ED",
    "#886F4C",
    "#34362D",
    "#B4A8BD",
    "#00A6AA",
    "#452C2C",
    "#636375",
    "#A3C8C9",
    "#FF913F",
    "#938A81",
    "#575329",
    "#00FECF",
    "#B05B6F",
    "#8CD0FF",
    "#3B9700",
    "#04F757",
    "#C8A1A1",
    "#1E6E00",
    "#7900D7",
    "#A77500",
    "#6367A9",
    "#A05837",
    "#6B002C",
    "#772600",
    "#D790FF",
    "#9B9700",
    "#549E79",
    "#FFF69F",
    "#201625",
    "#72418F",
    "#BC23FF",
    "#99ADC0",
    "#3A2465",
    "#922329",
    "#5B4534",
    "#FDE8DC",
];

const defaultIPalette = [
    "#3288bd",
    "#66c2a5",
    "#abdda4",
    "#e6f598",
    "#fee08b",
    "#fdae61",
    "#f46d43",
    "#d53e4f",
];

export default DataStore;
export { hexToRGB, defaultPalette, RGBToHex };
