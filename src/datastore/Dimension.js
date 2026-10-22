class Dimension {
    /**
     * @class
     * @param {DataStore} parent
     */
    constructor(parent) {
        this.filterBuffer = new SharedArrayBuffer(parent.size);
        this.parent = parent;
        this.noClear = false;
        this.filterArray = new Uint8Array(this.filterBuffer);
        this.filterArguments = null;
        this.filterIndexes = null;
        this.filterMethod = null;
    }

 

    _applyStateTransition(index, nextState) {
        const prevState = this.filterArray[index];
        if (prevState === nextState) {
            return;
        }
        const wasFiltered = prevState === 1 || prevState === 3;
        const isFiltered = nextState === 1 || nextState === 3;
        this.filterArray[index] = nextState;
        if (wasFiltered === isFiltered) {
            return;
        }
        const parentFilter = this.parent.filterArray;
        if (isFiltered) {
            if (++parentFilter[index] === 1) {
                this.parent.filterSize--;
            }
        } else if (--parentFilter[index] === 0) {
            this.parent.filterSize++;
        }
    }

    /**
     * Applies a filter result to this dimension's local states and the parent ref-counts.
     * Every filter method should end here, so the 0/1/2/3 state machine lives in one place:
     * background-hidden rows (2) are promoted to 3, rows already at 3 are left alone, and
     * 0 <-> 1 follows `include`.
     *
     * @param {((i: number) => boolean) | ArrayLike<number>} include - a per-row predicate, or a
     * mask (e.g. a Uint8Array produced by a worker) where non-zero keeps the row
     */
    commitFilter(include) {
        const size = this.parent.size;
        const isMask = typeof include !== "function";
        if (isMask && include.length < size) {
            throw new Error(`filter mask has ${include.length} rows, expected ${size}`);
        }
        const local = this.filterArray;
        const parentFilter = this.parent.filterArray;
        let passingDelta = 0;
        for (let i = 0; i < size; i++) {
            const prevState = local[i];
            // local brush only mutates visible rows; background rows are updated when the
            // background filter is set or cleared
            if (prevState === 3) {
                continue;
            }
            // with a local filter active, background-hidden rows are excluded too
            if (prevState === 2) {
                local[i] = 3;
                if (++parentFilter[i] === 1) passingDelta--;
                continue;
            }
            const keep = isMask ? include[i] !== 0 : include(i);
            if (keep) {
                if (prevState === 1) {
                    local[i] = 0;
                    if (--parentFilter[i] === 0) passingDelta++;
                }
            } else if (prevState === 0) {
                local[i] = 1;
                if (++parentFilter[i] === 1) passingDelta--;
            }
        }
        this.parent.filterSize += passingDelta;
    }

    /**
     * Filters the data based on a predicate function.
     * The function should return false if the row should be filtered out.
     */
    filterPredicate(args, columns) {
        return this.commitFilter(args.predicate);
    }
    /**
     * 
     * @param {*} noClear - if true, then the filter will not be cleared 
     * when a removeAllFilters is called on the parent DataStore. Default is false.
     */
    setNoClear(noClear){
        this.noClear = noClear;
    }

    removeFilter(notify = true) {
        if (!this.filterMethod) {
            return;
        }
        const localFilter = this.filterArray;
        const len = this.parent.size;
        this.filterArguments = undefined;
        this.filterColumns = undefined;
        this.filterMethod = undefined;
        for (let i = 0; i < len; i++) {
            // clear local contribution only: 1->0 and 3->2
            if (localFilter[i] === 1) {
                this._applyStateTransition(i, 0);
            } else if (localFilter[i] === 3) {
                this._applyStateTransition(i, 2);
            }
        }
   
        if (notify) {
            this.parent._callListeners("filtered", this);
        }
    }

    getLocalFilter() {
        return this.filterArray;
    }

    /**
     * Keeps only the rows in `indexSet`.
     * @param {Set<number>} indexSet
     */
    filterOnIndex(indexSet) {
        return this.commitFilter((i) => indexSet.has(i));
    }

    /**
     * Filters  the data
     * @param {string} method - The name of the filter method.
     * if 'filterPredicate' then the args should contain a predicate function
     * @param {string[]} columns - a list of column ids used in the filtering
     * @param {object|string} args - any extra arguments for filtering
     * @param {boolean} [notify=true] - notify any listeners in the dataStore that the
     * data has been filtered
     */
    filter(method, columns, args, notify = true) {
        const t = performance.now();
        this.filterColumns = columns;
        this.filterArguments = args;
        this.filterMethod = method;
        this[method](args, columns);
        if (notify) {
            this.parent._callListeners("filtered", this);
        }
        if (performance.now() - t > 100) {
            console.warn(`method${method}: ${performance.now() - t}`);
        } else {
            console.log(`method${method}: ${performance.now() - t}`);
        }
    }

    async getValueSet(column) {
        const col = this.parent.columnIndex[column];
        const data = col.data;
        const parentFilter = this.parent.filterArray;
        const filter = this.filterArray;
        const len = this.parent.size;
        const set = new Set();
        for (let i = 0; i < len; i++) {
            if (filter[i] === 0 && parentFilter[i] === 0) {
                if (col.values) {
                    set.add(col.values[data[i]]);
                } else set.add(data[i]);
            }
        }
        return set;
    }

    /**
     * updates the filter if data in the supplied columns has been altered
     * does not propagate (call any listeners)
     * @param {string[]} columns - a list of columns whose data has changes
     * @returns {boolean}  True if refiltering has taken place
     */
    reFilterOnDataChanged(columns) {
        if (this.filterMethod && this.filterColumns) {
            for (const c of this.filterColumns) {
                if (typeof c === "string" && columns.indexOf(c) !== -1) {
                    this.filter(
                        this.filterMethod,
                        this.filterColumns,
                        this.filterArguments,
                        false,
                    );
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Experimental - should be invoked if the datastore has increased in size
     * Will increase size of local filters
     */
    updateSize() {
        const newBuff = new SharedArrayBuffer(this.parent.size);
        const newArr = new Uint8Array(newBuff);
        newArr.set(this.filterArray);
        this.filterBuffer = newBuff;
        this.filterArray = newArr;
        if (this.bgfData) {
            const bd = this.bgfData;
            this.setBackgroundFilter(
                bd.column,
                bd.cats?.length ? bd.cats : bd.cat,
            );
        }
        if (this.filterMethod) {
            this.filter(
                this.filterMethod,
                this.filterColumns,
                this.filterArguments,
                false,
            );
        }
    }

    /**
     * sets a permanent filter on the chart
     * @param {string} column - The column of the filter
     * @param {string|string[]} cat - category/category list in the column to filter on
     */
    setBackgroundFilter(column, cat) {
        const col = this.parent.columnIndex[column];
        if (!col?.values) {
            return;
        }
        const cats = (Array.isArray(cat) ? cat.slice(0) : [cat]).filter(
            (x) => x !== undefined && x !== null,
        );
        if (cats.length === 0) {
            this.clearBackgroundFilter();
            return;
        }
        const indices = new Set();
        for (const c of cats) {
            const ci = col.values.indexOf(c);
            if (ci !== -1) {
                indices.add(ci);
            }
        }
        const selectedCats = new Set(cats);
        const isMultitext = col.datatype === "multitext";
        const multitextCapacity = Math.max(col.stringLength || 1, 1);
        const multitextDelimiter = col.delimiter
            ? col.delimiter
            : col.values.some((value) => value?.includes(";"))
              ? ";"
              : ",";
        const data = col.data;
        this.bgfData = {
            column: column,
            cat: Array.isArray(cat) ? cat[0] : cat,
            cats: cats,
        };
        for (let i = 0; i < this.parent.size; i++) {
            const prevState = this.filterArray[i];
            let nextState = prevState;
            let rowMatches = false;
            if (isMultitext && multitextCapacity > 1) {
                const start = i * multitextCapacity;
                for (let offset = 0; offset < multitextCapacity; offset += 1) {
                    const valueIndex = data[start + offset];
                    if (
                        valueIndex === undefined ||
                        valueIndex === null ||
                        valueIndex === 65535
                    ) {
                        break;
                    }
                    if (indices.has(valueIndex)) {
                        rowMatches = true;
                        break;
                    }
                }
            } else if (isMultitext) {
                const valueIndex = data[i];
                if (indices.has(valueIndex)) {
                    rowMatches = true;
                } else {
                    const value = col.values[valueIndex];
                    if (typeof value === "string") {
                        const items = value
                            .split(multitextDelimiter)
                            .map((item) => item.trim())
                            .filter((item) => item.length > 0);
                        rowMatches = items.some((item) => selectedCats.has(item));
                    }
                }
            } else {
                rowMatches = indices.has(data[i]);
            }

            if (rowMatches) {
                if (prevState === 2 || prevState === 3) {
                    nextState -= 2;
                }
            } else {
                if (prevState === 0) {
                    nextState = 2;
                } else if (prevState === 1) {
                    nextState = 3;
                }
            }
            this._applyStateTransition(i, nextState);
        }
        if (this.filterMethod) {
            this[this.filterMethod](this.filterArguments, this.filterColumns);
        }
        this.parent._callListeners("filtered", this);
    }

    clearBackgroundFilter() {
        if (!this.bgfData) {
            return;
        }
        this.bgfData = null;
        for (let i = 0; i < this.parent.size; i++) {
            if (this.filterArray[i] > 1) {
                this._applyStateTransition(i, this.filterArray[i] - 2);
            }
        }
        if (this.filterMethod) {
            this[this.filterMethod](this.filterArguments, this.filterColumns);
        }
        this.parent._callListeners("filtered", this);
    }

    destroy(notify = true) {
        this.removeFilter(notify);
        this.filterBuffer = null;
        this.localFilter = null;
        const index = this.parent.dimensions.indexOf(this);
        // guard against double-destroy: splice(-1, 1) would remove an unrelated dimension
        if (index !== -1) this.parent.dimensions.splice(index, 1);
    }
}
/** @type {{[k: string]: undefined | typeof Dimension}} */
Dimension.types = { base_dimension: Dimension };
export default Dimension;
