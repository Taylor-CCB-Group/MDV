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
     * Filters the data based on a predicate function.
     * The function should return false if the row should be filtered out.
     */
    filterPredicate(args, columns) {
        const predicate = args.predicate;
        // could we pass in optional filteredIndices and use those rather than processing entire array?
        // it's not as simple as resetting the localFilter because of parent.filterSize mutation
        // if we guarantee that users of a particular Dimension instance always go through this method,
        // then we can fairly easily add some other properties for it to reference
        // - this will make things like a/b comparison of predicate vs other versions of filters harder to test.
        // - suspect preferred approach will be a new strategy for filter evaluation that doesn't attempt to use this class
        const localFilter = this.filterArray;
        for (let i = 0; i < this.parent.size; i++) {
            // try ... catch to handle errors in the predicate
            // we could probably pass something monadic in args to handle this
            // let value = false;
            // try {
            //     value = !predicate(i);
            // } catch (e) {
            //     // console.error('Error in evaluating filterPredicate', e);
            // }
            //not interested in rows already in background filter
          
            const shouldExclude = !predicate(i);
            const prevState = localFilter[i];
            // When a local filter is active, background-hidden rows must also
            // contribute to exclusion scope (promote 2 -> 3 once).
            if (prevState === 2) {
                this._applyStateTransition(i, 3);
                continue;
            }
            // local brush should only mutate visible points
            if (prevState === 3) {
                continue;
            }
            if (shouldExclude) {
                if (prevState === 0) {
                    this._applyStateTransition(i, 1);
                }
            } else if (prevState === 1) {
                this._applyStateTransition(i, 0);
            }
        }
        // xxx: why would we not notify listeners?
        //I'm a bit dubious about the general pattern
        //- using filter(methodName...) for now to be more in line with other things
        //this.parent._callListeners("filtered", this);
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

    //needs to be removed - produces strange behaviour << review
    filterOnIndex(indexSet) {
        const filter = this.parent.filterArray;
        const len = this.parent.size;
        const localFilter = this.filterArray;
        const parent = this.parent;
        for (let i = 0; i < len; i++) {
            if (indexSet.has(i)) {
                if (localFilter[i] === 1) {
                    if (--filter[i] === 0) {
                        parent.filterSize++;
                    }
                }
                localFilter[i] = 0;
            } else {
                if (localFilter[i] === 0) {
                    if (++filter[i] === 1) {
                        parent.filterSize--;
                    }
                }
                localFilter[i] = 1;
            }
        }
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
            this.clearBackGroundFilter();
            return;
        }
        const indices = new Set();
        for (const c of cats) {
            const ci = col.values.indexOf(c);
            if (ci !== -1) {
                indices.add(ci);
            }
        }
        const data = col.data;
        this.bgfData = {
            column: column,
            cat: Array.isArray(cat) ? cat[0] : cat,
            cats: cats,
        };
        for (let i = 0; i < this.parent.size; i++) {
            const prevState = this.filterArray[i];
            let nextState = prevState;
            const hasLocalFilter = Boolean(this.filterMethod);
            if (indices.has(data[i])) {
                if (prevState === 2 || prevState === 3) {
                    nextState -= 2;
                }
            } else {
                if (prevState === 0) {
                    nextState = hasLocalFilter ? 3 : 2;
                } else if (prevState === 1) {
                    nextState = 3;
                }
            }
            this._applyStateTransition(i, nextState);
        }
    }

    clearBackGroundFilter() {
        if (!this.bgfData) {
            return;
        }
        this.bgfData = null;
        for (let i = 0; i < this.parent.size; i++) {
            if (this.filterArray[i] > 1) {
                this._applyStateTransition(i, this.filterArray[i] - 2);
            }
        }
    }

    destroy(notify = true) {
        this.removeFilter(notify);
        this.filterBuffer = null;
        this.localFilter = null;
        this.parent.dimensions.splice(this.parent.dimensions.indexOf(this), 1);
    }
}
/** @type {{[k: string]: undefined | typeof Dimension}} */
Dimension.types = { base_dimension: Dimension };
export default Dimension;
