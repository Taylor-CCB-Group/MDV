import Dimension from "./Dimension.js";
import {
    MULTITEXT_EMPTY_INDEX,
    getMultitextCapacity,
    getMultitextValueItems,
} from "@/lib/multitext";

function rowHasAnyMultitextItem(data, start, capacity, valueItemsByIndex, selectedItems) {
    for (let offset = 0; offset < capacity; offset++) {
        const valueIndex = data[start + offset];
        if (valueIndex === MULTITEXT_EMPTY_INDEX) {
            break;
        }
        const items = valueItemsByIndex[valueIndex] || [];
        for (const item of items) {
            if (selectedItems.has(item)) {
                return true;
            }
        }
    }
    return false;
}

function rowHasAllMultitextItems(data, start, capacity, valueItemsByIndex, selectedItems) {
    const remaining = selectedItems.slice();
    for (let offset = 0; offset < capacity; offset++) {
        const valueIndex = data[start + offset];
        if (valueIndex === MULTITEXT_EMPTY_INDEX) {
            break;
        }
        const items = valueItemsByIndex[valueIndex] || [];
        for (const item of items) {
            const index = remaining.indexOf(item);
            if (index !== -1) {
                remaining.splice(index, 1);
                if (remaining.length === 0) {
                    return true;
                }
            }
        }
    }
    return false;
}

function createMultitextCategoryPredicate(column, category) {
    const data = column.data;
    const capacity = getMultitextCapacity(column);
    const valueItemsByIndex = column.values.map((_, index) =>
        getMultitextValueItems(column, index),
    );

    if (typeof category === "string") {
        const exactIndex = column.values.indexOf(category);
        if (exactIndex !== -1) {
            return (rowIndex) => {
                const start = rowIndex * capacity;
                for (let offset = 0; offset < capacity; offset++) {
                    if (data[start + offset] === exactIndex) {
                        return true;
                    }
                }
                return false;
            };
        }

        const selectedItems = new Set([category]);
        return (rowIndex) =>
            rowHasAnyMultitextItem(
                data,
                rowIndex * capacity,
                capacity,
                valueItemsByIndex,
                selectedItems,
            );
    }

    const selectedItems = category.filter(
        (item) => typeof item === "string" && item.length > 0,
    );

    if (selectedItems.length === 0) {
        return () => false;
    }

    if (category.operand !== "and") {
        const selectedItemSet = new Set(selectedItems);
        return (rowIndex) =>
            rowHasAnyMultitextItem(
                data,
                rowIndex * capacity,
                capacity,
                valueItemsByIndex,
                selectedItemSet,
            );
    }

    if (new Set(selectedItems).size === 1) {
        const selectedItemSet = new Set([selectedItems[0]]);
        return (rowIndex) =>
            rowHasAnyMultitextItem(
                data,
                rowIndex * capacity,
                capacity,
                valueItemsByIndex,
                selectedItemSet,
            );
    }

    return (rowIndex) =>
        rowHasAllMultitextItems(
            data,
            rowIndex * capacity,
            capacity,
            valueItemsByIndex,
            selectedItems,
        );
}

class CategoryDimension extends Dimension {
    constructor(parent) {
        super(parent);
        this.worker = new Worker(
            new URL("./catWorker.js?v=11", import.meta.url),
        );
    }

    filterCategoriesPredicate(args, columns) {
        const category = args;
        const col = this.parent.columnIndex[columns[0]];
        const data = col.data;
        const vals = col.values;

        if (col.datatype === "multitext") {
            const predicate = createMultitextCategoryPredicate(col, category);
            return this.filterPredicate({ predicate });
        }

        if (typeof category === "string") {
            const ind = vals.indexOf(category);
            const predicate = (i) => data[i] === ind;
            return this.filterPredicate({ predicate });
        }

        const cats = new Set();
        for (const cat of category) {
            cats.add(vals.indexOf(cat));
        }

        const catsArr = Array.from(cats);
        // special cases for 1 or 2 categories; marginally faster than Set.has
        if (catsArr.length === 1) {
            const val = catsArr[0];
            const predicate = (i) => data[i] === val;
            return this.filterPredicate({ predicate });
        }
        if (catsArr.length === 2) {
            const val1 = catsArr[0];
            const val2 = catsArr[1];
            const predicate = (i) => data[i] === val1 || data[i] === val2;
            return this.filterPredicate({ predicate });
        }
        //this appears to be taking about twice as long as non-predicate version??
        const predicate = (i) => cats.has(data[i]);
        // const predicate = i => catsArr.includes(data[i]);
        return this.filterPredicate({ predicate });
    }
    filterCategories(args, columns) {
        if (window.predicateTest)
            return this.filterCategoriesPredicate(args, columns);
        const category = args;
        const parent = this.parent;
        const col = parent.columnIndex[columns[0]];
        const data = col.data;
        const filter = parent.filterArray;
        const localFilter = this.filterArray;

        const vals = col.values;
        const cats = new Set();
        const len = this.parent.size;

        if (col.datatype === "multitext") {
            return this.filterPredicate({
                predicate: createMultitextCategoryPredicate(col, category),
            });
        }

        if (typeof category === "string") {
            const ind = vals.indexOf(category);

            for (let i = 0; i < len; i++) {
                if (data[i] === ind) {
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
        } else {
            for (const cat of category) {
                cats.add(vals.indexOf(cat));
            }

            for (let i = 0; i < len; i++) {
                if (cats.has(data[i])) {
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
    }

    filterMultipleCategories(args, columns) {
        const categories = args;
        const parent = this.parent;
        const data = [];
        const indexes = [];
        for (let i = 0; i < columns.length; i++) {
            const col = parent.columnIndex[columns[i]];
            data.push(col.data);
            indexes.push(col.values.indexOf(args[i]));
        }
        const filter = parent.filterArray;
        const localFilter = this.filterArray;
        const len = parent.size;
        for (let i = 0; i < len; i++) {
            let con = true;
            for (let n = 0; n < data.length; n++) {
                if (data[n][i] !== indexes[n]) {
                    con = false;
                    break;
                }
            }
            if (con) {
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

    getSankeyData(callback, columns, config = {}) {
        const col1 = this.parent.columnIndex[columns[0]];
        const col2 = this.parent.columnIndex[columns[1]];
        config.values = col1.values;
        config.values2 = col2.values;

        config.method = config.method || "sankey";
        if (config.method === "proportion") {
            config.cats =
                config.category == null
                    ? null
                    : this.parent.columnIndex[columns[2]].buffer;
        }
        const t = performance.now();

        const action = (e) => {
            console.log(
                `calc sankey ${col1.name} ${col2.name} : ${performance.now() - t}`,
            );
            callback(e.data);
            this.worker.removeEventListener("message", action);
        };

        this.worker.addEventListener("message", action);
        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            col1.buffer,
            config,
            col2.buffer,
        ]);
    }

    getCategories(callback, column, config = {}) {
        const col = this.parent.columnIndex[column];
        config.values = col.values;
        config.datatype = col.datatype;
        config.stringLength = col.stringLength;
        const t = performance.now();
        const action = (e) => {
            //console.log(`calc categories ${col.name} : ${performance.now()-t}`);
            callback(e.data);
            this.worker.removeEventListener("message", action);
        };

        this.worker.addEventListener("message", action);
        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            col.buffer,
            config,
        ]);
    }

    filterUnique(args, columns) {
        const column = this.parent.columnIndex[columns[0]];
        const searchVals = args.split(" ");
        const predicate = (i) => {
            const v = column.getValue(i);
            return searchVals.every((s) => v.indexOf(s) !== -1);
        }
        return this.filterPredicate({ predicate });
    }
    destroy() {
        super.destroy();
        this.worker.terminate();
    }
}

Dimension.types["category_dimension"] = CategoryDimension;

export default CategoryDimension;
