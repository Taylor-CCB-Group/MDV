import Dimension from "./Dimension.js";
import { MULTITEXT_EMPTY_INDEX, getMultitextCapacity, getMultitextValueItems } from "@/lib/multitext";

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
    const valueItemsByIndex = column.values.map((_, index) => getMultitextValueItems(column, index));

    if (typeof category === "string") {
        const exactIndex = column.values.indexOf(category);
        const selectedItems = new Set([category]);
        if (exactIndex !== -1) {
            return (rowIndex) => {
                const start = rowIndex * capacity;
                for (let offset = 0; offset < capacity; offset++) {
                    if (data[start + offset] === exactIndex) {
                        return true;
                    }
                }
                return rowHasAnyMultitextItem(
                    data,
                    start,
                    capacity,
                    valueItemsByIndex,
                    selectedItems,
                );
            };
        }

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

    /**
     * Keeps rows whose category is `args` (a string) or one of `args` (a list).
     * Multitext columns match individual items; a list with `operand: "and"` requires all of them.
     */
    filterCategories(args, columns) {
        const col = this.parent.columnIndex[columns[0]];
        if (col.datatype === "multitext") {
            return this.commitFilter(createMultitextCategoryPredicate(col, args));
        }
        const data = col.data;
        // lookup table indexed by category code; faster than Set.has per row
        const keep = new Uint8Array(col.values.length);
        for (const cat of typeof args === "string" ? [args] : args) {
            const index = col.values.indexOf(cat);
            if (index !== -1) keep[index] = 1;
        }
        return this.commitFilter((i) => keep[data[i]] === 1);
    }

    /** @deprecated same as {@link filterCategories} */
    filterCategoriesPredicate(args, columns) {
        return this.filterCategories(args, columns);
    }

    /**
     * Keeps rows where each column `columns[n]` has the category `args[n]`.
     */
    filterMultipleCategories(args, columns) {
        const parent = this.parent;
        const data = [];
        const indexes = [];
        for (let i = 0; i < columns.length; i++) {
            const col = parent.columnIndex[columns[i]];
            data.push(col.data);
            indexes.push(col.values.indexOf(args[i]));
        }
        return this.commitFilter((i) => {
            for (let n = 0; n < data.length; n++) {
                if (data[n][i] !== indexes[n]) return false;
            }
            return true;
        });
    }

    filterCategoryPair(args, columns) {
        const [cat1, cat2] = args;
        const col1 = this.parent.columnIndex[columns[0]];
        const col2 = this.parent.columnIndex[columns[1]];

        if (!col1 || !col2) {
            return;
        }

        const predicate1 =
            col1.datatype === "multitext"
                ? createMultitextCategoryPredicate(col1, cat1)
                : (() => {
                      const index = col1.values.indexOf(cat1);
                      if (index === -1) {
                          return () => false;
                      }
                      return (rowIndex) => col1.data[rowIndex] === index;
                  })();

        const predicate2 =
            col2.datatype === "multitext"
                ? createMultitextCategoryPredicate(col2, cat2)
                : (() => {
                      const index = col2.values.indexOf(cat2);
                      if (index === -1) {
                          return () => false;
                      }
                      return (rowIndex) => col2.data[rowIndex] === index;
                  })();

        return this.filterPredicate({
            predicate: (rowIndex) => predicate1(rowIndex) && predicate2(rowIndex),
        });
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

    getCategoryHeatmap(callback, columns, config = {}) {
        const col1 = this.parent.columnIndex[columns[0]];
        const col2 = this.parent.columnIndex[columns[1]];

        config.method = "double_cat";
        config.values = col1.values;
        config.values2 = col2.values;
        config.datatype = col1.datatype;
        config.datatype2 = col2.datatype;
        config.stringLength = col1.stringLength;
        config.stringLength2 = col2.stringLength;

        const action = (e) => {
            const payload = e.data || {};
            const matrix = payload.matrix || [];
            const rowIndexes = payload.rowIndexes || [];
            const colIndexes = payload.colIndexes || [];
            const xLabels = rowIndexes.map((i) => config.values[i]);
            const yLabels = colIndexes.map((i) => config.values2[i]);
            const counts = new Array(matrix.length);
            let maxCount = 0;

            for (let x = 0; x < matrix.length; x++) {
                const row = matrix[x] || [];
                const countRow = new Array(row.length);
                for (let y = 0; y < row.length; y++) {
                    const c = row[y]?.c || 0;
                    countRow[y] = c;
                    if (c > maxCount) {
                        maxCount = c;
                    }
                }
                counts[x] = countRow;
            }

            callback({
                xLabels,
                yLabels,
                counts,
                maxCount,
                totalCells: xLabels.length * yLabels.length,
            });
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
