import type DataStore from "./DataStore";
import Dimension from "./Dimension";

/**
 * A dimension that includes values from a given Set.
 */
class ValueSetDimension extends Dimension {
    constructor(parent: DataStore) {
        super(parent);
    }
    /**
     * Filter the dimension to only include values from the given set.
     * @paraam {Set} valueSet
     * @param {string[]} columns
     */
    filterValueset(valueSet, columns) {
        if (valueSet.size === 0) return;
        // if it's strings, we should first convert valueSet to numbers representing the index of those strings in the target column
        // we'll determine that we should do this by checking if the column is a text-like column, which we do by seeing if it has a values array
        const column = columns[0];
        const parent = this.parent;
        const col = parent.columnIndex[column];
        if (col.values) {
            const newSet = new Set();
            valueSet.forEach(v => {
                const i = col.values.indexOf(v);
                if (i !== -1) {
                    newSet.add(i);
                }
            });
            valueSet = newSet;
        }
        const filter = parent.filterArray;
        const localFilter = this.filterArray;
        const len = this.parent.size;

        for (let i=0; i<len; i++) {
            // should be following similar logic to RangeDimension here... hopefully it's correct.
            if (!valueSet.has(col.data[i])) {
                if (localFilter[i] === 0) {
                    if (++filter[i] === 1) {
                        parent.filterSize--;
                    }
                }
                localFilter[i] = 1;
            } else {
                if (localFilter[i] === 1) {
                    if (--filter[i] === 0) {
                        parent.filterSize++;
                    }
                }
                localFilter[i] = 0;
            }
        }
        this.parent._callListeners("filtered", this);
    }
}
Dimension.types["valueset_dimension"] = ValueSetDimension;