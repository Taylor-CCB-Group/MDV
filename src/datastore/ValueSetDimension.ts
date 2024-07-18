import type DataStore from "./DataStore";
import Dimension from "./Dimension";

/**
 * A dimension that includes values from a given Set.
 */
class ValueSetDimension extends Dimension {
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
        const predicate = i => valueSet.has(col.data[i]);
        return this.filterPredicate({predicate}, columns);
    }
}
Dimension.types["valueset_dimension"] = ValueSetDimension;