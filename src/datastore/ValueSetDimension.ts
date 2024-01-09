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
     * todo: same method signature as other dimensions / base...
     */
    filterValueset(values, columns) {
        const column = columns[0];
        const parent = this.parent;
        const col = parent.columnIndex[column];
        const filter = parent.filterArray;
        const localFilter = this.filterArray;
        const len = this.parent.size;

        for (let i=0; i<len; i++) {
            // should be following similar logic to RangeDimension here... hopefully it's correct.
            if (!values.has(col.data[i])) {
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