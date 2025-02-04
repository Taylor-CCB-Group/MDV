import Dimension from "./Dimension.js";
export default class SortableDimension extends Dimension {
    constructor(parent) {
        super(parent);
        this.sortWorker = new Worker(
            new URL("./sortWorker.js?v=1", import.meta.url),
        );
        this.orderBuffer = new SharedArrayBuffer(this.parent.size * 4);
        this.order = new Uint32Array(this.orderBuffer);
        for (let i = 0; i < this.parent.size; i++) this.order[i] = i;
    }

    setSortOrder(newOrder) {
        this.order.set(newOrder);
    }

    /**
     * Sorts the dimension based on the columns and order
     * @param {Array} sortColumns - An array of objects with the following properties
     * @returns {Promise} - A promise that resolves to the new order
     */
    async sort(sortColumns) {
        return await new Promise((resolve, reject) => {
            const orderBuffer = this.orderBuffer;
            const cols = sortColumns.map((x) => {
                const c = this.parent.columnIndex[x.col];
                return {
                    datatype: c.datatype,
                    buffer: c.buffer,
                    values: c.values,
                    stringLength: c.stringLength,
                    desc: x.desc,
                };
            });
            this.sortWorker.onmessage = (e) => {
                resolve(e.data);
            };
            this.sortWorker.postMessage({ orderBuffer, columns: cols });
        });
    }

    destroy() {
        super.destroy();
        this.sortWorker.terminate();
    }
}

Dimension.types["sortable_dimension"] = SortableDimension;
