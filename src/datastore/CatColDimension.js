import { local } from "d3-selection";
import Dimension from "./Dimension.js";
class CatColDimension extends Dimension {
    constructor(parent) {
        super(parent);
        this.worker = new Worker(
            new URL("./catColWorker.js?v=7", import.meta.url),
        );
    }

    filterCatCol(args, columns) {
        const p = this.parent;
        const catCol = p.columnIndex[columns[0]];
        const catData = catCol.data;
        const colData = p.columnIndex[columns[1]].data;
        const thr = args.threshold || 0;
        const cat = catCol.values.indexOf(args.cat);
        const predicate = (i) => {
            const d = colData[i];
            return catData[i] === cat && !Number.isNaN(d) && d > thr;
        };
        this.filterPredicate({ predicate }, columns);
    }

    getAverages(callback, columns, config = {}) {
        config.method = config.method || "mean";
        const t = performance.now();
        const cIndex = this.parent.columnIndex;
        this.worker.onmessage = (e) => {
            console.log(`calc averages  ${performance.now() - t}`);
            callback(e.data);
        };

        config.values = cIndex[columns[0]].values;
        config.columns = columns.slice(1);
        config.y_datatype= cIndex[columns[0]].datatype;
        const colBuffers = [];
        for (let n = 1; n < columns.length; n++) {
            colBuffers.push([
                cIndex[columns[n]].buffer,
                cIndex[columns[n]].datatype,
            ]);
        }
        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            cIndex[columns[0]].buffer,
            colBuffers,
            config,
        ]);
    }

    destroy() {
        super.destroy();
        this.worker.terminate();
    }
}

Dimension.types["catcol_dimension"] = CatColDimension;

export default CatColDimension;
