import RangeDimension from "./RangeDimension.js";
import Dimension from "./Dimension";
class DensityDimension extends RangeDimension {
    constructor(parent) {
        super(parent);
        this.worker.terminate();
        this.worker = new Worker(
            new URL("./densityWorker.js?v=1", import.meta.url),
            {
                //fix 'Uncaught SyntaxError: Cannot use import statement outside a module' with vite
                type: "module",
            },
        );
    }

    getDensityContours(callback, columns, config = {}) {
        const t = performance.now();
        const cIndex = this.parent.columnIndex;
        this.worker.onmessage = (e) => {
            console.log(`calc density contours  ${performance.now() - t}`);
            callback(e.data);
        };

        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            [cIndex[columns[0]].buffer, cIndex[columns[0]].datatype],
            [cIndex[columns[1]].buffer, cIndex[columns[1]].datatype],
            cIndex[columns[2]].buffer,
            config,
        ]);
    }
}

Dimension.types["density_dimension"] = DensityDimension;

export default DensityDimension;
