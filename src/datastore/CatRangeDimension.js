import Dimension from "./Dimension.js";
import RangeDimension from "./RangeDimension.js";
class CatRangeDimension extends RangeDimension {
    constructor(column, parent) {
        super(column, parent);
        this.worker.terminate();
        this.worker = new Worker(
            new URL("./boxPlotWorker.js?v=2", import.meta.url),
        );
    }

    getKernalDensity(callback, columns, config = {}) {
        const t = performance.now();
        const cat = this.parent.columnIndex[columns[0]];
        const val = this.parent.columnIndex[columns[1]];
        const mm = this.parent.getMinMaxForColumn(columns[1]);
        //config.bandwidth=(mm[1]-mm[0])/10;
        config.type = "kernal_density";

        if (config.scaletrim) {
            config.scaletrim = val.quantiles[config.scaletrim];
        }
        config.values = cat.values;
        const action = (e) => {
            console.log(
                `calculate violin ${cat.name} : ${performance.now() - t}`,
            );
            callback(e.data);
        };
        this.worker.onmessage = action;

        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            cat.buffer,
            [val.buffer, val.datatype],
            config,
        ]);
    }

    getMultiBoxPlotData(callback, columns, config = {}) {
        const cIndex = this.parent.columnIndex;
        const t = performance.now();
        const buffs = columns.map((x) => [
            cIndex[x].buffer,
            c.Index[x].datatype,
        ]);
        config.analysis = "multi";
        config.cat = this.parent
            .getColumnValues(columns[0])
            .indexOf(config.category);

        const action = (e) => {
            console.log(
                `calculate multiboxplot cols:${columns.length} time: ${performance.now() - t}`,
            );
            callback(e.data);
        };
        this.worker.onmessage = action;
        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            buffs.slice(1),
            buffs[0],
            config,
        ]);
    }

    getBoxPlotData(callback, columns, config = {}) {
        const t = performance.now();
        const cat = this.parent.columnIndex[columns[0]];
        const val = this.parent.columnIndex[columns[1]];

        config.analysis = "boxplot";
        config.values = cat.values;

        const action = (e) => {
            console.log(
                `calculate boxplot ${cat.name} : ${performance.now() - t}`,
            );
            callback(e.data);
        };
        this.worker.onmessage = action;

        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            cat.buffer,
            [val.buffer, val.datatype],
            config,
        ]);
    }
}

Dimension.types["catrange_dimension"] = CatRangeDimension;

export default CatRangeDimension;
