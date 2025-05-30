import { createEl } from "../utilities/Elements";
import BaseChart from "./BaseChart";
import SVGChart from "./SVGChart";
import { WGL2DI } from "../webgl/WGL2DI";
import {debounce} from "../utilities/Utilities"
const color_scheme = [
    "#313695",
    "#4575B4",
    "#74ADD1",
    "#ABD9E9",
    "#E0F3F8",
    "#E0F3F8",
    "#FFFFBF",
    "#FEE090",
    "#FDAE61",
    "#F46D43",
    "#D73027",
    "#A50026",
].reverse();

class DeepToolsHeatMap extends SVGChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config, {
            x: { type: "band" },
            y: { custom: true },
        });
        this.dim = this.dataStore.getDimension("deeptools_dimension");
        const box = this._getContentDimensions();
        //create the div to house the webgl
        this.graphDiv = createEl(
            "div",
            {
                styles: {
                    position: "absolute",
                    left: `${box.left}px`,
                    top: `${box.top}px`,
                    width: `${box.width}px`,
                    height: `${box.height}px`,
                },
            },
            this.contentDiv,
        );
        const c = this.config;
        if (!c.color_scale) {
            c.color_scale = {
                range: [0, 255],
                log: false,
            };
        }
        //create the webgl
        this.app = new WGL2DI(this.graphDiv, { lock_x_axis: true });

        this.app.addHandler("object_clicked", (i, apos) => {
            const r = Math.floor(apos[1] / 20);
            const index = this.reverse_map[r];
            this.dataStore.dataHighlighted([index], this);
        });

        //load the correct data
        this.hm_config = this.dataStore.deeptools.maps[c.heatmap];
        this.x_scale.domain(this.hm_config.groups.slice(0));

        this.dataStore.loadBinaryData(this.hm_config.data).then((buff) => {
            const r = this.hm_config.rows;
            const c = this.hm_config.cols;
            this.map_index = new Uint32Array(buff, 0, r);
            const data = new Uint8Array(buff, r * 4, r * c);
            const len = data.length;
            this.dataLength = 0;
            //filter out 0 values

            const a = Date.now();
            for (let n = 0; n < data.length; n++) {
                if (data[n] !== 0) {
                    this.dataLength++;
                }
            }

            this.data = new SharedArrayBuffer(this.dataLength * 7);
            this.rowPositions = new Uint32Array(this.data, 0, this.dataLength);
            this.colPositions = new Uint16Array(
                this.data,
                this.dataLength * 4,
                this.dataLength,
            );
            this.values = new Uint8Array(
                this.data,
                this.dataLength * 6,
                this.dataLength,
            );

            this.displayData = new SharedArrayBuffer(this.dataLength * 11);
            this.x = new Float32Array(this.displayData, 0, this.dataLength);
            this.y = new Float32Array(
                this.displayData,
                this.dataLength * 4,
                this.dataLength,
            );
            this.colors = new Uint8Array(
                this.displayData,
                this.dataLength * 8,
                this.dataLength * 3,
            );

            this._debouncedColor = debounce((range) => {
                this.config.color_scale.range = [...range];
                this.updateMap(true);
            }, 200);

            //this could done in server side
            let z = 0;
            for (let n = 0; n < len; n++) {
                if (data[n] === 0) {
                    continue;
                }
                const col = n % c;
                const row = Math.floor(n / c);
                this.rowPositions[z] = this.map_index[row];
                this.colPositions[z] = col;
                this.values[z] = data[n];
                z++;
            }
            this.app.addSquares({
                x: this.x,
                y: this.y,
                colors: this.colors,
            });
            this.dim.setParameters({
                displayData: this.displayData,
                data: this.data,
                columns: c,
                rows: r,
                groups: this.hm_config.groups.length,
                length: z,
            });
            this.sortBy(c.sortBy);
        });
    }
    updateMap(colorOnly = false) {
        const c = this.hm_config.cols;
        const colorScale = this.getColorScale();
        this.dim.updateData(colorScale, colorOnly).then((data) => {
            this.total = data.total;
            this.app.squares.count = this.total;
            this.reverse_map = data.reverse_map;
            this.currentRows = data.currentRows;
            const r = this.dataStore.filterSize;
            const dim = this._getContentDimensions();
            if (!colorOnly) {
                this.app.x_scale =
                    dim.width /
                    (c * 20 + (this.hm_config.groups.length - 1) * 200);
                this.app.y_scale = dim.height / (r * 20);
                this.app.offset = [0, 0];
            }
            const gl = this.hm_config.groups.length;
            const gw = (c / gl) * 20;
            let x = 0;
            const base_color = colorScale.colors[0];
            this.app.removeAllRectangles();
            for (let n = 0; n < gl; n++) {
                this.app.addRectangle([x, 0], gw, r * 20, base_color);
                x += gw + 200;
            }
            this.app.refresh();
            this.updateAxis();
        });
    }

    sortBy(column) {
        if (column) {
            const sortCols = [{ col: column, desc: false }];
            this.dim.sort(sortCols).then(() => {
                this.updateMap();
            });
        } else {
            this.dim.setSortOrder(this.map_index);
            this.updateMap();
        }
    }

    setSize(x, y) {
        super.setSize(x, y);
        const dim = this._getContentDimensions();
        this.app.setSize(dim.width, dim.height);
        this.graphDiv.style.left = `${dim.left}px`;
        this.graphDiv.style.top = `${dim.top}px`;
        this.updateAxis();
    }

    onDataFiltered() {
        this.updateMap();
    }

    onDataHighlighted(data) {
        this.app.removeAllLines();
        const r = this.currentRows[data.indexes[0]];
        const cols = this.hm_config.cols;
        if (r) {
            this.app.addLine(
                [0, r * 20],
                [cols * 20 + 600, r * 20],
                [255, 255, 0],
            );
            this.app.addLine(
                [0, r * 20 + 20],
                [cols * 20 + 600, r * 20 + 20],
                [255, 255, 0],
            );
        }

        /* if (data.source!==this && r){
            this.app.offset=[0,-r*20 +300];
            this.app.y_scale=this.app.height/(30*20);
        */
        this.app.refresh();
    }

    getColorScale() {
        const r = this.config.color_scale.range;
        const conf = {
            datatype: "integer",
            asArray: true,
            overideValues: {
                colorLogScale: this.config.color_scale.log,
                colors: color_scheme,
                min: r[0],
                max: r[1],
                bins: 200,
            },
        };
        return {
            colors: this.dataStore.getColumnColors(
                { datatype: "integer" },
                conf,
            ),
            bins: 200,
            min: r[0],
            max: r[1],
        };
    }

    remove(notify = true) {
        this.dim.destroy(notify);
        super.remove();
    }

    getSettings() {
      
        let settings = super.getSettings();
        const c = this.config;
        const cs = c.color_scale;
        const sortCols = this.dataStore.getColumnList();
        sortCols.push({ name: "Default", field: "__default__" });
        settings = settings.concat([
            {
                type: "doubleslider",
                max: 255,
                min: 0,
                doc: this.__doc__,
                current_value: cs.range,
                label: "Color Scale",
                func: (x) => this._debouncedColor(x)
            },
            {
                type: "dropdown",
                label: "SortBy",
                current_value: c.sortBy || "__default__",
                values: [sortCols, "name", "field"],
                func: (x) => {
                    if (x === "__default__") {
                        c.sortBy = undefined;
                        this.sortBy();
                    } else {
                        c.sortBy = x;
                        this.sortBy(x);
                    }
                },
            },
        ]);
        return settings;
    }
}

BaseChart.types["deeptools_heatmap"] = {
    class: DeepToolsHeatMap,
    required: ["deeptools"],
    name: "DeepTools HeatMap",
    methodsUsingColumns: ["sortBy"],
    configEntriesUsingColumns: ["sortBy"],
    extra_controls: (dataSource) => {
        const values = Object.keys(dataSource.deeptools.maps).map((x) => ({
            name: x,
            value: x,
        }));
        return [
            {
                type: "dropdown",
                name: "map",
                label: "Map",
                values,
            },
        ];
    },

    init: (config, dataSource, ec) => {
        config.heatmap = ec["map"];
    },
};
