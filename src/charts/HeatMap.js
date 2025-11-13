import { select } from "d3-selection";
import { easeLinear } from "d3-ease";
import BaseChart from "./BaseChart";
import SVGChart from "./SVGChart.js";
import { scaleLinear } from "d3-scale";
import { getHierarchicalNodes } from "../utilities/clustering.js";
import { axisLeft } from "d3";
import { loadColumnData } from "@/datastore/decorateColumnMethod";

class HeatMap extends SVGChart {
    constructor(dataStore, div, config) {
        config.title = config.title || dataStore.getColumnName(config.param[0]);
        super(dataStore, div, config, {
            x: { type: "band" },
            y: { type: "band" },
            tx: { custom: true, size: 45 },
        });

        const p = this.config.param;
        //! we used this.fieldNames instead of `param[1, ...]`
        // to avoid mutating config - which is observable - in setFields - which may be in an action
        // also note that a lot of this was actually because of the way active links worked during development
        // this is less needed now (although we still have the reaction -> mutation issue)
        this.fieldNames = p.slice(1);
        const vals = this.dataStore.getColumnValues(p[0]);
        this.x_scale.domain(vals.slice(0));
        const c = this.config;
        c.method = c.method || "mean";
        //work out color scales
        c.color_scale = c.color_scale || { log: false };
        if (!c.color_legend) {
            c.color_legend = { display: true };
        }
        this.setColorFunction();
        this.cat_scale = scaleLinear().range([0, 40]);
        this.cat_axis_svg = this.svg.append("g");
        this.cat_axis_call = axisLeft(this.cat_scale).ticks(3);
        this.dim = this.dataStore.getDimension("catcol_dimension");
        if (c.cluster_rows) {
            this.margins.right = 70;
        }
        this.addToolTip();
        this.setFields(p.slice(1));
    }
    @loadColumnData
    setParams(params) {
        this.config.param = params;
        this.setFields(params.slice(1));
        this.config.title = this.dataStore.getColumnName(params[0]);
        this.x_scale.domain(this.dataStore.getColumnValues(params[0]));
        this.setColorFunction();
        this.onDataFiltered();
    }
    // @loadColumnData
    setFields(p) {
        this.fieldNames = p;
        const yLabels = this.fieldNames.map(f => this.dataStore.getColumnName(f));
        this.y_scale.domain(yLabels);
        // this.config.param = [this.config.param[0]].concat(p);
        this.onDataFiltered();
    }

    remove(notify = true) {
        this.dim.destroy(notify);
        super.remove();
    }

    getChartData() {
        //get the data
        const d = this.data.transpose;
        const arr = new Array(d.length + 1);
        const c_name = this.dataStore.getColumnName(this.config.param[0]);
        //the order is the same as in params
        const colnames = this.config.param
            .slice(1)
            .map((x) => this.dataStore.getColumnName(x));
        arr[0] = [c_name].concat(colnames).join("\t");
        const vals = this.dataStore.getColumnValues(this.config.param[0]);
        for (let n = 0; n < d.length; n++) {
            arr[n + 1] = [vals[d[n]._id]].concat(d[n]).join("\t");
        }
        return new Blob([arr.join("\n")], { type: "text/plain" });
    }

    setColorFunction() {
        const p = this.config.param;
        const conf = {
            useValue: true,
            overideValues: {
                colorLogScale: this.config.color_scale.log,
                colors: [
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
                ],
                max: 1,
                min: 0,
            },
        };
        this.colorFunction = this.dataStore.getColorFunction(p[1], conf);
        this.setColorLegend();
    }

    getColorLegend() {
        const cs = this.config.color_scale;
        const conf = {
            overideValues: {
                max: 1,
                min: 0,
                colorLogScale: cs.log,
                colors: [
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
                ],
            },
            name: "Scale",
        };
        return this.dataStore.getColorLegend(this.config.param[1], conf);
    }

    onDataFiltered(dim) {
        //no need to change anything
        if (this.dim === dim || this.isPinned) {
            return;
        }
        if (dim === "all_removed") {
            this.filter = [];
            this.resetButton.style.display = "none";
        }
        const config = { method: this.config.method, scaleVals: [] };
        // const p = this.config.param;
        const p = [this.config.param[0], ...this.fieldNames];
        const q =
            this.config.color_scale.trim === "none"
                ? null
                : this.config.color_scale.trim;
        for (let x = 1; x < p.length; x++) {
            const col = this.dataStore.columnIndex[p[x]];
            const [min, max] = this.dataStore.getMinMaxForColumn(col);
            config.scaleVals.push([
                q ? col.quantiles[q][0] : min,
                q ? col.quantiles[q][1] : max,
            ]);
        }

        this.dim.getAverages(
            (data) => {
                this.data = data;
                try {
                    this.clusterRows();
                    this.clusterColumns();
                    this.drawChart();
                } catch (e) {
                    console.error("Error drawing heatmap:", e);
                }
            },
            p,
            config,
        );
    }

    clusterRows() {
        if (this.config.cluster_rows) {
            this.rowClusterNodes = getHierarchicalNodes(this.data.averages);
            this.data.averages.sort((x, y) => x._order - y._order);
            this.margins.right = 70;
            this.y_scale.domain(
                this.rowClusterNodes.order.map((x) =>
                    this.dataStore.getColumnName(this.fieldNames[x]),
                ),
            );
        } else {
            this.rowClusterNodes = undefined;
            this.y_scale.domain(
                this.fieldNames
                    .map((x) => this.dataStore.getColumnName(x)),
            );
            this.margins.right = 10;
            this.data.averages.sort((x, y) => x._id - y._id);
            this.ry_axis_svg.selectAll("*").remove();
        }
    }

    clusterColumns() {
        const vals = this.dataStore.getColumnValues(this.config.param[0]);
        if (this.config.cluster_columns) {
            this.columnClusterNodes = getHierarchicalNodes(this.data.transpose);

            this.data.transpose.sort((x, y) => x._order - y._order);
            this.x_scale.domain(
                this.columnClusterNodes.order.map((x) => vals[x]),
            );
            this.setAxisSize("tx", 85);
        } else {
            this.columnClusterNodes = undefined;
            this.x_scale.domain(this.data.transpose.map((x) => vals[x._id]));
            this.setAxisSize("tx", 45);
            this.data.transpose.sort((x, y) => x._id - y._id);
            this.x_scale.domain(this.data.transpose.map((x) => vals[x._id]));
            this.tx_axis_svg.selectAll("path").remove();
        }
    }

    drawChart(tTime = 400) {
        const trans = select(this.contentDiv)
            .transition()
            .duration(tTime)
            .ease(easeLinear);
        const vals = this.dataStore.getColumnValues(this.config.param[0]);
        const colors = this.dataStore.getColumnColors(this.config.param[0]);
        const dim = this._getContentDimensions();
        const recWidth = dim.width / this.data.transpose.length;
        const catRecWidth = recWidth / 2;
        const cMax = Math.max(...this.data.catTotals);
        const ytHeight = 40;
        this.cat_scale.domain([cMax, 0]);
        const cUnits = cMax === 0 ? 0 : 40 / cMax;
        const cTotals = this.data.transpose.map((x) => [
            this.data.catTotals[x._id],
            x._id,
        ]);
        const recHeight = dim.height / (this.fieldNames.length);
        this.graph_area
            .selectAll(".heatmap-row")
            .data(this.data.averages)
            .join("g")
            .attr("transform", (d, i) => `translate(0,${i * recHeight})`)
            .attr("class", "heatmap-row")
            .selectAll(".heatmap-rect")
            .data((d) =>
                d.map((x, i) => {
                    const row_id = d._id;
                    const col_id = this.data.transpose[i]._id;
                    return {
                        av: this.data.transpose[i][row_id],
                        row_id: row_id,
                        col_id: col_id,
                        miss: false,
                    };
                }),
            )
            .join("rect")
            .attr("class", "heatmap-rect")
            .attr("x", (d, i) => i * recWidth)
            .attr("height", recHeight)
            .attr("width", recWidth)
            .on("click", (e, d) => {
                const row = this.fieldNames[d.row_id];
                this._callListeners("cell_clicked", {
                    row: row,
                    col: vals[d.col_id],
                    val: d.av,
                });
            })
            .on("mouseover pointermove", (e, d) => {
                const row =
                    this.dataStore.columnIndex[this.fieldNames[d.row_id]]
                        .name;

                this.showToolTip(
                    e,
                    `${vals[d.col_id]}<br>${row}<br>${d.av.toPrecision(2)}`,
                );
            })
            .on("mouseleave", () => {
                this.hideToolTip();
            })
            .transition(trans)
            .attr("fill", (d, i) => {
                if (d.miss) {
                    return "#DCDCDC";
                }
                return this.colorFunction(d.av);
            });

        const trans2 = select(this.contentDiv)
            .transition()
            .duration(tTime)
            .ease(easeLinear);

        this.tx_axis_svg
            .selectAll(".heatmap-catsize")
            .data(cTotals, (d, i) => {
                return d[1];
            })
            .join(
                (enter) =>
                    enter
                        .append("rect")
                        .attr("class", "heatmap-catsize")
                        .attr("height", 0)
                        .attr("y", 45)
                        .attr("x", (d, i) => catRecWidth / 2 + i * recWidth)
                        .attr("width", catRecWidth),
                (update) => update,
                (exit) =>
                    exit.call((exit) =>
                        exit
                            .transition(trans2)
                            .attr("height", 0)
                            .attr("y", 45)
                            .remove(),
                    ),
            )

            .transition(trans2)
            .attr("y", (d) => {
                return ytHeight - d[0] * cUnits + 5;
            })
            .attr("height", (d) => d[0] * cUnits)
            .attr("fill", (d, i) => {
                return colors[d[1]];
            })
            .attr("width", catRecWidth)
            .attr("x", (d, i) => catRecWidth / 2 + i * recWidth);

        this.cat_axis_svg.attr("transform", `translate(${dim.left},5)`);
        this.cat_axis_svg.transition().call(this.cat_axis_call);
        this.updateAxis();
        if (this.rowClusterNodes) {
            this.drawYTree(
                this.rowClusterNodes.nodes,
                this.fieldNames.length,
            );
        }
        if (this.columnClusterNodes) {
            this.drawXTree(
                this.columnClusterNodes.nodes,
                this.data.transpose.length,
                40,
                45,
            );
        }
    }

    setSize(x, y) {
        super.setSize(x, y);
        //this.updateAxis();
        this.drawChart();
    }

    getSettings() {
        const c = this.config;
        const settings = super.getSettings();
        return settings.concat([
            {
                type: "check",
                label: "Cluster Rows",
                current_value: c.cluster_rows,
                func: (v) => {
                    c.cluster_rows = v;
                    this.clusterRows();
                    this.drawChart();
                },
            },
            {
                type: "check",
                label: "Cluster Columns",
                current_value: c.cluster_columns,
                func: (v) => {
                    c.cluster_columns = v;
                    this.clusterColumns();
                    this.drawChart();
                },
            },

            {
                type: "radiobuttons",
                label: "Averaging Method",
                current_value: c.method,
                choices: [
                    ["Mean", "mean"],
                    ["Median", "median"],
                ],

                func: (v) => {
                    c.method = v;
                    this.onDataFiltered();
                },
            },
            {
                type: "radiobuttons",
                label: "Trim to Percentile",
                current_value: c.color_scale.trim || "none",
                choices: [
                    ["No Trim", "none"],
                    ["0.001", "0.001"],
                    ["0.01", "0.01"],
                    ["0.05", "0.05"],
                ],
                func: (v) => {
                    c.color_scale.trim = v;
                    this.onDataFiltered();
                },
            },

            {
                label: "Show Color Legend",
                type: "check",

                current_value: c.color_legend ? c.color_legend.display : true,
                func: (x) => {
                    c.color_legend.display = x;
                    this.setColorLegend();
                },
            },

            {
                type: "check",
                label: "Log Color Scale",
                current_value: c.color_scale.log,
                func: (v) => {
                    c.color_scale.log = v;
                    this.setColorFunction();
                    this.drawChart();
                },
            },
        ]);
    }
}

BaseChart.types["heat_map"] = {
    name: "Heat Map",
    class: HeatMap,
    // methodsUsingColumns: ["setFields"],
    params: [
        {
            type: "text",
            name: "Categories on x-axis",
        },
        {
            type: "_multi_column:number",
            name: "Fields on y axis",
        },
    ],
};

export default HeatMap;
