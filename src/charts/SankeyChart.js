import { select } from "d3-selection";
import { easeLinear } from "d3-ease";
import SVGChart from "./SVGChart.js";
import BaseChart from "./BaseChart";
import { sankey, sankeyJustify, sankeyLinkHorizontal } from "d3-sankey";

class SankeyChart extends SVGChart {
    constructor(dataStore, div, config) {
        const at = {
            y: {
                custom: true,
            },
            ry: {
                label: true,
            },
        };
        const l1 = dataStore.getColumnName(config.param[0]);
        const l2 = dataStore.getColumnName(config.param[1]);
        if (!config.label) {
            config.label = `${l1} vs. ${l2}`;
        }
        if (!config.axis) {
            config.axis = {};
        }
        if (!config.axis.y) {
            config.axis.y = {
                size: 25,
                textSize: 15,
                label: l1,
            };
        }
        if (!config.axis.ry) {
            config.axis.ry = {
                size: 25,
                textSize: 15,
                label: l2,
            };
        }
        super(dataStore, div, config, at);
        const c = this.config;
        c.type = "sankey_chart";
        this.dim = this.dataStore.getDimension("category_dimension");
        this.filter = [];
        this.addToolTip();
        this.onDataFiltered();
    }

    setSize(x, y) {
        super.setSize(x, y);
        this.drawChart();
    }

    remove(notify = true) {
        this.dim.destroy(notify);
        super.remove();
    }

    removeFilter() {
        this.dim.removeFilter();
        this.filter = [];
        this.nodeFilter = null;
        this.linkFilter = null;
        this.drawChart();
    }

    onDataFiltered(dim) {
        //no need to change anything

        if (this.dim === dim || this.isPinned) {
            return;
        }

        if (dim === "all_removed") {
            this.filter = [];
            this.nodeFilter = null;
            this.linkFilter = null;
            this.resetButton.style.display = "none";
        }

        this.dim.getSankeyData((data) => {
            this.sankeyData = data;
            this.drawChart();
        }, this.config.param);
    }

    drawChart() {
        this.updateAxis();
        const trans = select(this.contentDiv)
            .transition()
            .duration(400)
            .ease(easeLinear);
        const cdim = this._getContentDimensions();
        const p = this.config.param;
        const colors = [
            this.dataStore.getColumnColors(p[0]),
            this.dataStore.getColumnColors(p[1]),
        ];
        let nodePadding = 20 - this.sankeyData.minNodes;
        nodePadding = nodePadding < 1 ? 1 : nodePadding > 10 ? 10 : nodePadding;
        const san = sankey()
            .size([cdim.width, cdim.height])
            .nodeId((d) => d.id)
            .nodeWidth(20)
            .nodePadding(nodePadding)
            .nodeAlign(sankeyJustify);

        const graph = san(this.sankeyData);

        const links = this.graph_area
            .selectAll(".s-links")

            .data(
                graph.links /*,d=>{
            return d.source.id+d.target.id;
        }*/,
            )
            .join("path")
            .on("click", (e, d) => {
                this.filterMultipleCategories(d);
            })
            .on("mouseover pointermove", (e, d) => {
                const s = d.source;
                const sn = this.dataStore.getColumnValues(
                    this.config.param[s.param],
                )[s.ind];
                const t = d.target;
                const tn = this.dataStore.getColumnValues(
                    this.config.param[t.param],
                )[t.ind];
                this.showToolTip(e, `${sn}<br>${tn}<br>${d.trueValue}`);
            })
            .on("mouseleave", () => {
                this.hideToolTip();
            })
            .transition(trans)
            .attr("class", "s-links")
            .attr("d", sankeyLinkHorizontal())
            .attr("fill", "none")
            .attr("stroke", "currentColor")
            .attr("stroke-width", (d) => d.width)
            .attr("stoke-opacity", 0.5)
            .attr("opacity", (d) => {
                if (this.linkFilter) {
                    if (
                        d.source.id === this.linkFilter[0] &&
                        d.target.id === this.linkFilter[1]
                    ) {
                        return 0.8;
                    }
                    return 0.2;
                }

                if (this.nodeFilter) {
                    if (
                        d.source.id === this.nodeFilter ||
                        d.target.id === this.nodeFilter
                    ) {
                        return 0.8;
                    }
                    return 0.2;
                }

                return 0.5;
            });

        this.graph_area
            .selectAll(".s-nodes")

            .data(graph.nodes)
            .join("rect")
            .on("click", (e, d) => {
                const col = p[d.param];
                this.nodeFilter = (d.param === 0 ? "A" : "B") + d.ind;
                this.linkFilter = null;
                const values = this.dataStore.getColumnValues(col);
                this.filterCategories(col, values[d.ind]);
            })
            .on("mouseover pointermove", (e, d) => {
                const n = this.dataStore.getColumnValues(
                    this.config.param[d.param],
                )[d.ind];
                this.showToolTip(e, n);
            })
            .on("mouseleave", () => {
                this.hideToolTip();
            })
            .attr("class", "s-nodes")
            .attr("x", (d) => d.x0)
            .attr("y", (d) => d.y0)
            .attr("width", (d) => d.x1 - d.x0)
            .attr("height", (d) => d.y1 - d.y0)
            .attr("fill", (d) => {
                return colors[d.param][d.ind];
            });
        // .attr("opacity", 0.8);
    }

    getFilter() {
        const p = this.config.param;

        if (this.linkFilter) {
            const f = {};
            let i = Number.parseInt(this.linkFilter[0].substring(1));
            f[p[0]] = [this.dataStore.getColumnValues(p[0])[i]];
            i = Number.parseInt(this.linkFilter[1].substring(1));
            f[p[1]] = [this.dataStore.getColumnValues(p[1])[i]];
            return f;
        }

        if (this.nodeFilter) {
            const f = {};
            const pi = this.nodeFilter.startsWith("A") ? 0 : 1;
            const i = Number.parseInt(this.nodeFilter.substring(1));
            f[p[pi]] = [this.dataStore.getColumnValues(p[pi])[i]];
            return f;
        }
        return null;
    }

    filterMultipleCategories(d) {
        const p = this.config.param;
        const c1 = p[d.source.param];
        const v1 = this.dataStore.getColumnValues(c1)[d.source.ind];
        const c2 = p[d.target.param];
        const v2 = this.dataStore.getColumnValues(c2)[d.target.ind];
        this.linkFilter = [`A${d.source.ind}`, `B${d.target.ind}`];
        this.resetButton.style.display = "inline";
        this.drawChart();
        this.dim.filter("filterMultipleCategories", [c1, c2], [v1, v2]);
    }

    filterCategories(column, cat, append) {
        if (append) {
            this.filter.push(cat);
        } else {
            this.filter = [cat];
        }
        this.resetButton.style.display = "inline";
        this.drawChart();
        this.dim.filter("filterCategories", [column], this.filter);
    }

    pinChart() {
        this.isPinned = true;
    }

    unpinChart() {
        this.isPinned = false;
        this.onDataFiltered();
    }

    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        return settings;
    }
}

BaseChart.types["sankey_chart"] = {
    name: "Sankey Diagram",
    class: SankeyChart,
    params: [
        {
            type: "text",
            name: "First Group",
        },
        {
            type: "text",
            name: "Second Group",
        },
    ],
};

export default SankeyChart;
