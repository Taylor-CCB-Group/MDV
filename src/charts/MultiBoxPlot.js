import BaseChart from "./BaseChart";
import SVGChart from "./SVGChart";
import { defaultPalette } from "../datastore/DataStore.js";
import { select } from "d3-selection";
import { easeLinear } from "d3-ease";

class MultiBoxPlot extends SVGChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config, { x: { type: "band" }, y: {} });
        const c = this.config;
        c.type = "multi_box_plot";
        const names = [];
        this.dim = this.dataStore.getDimension("catrange_dimension");
        for (let n = 1; n < c.param.length; n++) {
            names.push(this.dataStore.getColumnName(c.param[n]));
        }
        this.x_scale.domain(names);
        this.names = names;
        this.y_scale.domain([1, 0]);
        this.addToolTip();
        this.onDataFiltered();
    }

    onDataFiltered(dim) {
        const config = {
            analysis: "multi",
            scaleVals: [],
            category: this.config.category,
        };
        const p = this.config.param;
        const q = this.config.percentile_trim;
        for (let x = 1; x < p.length; x++) {
            const col = this.dataStore.columnIndex[p[x]];
            const [min, max] = this.dataStore.getMinMaxForColumn(col);
            config.scaleVals.push([
                q ? col.quantiles[q][0] : min,
                q ? col.quantiles[q][1] : max,
            ]);
        }
        this.dim.getMultiBoxPlotData(
            (data) => {
                this.data = data;
                this.drawChart();
            },
            this.config.param,
            config,
        );
    }

    setSize(x, y) {
        super.setSize(x, y);
        this.drawChart();
    }

    remove(notify = true) {
        this.dim.destroy(notify);
        super.remove();
    }

    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        settings.push({
            type: "radiobuttons",
            label: "Trim  to Percentile",
            current_value: c.percentile_trim || "none",
            choices: [
                ["No Trim", "none"],
                ["0.001", "0.001"],
                ["0.01", "0.01"],
                ["0.05", "0.05"],
            ],
            func: (v) => {
                if (v === "none") {
                    c.percentile_trim = undefined;
                } else {
                    c.percentile_trim = v;
                }

                this.onDataFiltered();
            },
        });
        return settings;
    }

    drawChart() {
        this.y_scale.domain([this.data.max, this.data.min]);
        this.updateAxis();
        const trans = select(this.contentDiv)
            .transition()
            .duration(400)
            .ease(easeLinear);
        const yscale = this.y_scale;
        const cdim = this._getContentDimensions();
        const vWidth = cdim.width / (this.config.param.length - 1);
        const colors = defaultPalette;
        this.graph_area
            .selectAll(".boxplot-rect")
            .data(this.data)
            .join("rect")
            .attr("class", "boxplot-rect")
            .on("mouseover pointermove", (e, d) => {
                const h = `${this.names[d.id]}`;
                this.showToolTip(e, h);
            })
            .on("mouseleave", () => {
                this.hideToolTip();
            })
            .transition(trans)
            .attr("x", (d, i) => i * (vWidth - 3) + (i + 1) * 3)
            .attr("y", (d) => yscale(d.Q3))
            .attr("width", vWidth - 6 > 1 ? vWidth - 6 : 1)
            .attr("height", (d) => {
                const w = yscale(d.Q1) - yscale(d.Q3);
                return Number.isNaN(w) ? 0 : w;
            })
            .attr("fill", "none")
            .attr("stroke", (d, i) => colors[i])
            .attr("stroke-width", 3);

        this.graph_area
            .selectAll(".boxplot-whisker")
            .data(this.data)
            .join("path")
            .attr("class", "boxplot-whisker")
            .transition(trans)
            .attr("stroke", (d, i) => colors[i])
            .attr("fill", "none")
            .attr("stroke-width", 3)

            .attr("d", (d, i) => {
                if (Number.isNaN(d.med)) {
                    return "";
                }
                const es = i * (vWidth - 6) + (i + 1) * 6;
                const mm = i * (vWidth - 3) + (i + 1) * 3;
                const m = i * vWidth + 0.5 * vWidth;
                const max = yscale(d.max);
                const min = yscale(d.min);
                if (max == null || min == null) {
                    return "";
                }
                return `
             M${es} ${max} L${es + vWidth - 12} ${max}
             M${m} ${max} L${m} ${yscale(d.Q3)}
             M${mm} ${yscale(d.med)} L${mm + vWidth - 6} ${yscale(d.med)}
             M${m} ${min} L${m} ${yscale(d.Q1)}
             M${es} ${min} L${es + vWidth - 12} ${min}`;
            });
    }
}

BaseChart.types["multi_box_plot"] = {
    name: "Multi BoxPlot",
    class: MultiBoxPlot,
    allow_user_add: false,
    params: [
        {
            type: "_multi_column:number",
            name: "Fields",
        },
    ],
};

export default MultiBoxPlot;
