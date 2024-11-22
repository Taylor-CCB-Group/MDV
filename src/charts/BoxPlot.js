import ViolinPlot from "./ViolinPlot.js";
import { easeLinear } from "d3-ease";
import { select } from "d3-selection";
import BaseChart from "./BaseChart";

class BoxPlot extends ViolinPlot {
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        this.config.type = "box_plot";
    }

    onDataFiltered(dim) {
        if (this.isPinned) {
            return;
        }

        if (dim !== this.dim) {
            if (dim === "all_removed") {
                this.app.clearBrush();
                this.app.setFilter(false);
                this.resetButton.style.display = "none";
            }
            this.ticks = this.y_scale.ticks(20);
            this.dim.getBoxPlotData(
                (data) => {
                    this.data = data;
                    this.drawChart();
                },
                [this.config.param[0], this.valueField],
                { xPos: this.xPosBuff },
            );
        }
    }

    drawChart() {
        const trans = select(this.contentDiv)
            .transition()
            .duration(400)
            .ease(easeLinear);
        const catLen = this.data.length;
        const yscale = this.y_scale;
        const values = [];
        for (const i of this.data) {
            values.push(this.values[i.id]);
        }
        this.x_scale.domain(values);
        this.centerGraph();
        const cdim = this._getContentDimensions();
        const vWidth = cdim.width / catLen;
        const colors = this.dataStore.getColumnColors(this.config.param[0]);
        this.graph_area
            .selectAll(".boxplot-rect")
            .data(this.data, (d) => d.id)
            .join("rect")
            .attr("class", "boxplot-rect")
            .transition(trans)
            .attr("x", (d, i) => i * (vWidth - 3) + (i + 1) * 3)
            .attr("y", (d) => yscale(d.Q3))
            .attr("width", vWidth - 6)
            .attr("height", (d) => {
                const w = yscale(d.Q1) - yscale(d.Q3);
                return Number.isNaN(w) ? 0 : w;
            })
            .attr("fill", "none")
            .attr("stroke", (d, i) => colors[d.id])
            .attr("stroke-width", 3);

        this.graph_area
            .selectAll(".boxplot-whisker")
            .data(this.data, (d) => d.id)
            .join("path")
            .attr("class", "boxplot-whisker")
            .transition(trans)
            .attr("stroke", (d, i) => colors[d.id])
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

        this.app.refresh();
    }
}

BaseChart.types["box_plot"] = {
    class: BoxPlot,
    name: "Box Plot",
    params: [
        {
            type: "text",
            name: "Category (X axis)",
        },
        {
            type: "number",
            name: "Value (Y axis)",
        },
    ],
};

export default BoxPlot;
