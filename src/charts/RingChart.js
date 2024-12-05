import { select } from "d3-selection";
import { easeLinear } from "d3-ease";
import CategoryChart from "./CategoryChart.js";
import BaseChart from "./BaseChart";
import { pie, arc } from "d3-shape";

class RingChart extends CategoryChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config, { x: {} });
        this.config.type = "ring_chart";
        this.pie = pie()
            .sort(null)
            .value((d) => d[0]);
        this.arc = arc();
        //this will draw the chart
        this.onDataFiltered(null);
    }

    drawChart(tTime = 1000) {
        //can't get transitions to work properly
        const trans = select(this.contentDiv)
            .transition()
            .duration(tTime)
            .ease(easeLinear);

        const dim = this._getContentDimensions();
        //calculate the radius and center of the plot based on the
        //container dimensions
        const radius = Math.min(dim.width, dim.height) / 2;
        const x = dim.left + dim.width / 2;
        const y = dim.top + dim.height / 2;
        this.graph_area.attr("transform", `translate(${x},${y})`);
        this.arc.innerRadius(radius * 0.5).outerRadius(radius * 0.9);
        const fontSize = Math.floor(radius / 8);
        const vals = this.dataStore.getColumnValues(this.config.param);
        const colors = this.dataStore.getColumnColors(this.config.param);

        //only get categories with >1 count
        const data = this.pie(this.rowData.filter((x) => x[0] !== 0));
        this.graph_area
            .selectAll(".pie-arc")
            .data(data, (d, i) => {
                d.data[1];
            })
            .join("path")

            .on("click", (e, d) => {
                this.filterCategories(vals[d.data[1]], e.shiftKey);
            })

            .style("fill", (d, i) => {
                if (this.filter.length > 0) {
                    if (this.filter.indexOf(vals[d.data[1]]) === -1) {
                        return "lightgray";
                    }
                }
                return colors[d.data[1]];
            })
            .attr("class", "pie-arc")
            .attr("stroke", "white")
            .attr("d", this.arc)
            .style("stroke-width", "2px");

        //update the text
        this.graph_area
            .selectAll(".pie-text")
            .data(data, (d, i) => {
                d.data[1];
            })
            .join("text")
            .on("click", (e, d) => {
                this.filterCategories(vals[d.data[1]], e.shiftKey);
            })
            .attr("class", "pie-text")
            .text((d, i) => {
                //only show text if segment is big enough
                const angle = d.endAngle - d.startAngle;
                return angle > 0.4 ? vals[d.data[1]] : "";
            })
            .attr("x", (d) => {
                return this.arc.centroid(d)[0];
            })
            .attr("font-size", `${fontSize}px`)
            .attr("y", (d) => this.arc.centroid(d)[1])
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "central");
    }
}

BaseChart.types["ring_chart"] = {
    class: RingChart,
    name: "Pie Chart",
    params: [
        {
            type: "text",
            name: "Category",
        },
    ],
};

export default RingChart;
