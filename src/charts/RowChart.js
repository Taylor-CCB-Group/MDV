import { select } from "d3-selection";
import { easeLinear } from "d3-ease";
import CategoryChart from "./CategoryChart.js";
import BaseChart from "./BaseChart";
import { createEl } from "../utilities/Elements.js";
import WordCloud from "wordcloud";

class RowChart extends CategoryChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config, { x: {} });
        //redraw the chart
        this.onDataFiltered(null);
        this.labelg = this.graph_area.append("g");
        this.wordcloudCanvas = createEl(
            "canvas",
            {
                styles: {
                    position: "absolute",
                    top: "0px",
                    left: "0px",
                    width: "100%",
                    height: "100%",
                },
            },
            this.contentDiv,
        );
    }

    drawWordCloud(data) {
        const vals = this.dataStore.getColumnValues(this.config.param[0]);
        const colors = this.dataStore.getColumnColors(this.config.param[0]);
        const color = (word) => {
            const filtered =
                this.filter.length > 0 && this.filter.indexOf(word) === -1;
            return filtered ? "lightgray" : colors[vals.indexOf(word)];
        };
        const click = (item, dimension, e) =>
            this.filterCategories(item[0], e.shiftKey);
        const list = data.map((d) => {
            return [vals[d[1]], Math.log(d[0])];
        });
        const maxVal = Math.max(...list.map((d) => d[1]));
        const canvas = this.wordcloudCanvas;
        const w = (canvas.width = this.contentDiv.clientWidth);
        const h = (canvas.height = this.contentDiv.clientHeight);
        const weightFactor = (this.config.wordSize || 20) / maxVal;
        const p2 = Math.PI / 2;
        canvas.style.display = "block";
        this.graph_area.style.display = "none";
        //won't redraw if we change the theme...
        const backgroundColor =
            getComputedStyle(canvas).getPropertyValue("--main_panel_color");
        const options = {
            list,
            color,
            click,
            weightFactor,
            backgroundColor,
            minRotation: -p2,
            maxRotation: p2,
            rotationSteps: 2,
        };
        console.log("wordcloud options", options);
        WordCloud(canvas, options);
    }

    drawChart(tTime = 400) {
        const c = this.config;
        const trans = select(this.contentDiv)
            .transition()
            .duration(tTime)
            .ease(easeLinear);
        const chartWidth = this.width - this.margins.left - this.margins.right;
        const chartHeight =
            this.height - this.margins.bottom - this.margins.top;

        const colors = this.dataStore.getColumnColors(this.config.param[0]);
        const vals = this.dataStore.getColumnValues(this.config.param[0]);

        const units = chartWidth / this.maxCount;
        let data = this.rowData;
        if (!data) {
            //seen this happen at least once when it shouldn't have, likely related to other data-loading issues
            console.error(
                ">>> No data for row chart - probably a bug, needs tracking...",
            );
        }
        let maxCount = this.maxCount;

        if (c.exclude_categories) {
            const ex = new Set();
            for (const n of c.exclude_categories) {
                ex.add(vals.indexOf(n));
            }
            data = this.rowData.filter((x, i) => !ex.has(i));
            maxCount = data.reduce((a, b) => Math.max(a[0], b[0]));
        }

        if (this.config.wordcloud) {
            this.drawWordCloud(data);
            return;
        }
        this.wordcloudCanvas.style.display = "none";
        this.graph_area.style.display = "block";

        const nBars = data.length;

        const barHeight = (chartHeight - (nBars + 1) * 3) / nBars;

        let fontSize = Math.round(barHeight);
        fontSize = fontSize > 20 ? 20 : fontSize;

        this.x_scale.domain([0, this.maxCount]);
        this.updateAxis();

        this.labelg
            .selectAll(".row-bar")
            .data(data, (d) => d[1])
            .join("rect")
            .attr("class", "row-bar")
            .on("click", (e, d) => {
                this.filterCategories(vals[d[1]], e.shiftKey);
            })
            .on("mouseover", (e, d) => {
                if (!c.show_tooltip) return;
                this.showToolTip(e, `${vals[d[1]]}: ${d[0]}`);
            })
            .on("mouseout", () => this.hideToolTip())
            .transition(trans)
            .style("fill", (d) => {
                const i = d[1];
                if (this.filter.length > 0) {
                    if (this.filter.indexOf(vals[i]) === -1) {
                        return "lightgray";
                    }
                }
                return colors[i];
            })
            .attr("class", "row-bar")
            .attr("x", 0)
            .attr("width", (d) => d[0] * units)
            .attr("y", (d, i) => (i + 1) * 3 + i * barHeight)
            .attr("height", barHeight);

        this.graph_area
            .selectAll(".row-text")

            .data(data, (d) => d[1])
            .join("text")
            .on("click", (e, d) => {
                this.filterCategories(vals[d[1]], e.shiftKey);
            })
            .on("mouseover", (e, d) => {
                if (!c.show_tooltip) return;
                this.showToolTip(e, `${vals[d[1]]}: ${d[0]}`);
            })
            .on("mouseout", () => this.hideToolTip())
            .attr("class", "row-text")
            .transition(trans)
            .text((d) => (vals[d[1]] === "" ? "none" : vals[d[1]]))
            .attr("font-size", `${fontSize}px`)
            .attr("x", 5)
            .style("fill", "currentColor")
            .attr("y", (d, i) => (i + 1) * 3 + i * barHeight + barHeight / 2)
            //.attr("text-anchor", "middle")
            .attr("dominant-baseline", "central");
    }

    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        if (c.wordcloud) return this.getWordCloudSettings();
        const max = Math.max(this.data.length || 60);

        return settings.concat([
            {
                type: "spinner",
                label: "Max Rows",
                current_value: c.show_limit || max,
                max: max,
                func: (x) => {
                    c.show_limit = x;
                    this.updateData();
                    this.drawChart();
                },
            },
            {
                type: "check",
                label: "Hide zero values",
                current_value: c.filter_zeros,
                func: (x) => {
                    c.filter_zeros = x;
                    this.updateData();
                    this.drawChart();
                },
            },
            {
                type: "radiobuttons",
                label: "Sort Order",
                current_value: c.sort || "default",
                choices: [
                    ["Default", "default"],
                    ["Size", "size"],
                    ["Name", "name"],
                ],
                func: (v) => {
                    c.sort = v;
                    this.updateData();
                    this.drawChart();
                },
            },
            {
                type: "check",
                label: "Display as WordCloud",
                current_value: c.wordcloud,
                func: (x) => {
                    c.wordcloud = x;
                    this.drawChart();
                },
            },
            {
                type: "check",
                label: "Show tooltip",
                current_value: c.show_tooltip,
                func: (x) => {
                    c.show_tooltip = x;
                },
            },
        ]);
    }

    getWordCloudSettings() {
        const c = this.config;

        return [
            {
                type: "slider",
                label: "Word Size",
                current_value: c.wordSize || 100,
                min: 10,
                max: 100,
                func: (x) => {
                    c.wordSize = x;
                    this.drawChart();
                },
            },
        ];
    }
}

BaseChart.types["row_chart"] = {
    class: RowChart,
    name: "Row Chart",
    params: [
        {
            type: ["text", "multitext", "text16"],
            name: "Category",
        },
    ],
};

BaseChart.types["wordcloud"] = {
    class: RowChart,
    name: "Word Cloud",
    params: [
        {
            type: ["text", "multitext"],
            name: "Category",
        },
    ],
    init: (config, dataStore) => {
        config.wordcloud = true;
        config.wordSize = 20;
    },
};

export default RowChart;
