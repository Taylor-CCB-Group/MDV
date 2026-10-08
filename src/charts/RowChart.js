import { select } from "d3-selection";
import { easeLinear } from "d3-ease";
import CategoryChart from "./CategoryChart.js";
import BaseChart from "./BaseChart";
import { createEl } from "../utilities/Elements.js";
import WordCloud from "wordcloud";

function mulberry32(seed) {
    let state = seed >>> 0;
    return () => {
        state = (state + 0x6d2b79f5) | 0;
        let t = Math.imul(state ^ (state >>> 15), 1 | state);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

function monochromeInk(backgroundColor) {
    const probe = document.createElement("canvas").getContext("2d");
    if (!probe) return "#000";
    probe.fillStyle = "#000000";
    probe.fillStyle = backgroundColor || "#ffffff";
    const hex = probe.fillStyle;
    const red = Number.parseInt(hex.slice(1, 3), 16);
    const green = Number.parseInt(hex.slice(3, 5), 16);
    const blue = Number.parseInt(hex.slice(5, 7), 16);
    if ([red, green, blue].some(Number.isNaN)) return "#000";
    const luminance = (red * 299 + green * 587 + blue * 114) / 1000;
    return luminance > 160 ? "#000" : "#fff";
}

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
        this.wordcloudCanvas.addEventListener("mouseleave", () =>
            this.hideToolTip(),
        );
    }

    drawWordCloud(data) {
        const vals = this.dataStore.getColumnValues(this.config.param[0]);
        const colors = this.dataStore.getColumnColors(this.config.param[0]);
        const canvas = this.wordcloudCanvas;
        const backgroundColor =
            getComputedStyle(canvas).getPropertyValue("--main_panel_color");
        const ink = monochromeInk(backgroundColor);
        const color = (word) => {
            const filtered =
                this.filter.length > 0 && this.filter.indexOf(word) === -1;
            if (filtered) return "lightgray";
            if (this.config.black_and_white) return ink;
            return colors[vals.indexOf(word)];
        };
        const click = (item, dimension, e) =>
            this.filterCategories(item[0], e.shiftKey);
        const hover = (item, dimension, event) => {
            if (!item) {
                this.hideToolTip();
                return;
            }
            this.showToolTip(event, `${item[0]}: ${item[2]}`);
        };
        // log(count) is 0 when a category appears once, and wordcloud skips that size.
        const list = data.map((d) => [
            vals[d[1]],
            d[0] > 0 ? Math.log(d[0] + 1) : 0.25,
            d[0],
        ]);
        const maxVal = list.reduce(
            (max, entry) =>
                Number.isFinite(entry[1]) ? Math.max(max, entry[1]) : max,
            0,
        );
        const w = (canvas.width = this.contentDiv.clientWidth);
        const h = (canvas.height = this.contentDiv.clientHeight);
        const shortest = Math.min(w, h);
        const wordSize = this.config.wordSize || 20;
        const crowd = Math.min(1, Math.sqrt(36 / Math.max(list.length, 1)));
        const largestPx = ((Math.max(shortest, 0) * wordSize) / 100) * crowd;
        const weightFactor = maxVal > 0 ? largestPx / maxVal : 1;
        const gridSize = Math.max(2, Math.round(shortest / Math.max(64, list.length)));
        const p2 = Math.PI / 2;
        canvas.style.display = "block";
        this.graph_area.style.display = "none";
        const options = {
            list,
            color,
            click,
            hover,
            weightFactor,
            gridSize,
            shrinkToFit: true,
            backgroundColor,
            minRotation: -p2,
            maxRotation: p2,
            rotationSteps: 2,
        };
        console.log("wordcloud options", options);
        this._restoreWordcloudRandom?.();
        WordCloud.stop();
        if (this._wordcloudSeed == null) {
            this._wordcloudSeed = (Math.random() * 0x100000000) >>> 0;
        }
        const nativeRandom = Math.random;
        const seeded = mulberry32(this._wordcloudSeed);
        Math.random = seeded;
        const restoreRandom = () => {
            if (Math.random === seeded) Math.random = nativeRandom;
            canvas.removeEventListener("wordcloudstop", restoreRandom);
            canvas.removeEventListener("wordcloudabort", restoreRandom);
            if (this._restoreWordcloudRandom === restoreRandom) {
                this._restoreWordcloudRandom = null;
            }
        };
        this._restoreWordcloudRandom = restoreRandom;
        canvas.addEventListener("wordcloudstop", restoreRandom);
        canvas.addEventListener("wordcloudabort", restoreRandom);
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
        if (c.wordcloud) {
            return settings
                .filter((setting) => setting.label !== "Axis controls")
                .concat(this.getWordCloudSettings());
        }
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
        const columnValues = this.dataStore.getColumnValues(c.param[0]);
        const categoryCount = Math.max(
            1,
            columnValues?.length || this.data?.length || 1,
        );

        return [
            {
                type: "slider",
                label: "Word Size",
                current_value: c.wordSize || 20,
                min: 10,
                max: 100,
                func: (x) => {
                    c.wordSize = x;
                    this.drawChart();
                },
            },
            {
                type: "spinner",
                label: "Max words",
                current_value: Math.min(c.show_limit || categoryCount, categoryCount),
                min: 1,
                max: categoryCount,
                step: 1,
                func: (x) => {
                    const next = Math.min(Math.max(1, x || 1), categoryCount);
                    c.show_limit = next;
                    this.updateData();
                    this.drawChart();
                },
            },
            {
                type: "radiobuttons",
                label: "Sort Order",
                current_value: c.sort || "size",
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
                label: "Hide zero values",
                current_value: c.filter_zeros,
                func: (x) => {
                    c.filter_zeros = x;
                    this.updateData();
                    this.drawChart();
                },
            },
            {
                type: "check",
                label: "Black and white",
                current_value: c.black_and_white,
                func: (x) => {
                    c.black_and_white = x;
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
            type: ["text", "multitext", "text16"],
            name: "Category",
        },
    ],
    init: (config, dataStore) => {
        config.wordcloud = true;
        config.wordSize = 20;
        config.sort = "size";
    },
};

export default RowChart;
