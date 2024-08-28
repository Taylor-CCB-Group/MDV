import { select } from "d3-selection";
import { easeLinear } from "d3-ease";
import BaseChart from "./BaseChart.js";
import SVGChart from "./SVGChart.js";
import { scaleSqrt } from "d3-scale";
import { schemeReds } from "d3";
import { getColorLegendCustom } from "../utilities/Color.js";

class DotPlot extends SVGChart {
    constructor(dataStore, div, config) {
        config.title = config.title || dataStore.getColumnName(config.param[0]);
        super(dataStore, div, config, {
            x: { type: "band" },
            y: { type: "band" },
            ry: {},
        });
        const p = this.config.param;
        const yLabels = [];
        const c = this.config;
        //work out color scales
        c.color_scale = c.color_scale || { log: false };
        if (!c.color_legend) {
            c.color_legend = { display: true };
        }
        if (!c.fraction_legend) {
            c.fraction_legend = { display: true };
        }
        this.fractionScale = scaleSqrt().domain([0, 100]);
        for (let x = 1; x < p.length; x++) {
            yLabels.push(this.dataStore.getColumnName(p[x]));
        }
        this.x_scale.domain(yLabels);
        this.dim = this.dataStore.getDimension("catcol_dimension");
        this.addToolTip();
        this.onDataFiltered();
        this.colorScheme = schemeReds[8];
    }

    remove(notify = true) {
        this.dim.destroy(notify);
        super.remove();
    }

    removeFilter() {
        this.dim.removeFilter();
        this.filter = [];
        this.drawChart();
    }

    setColorFunction() {
        const p = this.config.param;
        const mm = this.data.mean_range;
        const conf = {
            useValue: true,
            overideValues: {
                colorLogScale: this.config.color_scale.log,
                colors: this.colorScheme,
                max: mm[1],
                min: mm[0],
            },
        };
        this.colorFunction = this.dataStore.getColorFunction(p[1], conf);
        this.setColorLegend();
    }

    getColorLegend() {
        const cs = this.config.color_scale;
        const mm = this.data.mean_range;
        const conf = {
            overideValues: {
                max: mm[1],
                min: mm[0],
                colorLogScale: cs.log,
                colors: this.colorScheme,
            },
            name: "Mean Expression",
        };
        return this.dataStore.getColorLegend(this.config.param[1], conf);
    }

    showFractionLegend() {
        const l = this.fractionLegend;
        const c = this.config;
        if (l) {
            c.fraction_legend.position = [l.offsetLeft, l.offsetTop];
            l.remove();
        }
        if (!c.fraction_legend.display) {
            this.nodeFractionLegend = undefined;
            return;
        }
        const pos = c.fraction_legend.position || [0, 0];
        this.fractionLegend = getColorLegendCustom(this.fractionScale, {
            label: "fraction",
            type: "circle",
        });
        this.contentDiv.append(this.fractionLegend);
        this.fractionLegend.style.top = `${pos[1]}px`;
        this.fractionLegend.style.left = `${pos[0]}px`;
    }

    getConfig() {
        const config = super.getConfig();
        const l = this.linkThicknessLegend;
        if (l) {
            config.fraction_legend.position = [l.offsetLeft, l.offsetTop];
        }
        return config;
    }

    onDataFiltered(dim) {
        //no need to change anything
        if (this.dim === dim || this.isPinned) {
            return;
        }
        if (dim === "all_removed") {
            this.dim.removeFilter();
            //this.drawChart();
            this.resetButton.style.display = "none";
        }
        const config = {
            method: "averages_simple",
            threshold: this.config.threshold_value,
        };

        this.dim.getAverages(
            (data) => {
                this.data = data;

                this.setColorFunction();
                this.drawChart();
            },
            this.config.param,
            config,
        );
    }

    filterCategories(cat, col) {
        this.dim.filter("filterCatCol", [this.config.param[0], col], {
            threshold: this.config.threshold_value,
            cat: cat,
        });
        this.drawChart();
        this.resetButton.style.display = "inline";
    }

    drawChart(tTime = 400) {
        const trans = select(this.contentDiv)
            .transition()
            .duration(tTime)
            .ease(easeLinear);
        const dim = this._getContentDimensions();
        const cWidth = dim.width / (this.config.param.length - 1);
        const fa = this.dim.filterMethod;
        this.setColorFunction();
        const vals = this.dataStore.getColumnValues(this.config.param[0]);
        const data = this.data.data.filter((x) => x.count !== 0);
        this.y_scale.domain(data.map((x) => vals[x.id]));
        this.updateAxis();
        const cHeight = dim.height / data.length;
        let r = (cWidth > cHeight ? cHeight : cWidth) / 2;
        r = r > 25 ? 25 : r;
        this.fractionScale.range([0, r]);
        this.showFractionLegend();
        const cyPos = cHeight / 2;
        this.graph_area
            .selectAll(".dotplot-row")
            .data(data, (d) => d.id)
            .join(
                (enter) =>
                    enter
                        .append("g")
                        .attr("class", "dotplot-row")
                        .attr(
                            "transform",
                            (d, i) => `translate(${-dim.width},${i * cHeight})`,
                        ),
                null,
                (exit) =>
                    exit
                        .transition(trans)
                        .attr(
                            "transform",
                            (d, i) => `translate(${dim.width},${i * cHeight})`,
                        )
                        .remove(),
            )
            .call((a) =>
                a
                    .transition(trans)
                    .attr("transform", (d, i) => `translate(0,${i * cHeight})`),
            ) //use call so can chain selectAll
            .selectAll(".dotplot-circle")
            .data(
                (d) => d.values,
                (d) => d.id,
            )
            .join((enter) =>
                enter
                    .append("circle")
                    .attr("class", "dotplot-circle")
                    .attr("stroke", "black")
                    .on("click", (e, d) => {
                        this.filterCategories(vals[d.cat_id], d.id);
                    }),
            )
            .attr("stroke-width", (d) => {
                if (fa) {
                    if (
                        this.dim.filterArguments.cat === vals[d.cat_id] &&
                        this.dim.filterColumns[1] === d.id
                    ) {
                        return 4;
                    }
                }
                return 1;
            })
            .transition(trans)

            .attr("cx", (d, i) => i * cWidth + 0.5 * cWidth)
            .attr("cy", cyPos)
            .attr("r", (d) => this.fractionScale(d.frac))
            .attr("fill", (d, i) => {
                return this.colorFunction(d.mean);
            });
    }

    setSize(x, y) {
        super.setSize(x, y);
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

BaseChart.types["dot_plot"] = {
    name: "Dot Plot",
    class: DotPlot,
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

export default DotPlot;
