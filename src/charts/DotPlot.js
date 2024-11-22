import { select } from "d3-selection";
import { easeLinear } from "d3-ease";
import BaseChart from "./BaseChart";
import SVGChart from "./SVGChart.js";
import { scaleSqrt } from "d3-scale";
import { schemeReds } from "d3";
import { getColorLegendCustom } from "../utilities/Color.js";
import decorateColumnMethod, { loadColumnData } from "@/datastore/decorateColumnMethod";
import { serialiseConfig, serialiseQueries } from "./chartConfigUtils";

class DotPlot extends SVGChart {
    constructor(dataStore, div, config) {
        config.title = config.title || dataStore.getColumnName(config.param[0]);
        super(dataStore, div, config, {
            x: { type: "band" },
            y: { type: "band" },
            ry: {},
        });
        //! this.config is not the object that was passed in, it's been processed by super()
        const c = this.config; 
        
        // ! we shouldn't need to call this... but "methodUsingColumns" is applied too late
        // ! and `@loadColumnData` is not not working as expected
        decorateColumnMethod('setFields', this, dataStore.name); //dataStore.name is ignored

        this.addToolTip();
        this.dim = this.dataStore.getDimension("catcol_dimension");
        this.colorScheme = schemeReds[8];
        this.setFields(c.param.slice(1));
        //work out color scales
        c.color_scale = c.color_scale || { log: false };
        if (!c.color_legend) {
            c.color_legend = { display: true };
        }
        if (!c.fraction_legend) {
            c.fraction_legend = { display: true };
        }
        this.fractionScale = scaleSqrt().domain([0, 100]);
    }

    // todo figure out why annotation version isn't fully working - need methodsUsingColumns as well???
    //if called with objects representing queries, this will cause the method to be called with column names
    //and again whenever the column data changes
    //but when we clone a chart like that, at the moment the old chart stops responding to changes
    // @loadColumnData 
    setFields(fieldNames) {
        const cm = window.mdv.chartManager;
        //! we don't want to mutate the config object... 
        // we want to have a special value which signifies that it should use this behaviour.
        // then when we save state, it will have the appropriate value.
        //this.config.param = [p0, ...fieldNames]; //first is the category column
        this.fieldNames = fieldNames;
        // await cm.loadColumnSetAsync(fieldNames, this.dataStore.name);
        const yLabels = fieldNames.map(f => this.dataStore.getColumnName(f));
        this.x_scale.domain(yLabels);
        this.onDataFiltered();
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
        const f = this.fieldNames[0];
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
        this.colorFunction = this.dataStore.getColorFunction(f, conf);
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
        return this.dataStore.getColorLegend(this.fieldNames[0], conf);
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
        // it looks as though individual charts may be on the hook to serialise their queries
        // because only they really know how they relate to order of params in the config
        const fieldQuery = serialiseQueries(this)['setFields'];
        // if we have a way of interpreting this we'd be good to go
        if (fieldQuery) {
            console.log('DotPlot fieldQuery', fieldQuery);
            config.param = [config.param[0], ...fieldQuery];
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
        const p = [this.config.param[0], ...this.fieldNames];
        this.dim.getAverages(
            (data) => {
                this.data = data;

                this.setColorFunction();
                this.drawChart();
            },
            p,
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
        const cWidth = dim.width / (this.fieldNames.length);
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
                    })
                    .on("mouseover pointermove", (e, d) => {
                        // ['id', 'total', 'count', 'frac', 'mean', 'cat_id']
                        //const tip = { category: vals[d.cat_id], value: d.id, fraction: d.frac };
                        const category = this.dataStore.columnIndex[this.config.param[0]].name;
                        const id = this.dataStore.columnIndex[d.id].name;
                        this.showToolTip(e, `${d.id}<br />${category}: ${vals[d.cat_id]}<br>percentage: ${d.frac.toFixed(0)}%`);
                    })
                    .on("mouseout", () => {
                        this.hideToolTip();
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
            {
                // perhaps the GuiType should be more aligned with params type - i.e. _multi_column:number
                // * we should then be able to expose all params in the settings in a more consistent way *
                // then we wouldn't want current_value to be this.fieldNames, but the entries in config.param
                // corresponding to fields... as defined in BaseChart.types["dot_plot"].params
                // - param is a flat array - so we'd need to mimic the behaviour of spreading values from here
                // there's a more general question of whether settings operates on the mobx mutable config object
                // ... in many cases we could avoid having a `func`...
                type: "multicolumn",
                label: "Fields on x axis",
                // this is more of a nuisance than type: "_multi_column:number"
                columnSelection: {
                    filter: ["double", "integer", "int32"]
                },
                current_value: this.fieldNames,
                func: (v) => {
                    // given that this is "multicolumn", we should be able to assume that v is an array
                    this.setFields(Array.isArray(v) ? v : [v]);
                },
            },
        ]);
    }
}

BaseChart.types["dot_plot"] = {
    name: "Dot Plot",
    class: DotPlot,
    methodsUsingColumns: ["setFields"],
    params: [
        {
            type: "text",
            name: "Categories on y-axis",
        },
        {
            type: "_multi_column:number",
            name: "Fields on x axis",
        },
    ],
};

export default DotPlot;