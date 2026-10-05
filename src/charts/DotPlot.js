import { select } from "d3-selection";
import { easeLinear } from "d3-ease";
import BaseChart from "./BaseChart";
import SVGChart from "./SVGChart.js";
import { scaleSqrt } from "d3-scale";
import { schemeReds } from "d3";
import { getHierarchicalNodes } from "../utilities/clustering.js";
import { loadColumnData } from "@/datastore/decorateColumnMethod";
import { loadColumn } from "@/dataloaders/DataLoaderUtil";
import { buildColorLegendSpec } from "@/react/legend/color_legend/buildColorLegendSpec";
import FractionLegend from "@/react/components/legend/FractionLegend";
import LegendWrapper from "@/react/components/legend/LegendWrapper";
import {
    canvasMeasureText,
    layoutXLabels,
    stripSharedSubgroupSuffix,
} from "./dotplot/xAxisLabels";

class DotPlot extends SVGChart {
    constructor(dataStore, div, config) {
        const hasCustomXAxis = Boolean(config.axis?.x);
        const hasCustomYAxis = Boolean(config.axis?.y);
        config.title = config.title || dataStore.getColumnName(config.param[0]);
        const cf = config;
        if (!cf.background_filter) {
            cf.background_filter = {
                column: null,
                categories: [],
            };
        } else if (cf.background_filter.column === "__none__") {
            cf.background_filter.column = null;
        } else if (!cf.background_filter.categories) {
            cf.background_filter.categories = cf.background_filter.category
                ? [cf.background_filter.category]
                : [];
        }
       
        super(dataStore, div, config, {
            x: { type: "band" },
            y: { type: "band" },
            ry: {},
        });
        //! this.config is not the object that was passed in, it's been processed by super()
        const c = this.config; 
        
        this.addToolTip();
        this.dim = this.dataStore.getDimension("catcol_dimension");
        this.colorScheme = schemeReds[8];
        console.log('setting fields for dot plot', c.param.slice(1));
        
        //this works ok because the method is decorated... which I think means if there's a query object,
        //`setFields` won't actually be called until the link is properly initialised
        // this.setParams(c.param);
        //work out color scales
        c.color_scale = c.color_scale || { log: false };
        c.color_scale.palette = c.color_scale.palette || "reds";
      
        if (!c.color_legend) {
            c.color_legend = { display: true };
        }
        if (!c.fraction_legend) {
            c.fraction_legend = { display: true };
        }
        this.fractionLegendWrapper = new LegendWrapper(this);
        c.y_axis_order = c.y_axis_order || "data";
        this.fractionScale = scaleSqrt().domain([0, 100]);
        if (!hasCustomXAxis) {
            c.axis.x.rotate_labels = true;
            this.setAxisSize("x", 60);
        }
        // Fit and truncate labels by default
        c.axis.x.fit_labels ??= true;
        c.axis.x.truncate_labels ??= true;

        // Only thin labels when user explicitly allows it
        c.axis.x.thin_labels ??= false;
        if (!hasCustomYAxis) {
            this.setAxisSize("y", 110);
        }
        this.mobxAutorun(() => {
            this._setFields(this.config.param.slice(1));
        });
        this.setBackgroundFilter();
    }

    updateColorScheme() {
        const palette = this.config.color_scale?.palette || "reds";
        if (palette === "blue_white_red") {
            this.colorScheme = [
                "#313695",
                "#4575B4",
                "#74ADD1",
                "#ABD9E9",
                "#E0F3F8",
                "#FFFFFF",
                "#FEE090",
                "#FDAE61",
                "#F46D43",
                "#D73027",
                "#A50026",
            ];
            return;
        }
        this.colorScheme = schemeReds[8];
    }
    
    // removing this decorator as it muddies the waters with setParams.
    // - both methods having the annotation didn't particularly cause an issue..
    //   as long as we had the workaround in `setParams()` to manually look for activeQuery
    //   but that was itself a hack.
    // @loadColumnData
    _setFields(fieldNames) {
        //! we don't want to mutate the config object... 
        // we want to have a special value which signifies that it should use this behaviour.
        // then when we save state, it will have the appropriate value.
        //this.config.param = [p0, ...fieldNames]; //first is the category column
        this.fieldNames = fieldNames;
        this.colorLegendRange = undefined;
        // await cm.loadColumnSetAsync(fieldNames, this.dataStore.name);
        this.x_scale.domain(fieldNames);
        this.onDataFiltered();
    }

    /**
     * Maps each field to its x-axis label, without the `(subgroup)` suffix
     * when every field shares one.
     */
    getXAxisLabels(fieldIds) {
        const names = stripSharedSubgroupSuffix(
            fieldIds.map((f) => ({
                field: f,
                name: this.dataStore.getColumnName(f) ?? f,
                inSubgroup: Boolean(this.dataStore.columnIndex[f]?.subgroup),
            })),
        );
        return new Map(fieldIds.map((f, i) => [f, names[i]]));
    }

    /**
     * Fits the x labels to the chart's current size (see `layoutXLabels`) and
     * sets the bottom margin to match.
     * @returns the labels to display, keyed by field
     */
    fitXAxisLabels(fieldIds, labels) {
        const ax = this.config.axis.x;
        const layout = layoutXLabels(
            {
                labels: fieldIds.map((f) => labels.get(f)),
                plotLeft: this.margins.left,
                plotWidth: this.width - this.margins.left - this.margins.right,
                chartHeight: this.height,
                fontSize: ax.tickfont,
                titleSize: ax.label ? (ax.textsize ?? ax.textSize ?? 13) + 6 : 0,
                truncate: ax.truncate_labels,
                thin: ax.thin_labels,
            },
            canvasMeasureText,
        );
        this.xLabelAngle = layout.angle;
        this.xLabelFontSize = layout.fontSize;
        // keep the configured height, which applies again when auto-fit is off
        const configuredSize = ax.size;
        this.setAxisSize("x", layout.margin);
        ax.size = configuredSize;
        return new Map(fieldIds.map((f, i) => [f, layout.labels[i]]));
    }

    updateAxis() {
        super.updateAxis();
        // full name on hover, as labels may be shortened
        this.x_axis_svg
            .selectAll(".tick")
            .on("mouseover pointermove", (e, id) => {
                this.showToolTip(e, this.dataStore.getColumnName(id) ?? id);
            })
            .on("mouseout", () => this.hideToolTip());
        const ax = this.config.axis.x;
        if (!ax.fit_labels || this.xLabelAngle === undefined) {
            return;
        }
        const text = this.x_axis_svg
            .selectAll(".tick text")
            .attr("font-size", this.xLabelFontSize);
        if (this.xLabelAngle === 0) {
            text.style("text-anchor", "middle")
                .attr("dx", null)
                .attr("dy", ".71em")
                .attr("transform", null);
        } else {
            text.style("text-anchor", "end")
                .attr("dx", "-.8em")
                .attr("dy", ".10em")
                .attr("transform", "rotate(-45)");
        }
    }

    @loadColumnData
    setParams(params) {
        const oldTitle = this.config.title;
        const oldCategoryName = this.dataStore.getColumnName(this.config.param[0]);
        const newCategoryName = this.dataStore.getColumnName(params[0]);

        // Preserve user-customized titles; update only when title is still auto-derived.
        if (!oldTitle || oldTitle === oldCategoryName) {
            this.config.title = newCategoryName;
        }

        this.config.param = params;
    }

    // prefering to setFields in an autorun so there's less potential for confusion
    // @loadColumnData
    // setParams(params) {
    //     this.config.param = params;
    //     this.setFields(params.slice(1));
    // }

    remove(notify = true) {
        this.fractionLegendWrapper.unmount();
        this.dim.destroy(notify);
        super.remove();
    }

    removeFilter() {
        this.dim.removeFilter();
        this.filter = [];
        this.drawChart();
    }

   
    async setBackgroundFilter(){
        const bf = this.config.background_filter;
        const categories = Array.isArray(bf?.categories)
            ? bf.categories.filter((x) => x !== undefined && x !== null)
            : [];
        if (!bf || !bf.column || bf.column === "__none__" || categories.length === 0) {
            this.dim.clearBackgroundFilter();
        }
        else{
            await loadColumn(this.dataStore.name, bf.column);
            this.dim.setBackgroundFilter(bf.column, categories);
        }
        this.onDataFiltered();
    }

    

    getColorScaleRange() {
        if (!this.data?.mean_range) {
            return undefined;
        }
        if (!this.colorLegendRange) {
            this.colorLegendRange = this.data.mean_range.slice();
        }
        return this.colorLegendRange;
    }

    setColorFunction(updateLegend = true) {
        if (!this.data?.mean_range || !this.fieldNames?.length) {
            return;
        }
        this.updateColorScheme();
        const f = this.fieldNames[0];
        const mm = this.getColorScaleRange();
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
        if (updateLegend) {
            this.setColorLegend();
        }
    }

    getColorLegendSpec() {
        const cs = this.config.color_scale;
        const mm = this.getColorScaleRange();
        const conf = {
            overideValues: {
                max: mm[1],
                min: mm[0],
                colorLogScale: cs.log,
                colors: this.colorScheme,
            },
            name: "Mean Expression",
        };
        return buildColorLegendSpec(this.dataStore, this.fieldNames[0], conf);
    }

    showFractionLegend() {
        const l = this.fractionLegend;
        const c = this.config;
        if (l) {
            c.fraction_legend.position = [l.offsetLeft, l.offsetTop];
            c.fraction_legend.pos = c.fraction_legend.position;
            this.fractionLegendWrapper.unmount();
        }
        if (!c.fraction_legend.display) {
            this.fractionLegend = undefined;
            this.fractionLegendWrapper.unmount();
            return;
        }
        const pos = c.fraction_legend.position || c.fraction_legend.pos || [0, 0];
        const tickValues = this.fractionScale
            .ticks(4)
            .filter((x) => x > 0);
        const values = tickValues.length
            ? tickValues
            : [this.data?.frac_range?.[1] || 1];
        const spec = {
            label: "fraction",
            width: 100,
            maxHeight: Math.max(80, this.contentDiv.clientHeight - pos[1] - 8),
            items: values.map((value) => ({
                key: String(value),
                label: String(value),
                radius: this.fractionScale(value),
            })),
        };
        this.fractionLegendWrapper.render(
            spec,
            { left: `${pos[0]}px`, top: `${pos[1]}px` },
            FractionLegend,
            {
                dragHandle: ".legend-drag-handle",
                resizable: true,
            },
        );
        this.fractionLegend =
            this.fractionLegendWrapper.getWrapperElement() ?? undefined;
    }

    getConfig() {
        const config = super.getConfig();
        const l = this.fractionLegend;
        if (l) {
            config.fraction_legend.position = [l.offsetLeft, l.offsetTop];
            config.fraction_legend.pos = config.fraction_legend.position;
        }
        return config;
    }
  
    onDataFiltered(dim) {
        //no need to change anything
        //(we may want to review what isPinned means in relation to dynamic column selection - would be useful to have a way of pinning that)
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
                this.clusterRows();
                this.clusterColumns();
                //perhaps normalise should be optional
                this.fractionScale.domain([0, this.data.frac_range[1]]);
                this.setColorFunction(!this.colorLegendRange || !this.legend);
                this.drawChart();
            },
            p,
            config,
        );
    }

    clusterRows() {
        const data = this.data?.data?.filter((x) => x.count !== 0) || [];
        if (this.config.cluster_rows && data.length > 1) {
            const rowVectors = data.map((row) => {
                const v = row.values.map((item) => item.mean);
                v._id = row.id;
                return v;
            });
            this.rowClusterNodes = getHierarchicalNodes(rowVectors);
            this.rowOrder = this.rowClusterNodes.order;
            this.setAxisSize("ry", 55);
            return;
        }
        this.rowClusterNodes = undefined;
        this.rowOrder = undefined;
        this.setAxisSize("ry", 25);
        this.ry_axis_svg.selectAll("*").remove();
    }

    clusterColumns() {
        const data = this.data?.data?.filter((x) => x.count !== 0) || [];
        if (!data.length) {
            this.columnClusterNodes = undefined;
            this.columnOrder = undefined;
            this.setAxisSize("tx", 25);
            this.tx_axis_svg.selectAll("path").remove();
            return;
        }

        const colIds = data[0].values.map((v) => v.id);
        if (this.config.cluster_columns && colIds.length > 1) {
            const columnVectors = colIds.map((colId, idx) => {
                const v = data.map((row) => row.values[idx].mean);
                v._id = colId;
                return v;
            });
            this.columnClusterNodes = getHierarchicalNodes(columnVectors);
            this.columnOrder = this.columnClusterNodes.order;
            this.setAxisSize("tx", 85);
            return;
        }
        this.columnClusterNodes = undefined;
        this.columnOrder = undefined;
        this.setAxisSize("tx", 25);
        this.tx_axis_svg.selectAll("path").remove();
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
        if (!this.data?.data || !this.data?.mean_range || !this.fieldNames?.length) {
            return;
        }
        const trans = select(this.contentDiv)
            .transition()
            .duration(tTime)
            .ease(easeLinear);
        let xFieldIds = this.fieldNames;
        if (this.config.cluster_columns && this.columnOrder) {
            xFieldIds = this.columnOrder;
        }
        // sets the bottom margin, so must run before reading content dimensions
        const ax = this.config.axis.x;
        let xLabels = this.getXAxisLabels(xFieldIds);
        if (ax.fit_labels) {
            xLabels = this.fitXAxisLabels(xFieldIds, xLabels);
        } else {
            this.setAxisSize("x", ax.size);
        }
        this.x_scale.domain(xFieldIds);
        this.x_axis_call.tickFormat((id) => xLabels.get(id));
        const dim = this._getContentDimensions();
        const cWidth = dim.width / (this.fieldNames.length);
        const fa = this.dim.filterMethod;
        this.setColorFunction(false);
        // is this something we need to review to make sure that changing the category column behaves as expected?
        const vals = this.dataStore.getColumnValues(this.config.param[0]);
        const data = this.data.data.filter((x) => x.count !== 0);
        let orderedData = this.config.y_axis_order === "alphabetical"
            ? data.slice(0).sort((a, b) =>
                String(vals[a.id]).localeCompare(String(vals[b.id]), undefined, {
                    sensitivity: "base",
                    numeric: true,
                }),
            )
            : data;
        if (this.config.cluster_rows && this.rowOrder) {
            const rowById = new Map(data.map((row) => [row.id, row]));
            orderedData = this.rowOrder
                .map((id) => rowById.get(id))
                .filter((row) => row !== undefined);
        }

        const orderedDataWithColumns = orderedData.map((row) => {
            if (!this.config.cluster_columns || !this.columnOrder) {
                return row;
            }
            const valueById = new Map(row.values.map((v) => [v.id, v]));
            return {
                ...row,
                values: this.columnOrder
                    .map((id) => valueById.get(id))
                    .filter((v) => v !== undefined),
            };
        });
        this.y_scale.domain(orderedData.map((x) => vals[x.id]));
        this.updateAxis();
        if (orderedDataWithColumns.length === 0) {
            this.graph_area.selectAll(".dotplot-row").remove();
            this.showFractionLegend();
            return;
        }
        const cHeight = dim.height / orderedDataWithColumns.length;
        let r = (cWidth > cHeight ? cHeight : cWidth) / 2;
        r = r > 25 ? 25 : r;
        this.fractionScale.range([0, r]);
        this.showFractionLegend();
        const cyPos = cHeight / 2;
        this.graph_area
            .selectAll(".dotplot-row")
            .data(orderedDataWithColumns, (d) => d.id)
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
                        this.showToolTip(e, `${id}<br />${category}: ${vals[d.cat_id]}<br>percentage: ${d.frac.toFixed(0)}%`);
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
            .attr("r", (d) => this.fractionScale(d.frac)) // 
            .attr("fill", (d, i) => {
                return this.colorFunction(d.mean);
            });

        if (this.rowClusterNodes) {
            this.drawYTree(this.rowClusterNodes.nodes, orderedDataWithColumns.length);
        }
        if (this.columnClusterNodes) {
            this.drawXTree(this.columnClusterNodes.nodes, xFieldIds.length, 40, 45);
        }
    }

    setSize(x, y) {
        super.setSize(x, y);
        this.drawChart();
    }

    getSettings() {
        const c = this.config;
        const settings = super.getSettings();
        const textColumns = this.dataStore.getColumnList("string");
        const textColumnsWithNone = textColumns.slice(0);
        textColumnsWithNone.push({ name: "None", field: null });
        const bfColumn = c.background_filter.column || null;
        const bgCurrentCats = c.background_filter?.categories
            || (c.background_filter?.category ? [c.background_filter.category] : []);

        return settings.concat([
            {
                type: "dropdown",
                label: "Background Filter Column",
                current_value: bfColumn,
                values: [textColumnsWithNone, "name", "field"],
                func: (x) => {
                    c.background_filter.column = x;
                    //clear irrelevant category selections when changing column
                    c.background_filter.categories = [];
                    c.background_filter.category = undefined;
                    
                    this.setBackgroundFilter();
                },
            },
            {
                type: "category_selection",
                label: "Background Filter Categories",
                current_value: bgCurrentCats,
                sourceColumn: () => {
                    const col = c.background_filter.column;
                    return col || undefined;
                },
                getCurrentValue: () => c.background_filter.categories || [],
                func: (vals) => {
                    c.background_filter.categories = vals;
                    c.background_filter.category = vals[0];
                    void this.setBackgroundFilter();
                },
            },
            {
                type: "check",
                label: "Auto-fit X axis labels",
                current_value: c.axis.x.fit_labels,
                func: (v) => {
                    c.axis.x.fit_labels = v;
                    this.drawChart();
                },
            },
            {
                type: "check",
                label: "Shorten long X axis labels",
                current_value: c.axis.x.truncate_labels,
                func: (v) => {
                    c.axis.x.truncate_labels = v;
                    this.drawChart();
                },
            },
            {
                type: "check",
                label: "Hide overlapping X axis labels",
                current_value: c.axis.x.thin_labels,
                func: (v) => {
                    c.axis.x.thin_labels = v;
                    this.drawChart();
                },
            },
            {
                type: "radiobuttons",
                label: "Y-axis order",
                current_value: c.y_axis_order || "data",
                choices: [
                    ["Data order", "data"],
                    ["Alphabetical", "alphabetical"],
                ],
                func: (v) => {
                    c.y_axis_order = v;
                    this.drawChart();
                },
            },
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
                type: "radiobuttons",
                label: "Color Palette",
                current_value: c.color_scale.palette || "reds",
                choices: [
                    ["Reds", "reds"],
                    ["Blue-White-Red", "blue_white_red"],
                ],
                func: (v) => {
                    c.color_scale.palette = v;
                    this.setColorFunction();
                    this.drawChart();
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
                label: "Show Fraction Legend",
                type: "check",
                current_value: c.fraction_legend.display,
                func: (x) => {
                    c.fraction_legend.display = x;
                    this.showFractionLegend();
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
    // methodsUsingColumns: ["setFields"], //no longer necessary, @loadColumnData decorator will handle this
    params: [
        {
            type: ["text","text16"],
            name: "Categories on y-axis",
        },
        {
            type: "_multi_column:number",
            name: "Fields on x axis",
            //maybe we could have some information here about which method this corresponds to
        },
    ],
};

export default DotPlot;
