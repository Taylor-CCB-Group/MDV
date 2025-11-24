import WGLScatterPlot from "./WGLScatterPlot.js";
import BaseChart from "./BaseChart";
import { geoPath } from "d3-geo";
import { scaleLinear } from "d3-scale";
import { createEl } from "../utilities/Elements.js";
import "../datastore/DensityDimension.js";
import { loadColumnData } from "@/datastore/decorateColumnMethod";

class DensityScatterPlot extends WGLScatterPlot {
    @loadColumnData
    setParams(p) {
        // super.setParams(p);
        this.config.param = p;
        // this.afterAppCreation();
        this._updateParamScale(p);
        this.onDataFiltered();
        this.drawChart();
        super.drawChart();
    }
    _updateParamScale(p) {
        
        const mm1 = this.dataStore.getMinMaxForColumn(p[0]);
        const mm2 = this.dataStore.getMinMaxForColumn(p[1]);
        const w = mm1[1] - mm1[0];
        const h = mm2[1] - mm2[0];
        this.whRatio = w / h;
        this.catKeys = [-1, -1];
        // this.data = [null, null];
        this.orig_y_scale = this.y_scale.domain();
        this.orig_x_scale = this.x_scale.domain();
        this.or_y_scale = [0, 400 / this.whRatio];
        this.or_x_scale = [0, 400];
    }
    afterAppCreation() {
        //make sure svg is on top of scatter plot and events pass through
        this.contentDiv.prepend(this.app.div_container);
        this.svg.style("position", "absolute").style("pointer-events", "none");
        const c = this.config;
        c.contour_bandwidth = c.contour_bandwidth || 5;
        c.contour_intensity = c.contour_intensity || 0.7;
        c.contour_opacity = c.contour_opacity || 1;
        this.data = [null, null];
        const p = c.param;
        this._updateParamScale(p);
        this.app.addHandler(
            "pan_or_zoom",
            (offset, x_scale, y_scale) => {
                this._updateScale(this.app.getRange());
                this._rescaleSVG();
            },
            c.id,
        );

        return super.afterAppCreation();
    }

    getDimension() {
        return this.dataStore.getDimension("density_dimension");
    }

    onDataFiltered(dim) {
        const c = this.config;
        const p = c.param;
        super.onDataFiltered(dim);
        const vals = this.dataStore.getColumnValues(c.param[2]);
        try {
            this.catKeys = [
                vals.indexOf(c.category1),
                vals.indexOf(c.category2),
            ];
            if (this.catKeys[0] === -1 && this.catKeys[1] === -1) {
                this.data = [null, null];
                return;
            }
        } catch (e) {
            this.data = [null, null];
            return;
        }

        const config = {
            categories: this.catKeys,
            yscale: [this.orig_y_scale, this.or_y_scale],
            xscale: [this.orig_x_scale, this.or_x_scale],
            bandwidth: this.config.contour_bandwidth,
        };
        this.dim.getDensityContours(
            (data) => {
                this.data = data;
                this.drawChart();
                this._rescaleSVG();
            },
            p,
            config,
        );
    }

    setSize(x, y) {
        super.setSize(x, y);
        this._rescaleSVG();
    }

    _drawContours(color, colorScale, i) {
        const dpc = `dp-poly-${i}`;
        const c = this.config;
        this.graph_area
            .selectAll(`.${dpc}`)
            .data(this.data[i])
            .join("path")
            .attr("class", dpc)
            .attr("stroke", c.contour_fill ? "none" : color)
            .style("opacity", c.contour_opacity)
            .attr("stroke-linejoin", "round")
            .attr("d", geoPath())
            .attr("fill", (d) => {
                if (c.contour_fill) {
                    return colorScale(d.value);
                }
                return "none";
            });
    }

    removeCategory(c) {
        delete this.config[`category${c}`];
        this.data[c - 1] = null;
        this.catKeys[c - 1] = -1;
        this.graph_area.selectAll(`.dp-poly-${c - 1}`).remove();
    }

    changeContourParameter(newCat) {
        this.config.param[2] = newCat;
        this.removeCategory(1);
        this.removeCategory(2);
        this.drawChart();
    }

    drawChart(tTime = 400) {
        const c = this.config;
        for (const i in this.catKeys) {
            const ck = this.catKeys[i];
            if (ck === -1) {
                this.graph_area.selectAll(`.dp-poly-${i}`).remove();
            } else {
                const color = this.dataStore.getColumnColors(c.param[2])[ck];
                const colorScale = scaleLinear(
                    [0, 1 - c.contour_intensity],
                    ["white", color],
                );
                this._drawContours(color, colorScale, i);
            }
        }
    }

    getContextMenu() {
        const cm = super.getContextMenu();
        if (this.dataStore.regions && !this.config.viv) {
            cm.push({
                text: "Set as Default Image",
                icon: "fas fa-file-image",
                func: () => {
                    const c = this.config;
                    if (c.background_image) {
                        this.dataStore.regions.all_regions[
                            c.region
                        ].default_image = c.background_image.name;
                        this.dataStore.dirtyMetadata.add("regions");
                    }
                },
            });
        }
        return cm;
    }

    _rescaleSVG() {
        const xr = this.x_scale.domain();
        const xlen = xr[1] - xr[0];
        const jj = this.x_scale.range()[1] / this.or_x_scale[1];
        const nxs = ((this.orig_x_scale[1] - this.orig_x_scale[0]) / xlen) * jj;
        const nsx = scaleLinear().domain(xr).range(this.or_x_scale);
        const nxp = nsx(xr[0]) - nsx(this.orig_x_scale[0]) * jj;
        const yr = this.y_scale.domain();
        const ylen = yr[0] - yr[1];
        const ll = this.y_scale.range()[1] / this.or_y_scale[1];
        const nys = ((this.orig_y_scale[0] - this.orig_y_scale[1]) / ylen) * ll;
        const nsy = scaleLinear().domain(yr).range(this.or_y_scale);
        const nyp = nsy(yr[0]) - nsy(this.orig_y_scale[0]) * ll;
        const x = this.margins.left;
        const y = this.margins.top;
        this.graph_area.attr(
            "transform",
            `translate(${x - nxp},${y - nyp}) scale(${nxs},${nys})  `,
        );
    }

    getSettings() {
        const s = super.getSettings();
        const c = this.config;
        const cols = this.dataStore.getColumnList("text");
        // returning undefined e.g. for numeric 'DAPI' column...
        // changing this at least means it doesn't crash (need to review 'Contour Category' stuff though...)
        let cats = this.dataStore.getColumnValues(c.param[2]) || [];
        cats = cats.map((x) => {
            return { t: x };
        });

        cats.push({ t: "None" });

        const contourSettings = [
            {
                type: "dropdown",
                label: "Contour Parameter",
                current_value: this.dataStore.getColumnName(c.param[2]),
                values: [cols, "name", "field"],
                func: (x) => {
                    this.changeContourParameter(x);
                },
                onchange: (controls, x) => {
                    const cats = this.dataStore.getColumnValues(x);
                    const cs = cats.slice(0);
                    cs.push("None");
                    for (const n of [
                        "Contour Category 1",
                        "Contour Category 2",
                    ]) {
                        controls[n].innerHTML = "";
                        for (const c of cs) {
                            createEl(
                                "option",
                                {
                                    text: c,
                                    value: c,
                                },
                                controls[n],
                            );
                        }
                        controls[n].value = "None";
                    }
                },
            },

            {
                type: "dropdown",
                label: "Contour Category 1",
                current_value: c.category1 || "None",
                values: [cats, "t", "t"],
                func: (x) => {
                    if (x === "None") {
                        this.removeCategory(1);
                    } else {
                        c.category1 = x;
                    }
                    this.onDataFiltered();
                },
            },

            {
                type: "dropdown",
                label: "Contour Category 2",
                current_value: c.category2 || "None",
                values: [cats, "t", "t"],
                func: (x) => {
                    if (x === "None") {
                        this.removeCategory(2);
                    } else {
                        c.category2 = x;
                    }
                    this.onDataFiltered();
                },
            },
            {
                type: "slider",
                max: 25,
                min: 1,

                doc: this.__doc__,
                current_value: c.contour_bandwidth,
                label: "KDE Bandwidth",
                continuous: true,
                func: (x) => {
                    c.contour_bandwidth = x;
                    this.onDataFiltered();
                },
            },
            {
                label: "Fill Contours",
                type: "check",

                current_value: c.contour_fill,
                func: (x) => {
                    c.contour_fill = x;
                    this.drawChart();
                },
            },

            {
                type: "slider",
                max: 1,
                min: 0,

                doc: this.__doc__,
                current_value: c.contour_intensity,
                continuous: true,
                label: "Fill Intensity",
                func: (x) => {
                    c.contour_intensity = x;
                    this.drawChart();
                },
            },
            {
                type: "slider",
                max: 1,
                min: 0,

                doc: this.__doc__,
                current_value: c.contour_opacity,
                continuous: true,
                label: "Contour opacity",
                func: (x) => {
                    c.contour_opacity = x;
                    this.drawChart();
                },
            },
        ];
        s.push({
            type: "folder",
            label: "Contour Settings",
            current_value: contourSettings
        });
        return s;
    }
}

BaseChart.types["density_scatter_plot"] = {
    name: "Density Scatter Plot",
    class: DensityScatterPlot,
    params: [
        {
            type: "number",
            name: "X axis",
        },
        {
            type: "number",
            name: "Y axis",
        },
        {
            type: "text",
            name: "Category Column",
        },
    ],
};

BaseChart.types["image_scatter_plot"] = {
    name: "Centroid Plot",
    class: DensityScatterPlot,
    required: ["regions"],
    extra_controls: (ds) => {
        const vals = [];
        for (const x in ds.regions.all_regions) {
            vals.push({ name: x, value: x });
        }
        return [
            {
                type: "dropdown",
                name: "region",
                label: ds.getColumnName(ds.regions.region_field),
                values: vals,
            },
        ];
    },
    init: (config, ds, ec) => {
        const r = ds.regions;
        const sr = r.all_regions[ec.region];
        config.color_by = r.default_color;
        // trouble... when default_color is not categorical bad things happen.
        // for now, we should be able to check whether default_color is categorical and if not, we use the region field.
        // we need to do this during chart creation for it to work with saved configs.
        config.param = r.position_fields.concat([r.default_color]);
        config.background_filter = {
            column: r.region_field,
            category: ec.region,
        };
        config.color_legend = { display: false };
        config.region = ec.region;
        config.roi = sr.roi;
        config.json = sr.json;
        config.background_image = sr.images[sr.default_image];
        config.title =
            ec.region + (sr.default_image ? `-${sr.default_image}` : "");
    },
};

export default DensityScatterPlot;
