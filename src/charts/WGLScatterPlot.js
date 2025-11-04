import { WGL2DI } from "../webgl/WGL2DI.js";
import WGLChart from "./WGLChart.js";
import BaseChart from "./BaseChart";
import { BaseDialog } from "../utilities/Dialog.js";
import { createEl } from "../utilities/Elements.js";
import { getProjectURL } from "../dataloaders/DataLoaderUtil.ts";

class OffsetDialog extends BaseDialog {
    constructor(scatterPlot) {
        super(
            { title: "Alter Position", width: 300, doc: scatterPlot.__doc__ },
            scatterPlot,
        );
    }

    init(scatterPlot) {
        this.plot = scatterPlot;
        this.ds = this.plot.dataStore;
        createEl("div", { text: "Offset" }, this.dialog);
        this.bdiv = createEl(
            "div",
            {
                styles: {
                    display: "flex",
                    justifyContent: "center",
                    alignItems: "center",
                    position: "relative",
                    height: "60px",
                },
            },
            this.dialog,
        );
        for (const pos of [
            ["top", "up"],
            ["bottom", "down"],
            ["right", "right"],
            ["left", "left"],
        ]) {
            this._addPButton(pos);
        }
        createEl("div", { text: "Rotation" }, this.dialog);
        this.rDiv = createEl(
            "div",
            {
                styles: {
                    height: "30px",
                    position: "relative",
                },
            },
            this.dialog,
        );
        for (const d of ["left", "right"]) {
            this._addRButton(d);
        }

        const dd = createEl("div", {}, this.dialog);
        const cname = this.ds.columnIndex[this.ds.offsets.groups].name;
        createEl("div", { text: cname }, dd);
        this.groupSelect = createEl("select", {}, dd);
        const bb = createEl(
            "button",
            { classes: ["ciview-button-sm"], text: "Reset" },
            dd,
        );
        bb.addEventListener("click", (e) => {
            this.ds.resetColumnOffsets(
                this.filter,
                this.groupSelect.value,
                true,
            );
            this._showOffsets(this.groupSelect.value);
        });
        const vals = this.ds.columnIndex[this.ds.offsets.groups].values;
        for (const v of vals) {
            createEl(
                "option",
                {
                    text: v,
                    value: v,
                },
                this.groupSelect,
            );
        }
        this.filter = this.plot.config.background_filter
            ? this.plot.config.background_filter.category
            : null;

        this.groupSelect.addEventListener("change", () => {
            this._showOffsets(this.groupSelect.value);
        });
        this._addOffsetLabels();
        this._showOffsets(vals[0]);
    }

    _alterOffset(dir) {
        const ds = this.plot.dataStore;
        let x = 0;
        let y = 0;
        const am = Number.parseFloat(this.offsetAm.value);
        if (Number.isNaN(am)) {
            return;
        }
        x += dir === "right" ? am : 0;
        x -= dir === "left" ? am : 0;
        y += dir === "up" ? am : 0;
        y -= dir === "down" ? am : 0;
        const gr = this.groupSelect.value;
        ds.setColumnOffset(
            { group: gr, offsets: [x, y], filter: this.filter },
            true,
        );
        this._showOffsets(gr);
    }

    _addOffsetLabels() {
        const ld = createEl(
            "div",
            {
                styles: {
                    position: "absolute",
                    top: "20px",
                    left: "20px",
                    display: "flex",
                    justifyContent: "center",
                    gap: "5px",
                },
            },
            this.bdiv,
        );
        this.xoffset = createEl("span", {}, ld);
        this.yoffset = createEl("span", {}, ld);
        createEl("span", { text: "step:" }, ld);
        this.offsetAm = createEl(
            "input",
            {
                styles: {
                    width: "30px",
                },
                value: 1,
            },
            ld,
        );

        const rd = createEl(
            "div",
            {
                styles: {
                    position: "absolute",
                    top: "0px",
                    justifyContent: "center",
                    left: "25px",
                    display: "flex",
                    gap: "5px",
                },
            },
            this.rDiv,
        );

        this.angle = createEl("span", {}, rd);
        createEl("span", { text: "step:" }, rd);
        this.offsetAn = createEl(
            "input",
            {
                styles: {
                    width: "30px",
                },
                value: 1,
            },
            rd,
        );
    }

    _showOffsets(group) {
        let xy = [0, 0];
        let angle = 0;
        const vs = this.ds.offsets.values[this.filter || "all"];
        if (vs?.[group]) {
            xy = vs[group].offsets;
            angle = vs[group].rotation;
        }

        this.xoffset.textContent = `x:${xy[0].toFixed(1)}`;
        this.yoffset.textContent = `y:${xy[1].toFixed(1)}`;
        this.angle.textContent = `angle:${angle.toFixed(1)}`;
    }

    _addPButton(position) {
        const styles = {
            position: "absolute",
            textAlign: "center",
            cursor: "pointer",
        };
        styles[position[0]] = "0px";
        const b = createEl(
            "i",
            {
                classes: ["fas", `fa-arrow-${position[1]}`],
                styles: styles,
            },
            this.bdiv,
        );
        b.addEventListener("click", () => {
            this._alterOffset(position[1]);
        });
    }
    _addRButton(dir) {
        const styles = {
            position: "absolute",
            textAlign: "center",
        };
        styles[dir] = "5px";
        const j = createEl(
            "i",
            {
                classes: ["fas", `fa-angle-double-${dir}`],
                styles: styles,
                cursor: "pointer",
            },
            this.rDiv,
        );
        j.addEventListener("click", (e) => {
            const am = Number.parseFloat(this.offsetAn.value);
            if (Number.isNaN(am)) {
                return;
            }
            const angle = dir === "left" ? -am : am;
            this.plot.dataStore.setColumnOffset(
                {
                    group: this.groupSelect.value,
                    rotation: angle,
                    filter: this.filter,
                },
                true,
            );
            this._showOffsets(this.groupSelect.value);
        });
    }
}

class WGLScatterPlot extends WGLChart {
    constructor(dataStore, div, config) {
        const x_name = dataStore.getColumnName(config.param[0]);
        const y_name = dataStore.getColumnName(config.param[1]);
        if (!config.axis) {
            config.axis = {
                x: { size: 30, label: x_name, textsize: 13 },
                y: { size: 45, label: y_name, textsize: 13 },
            };
        }
        if (!config.title) {
            config.title = `${x_name} x ${y_name}`;
        }
        super(dataStore, div, config, { x: {}, y: {} });

        this.x = this.config.param[0];
        this.y = this.config.param[1];
        this.dim = this.getDimension();
        const bf = config.background_filter;
        if (bf) {
            this.dim.setBackgroundFilter(bf.column, bf.category);
        }
        this.minMaxX = this.dataStore.getMinMaxForColumn(this.x);
        this.minMaxY = this.dataStore.getMinMaxForColumn(this.y);
        this.type = "wgl_scatter_plot";
        const c = this.config;
        c.brush = c.brush || "poly";

        const appConf = { brush: c.brush };

        this.app = new WGL2DI(this.graphDiv, appConf);

        this.app.addHandler("zoom_stopped", (data) => this.handlePanZoom(data));
        this.app.addHandler("panning_stopped", (data) =>
            this.handlePanZoom(data),
        );

        this.centerGraph();

        const colorFunc = this.afterAppCreation();
        this.app.setBackGroundColor(c.background_color);

        this.app.addHandler("brush_stopped", (range, is_poly) => {
            this.resetButton.style.display = "inline";
            this.app.setFilter(true);
            if (!is_poly) {
                this._createFilter(range);
            } else {
                this._createPolyFilter(range);
            }
        });

        if (c.background_image) {
            c.image_opacity = c.image_opacity == null ? 1 : c.image_opacity;
            this.addBackgroundImage(c.background_image);
        }
        const cx = this.dataStore.columnIndex[this.x];
        const cy = this.dataStore.columnIndex[this.y];
        //will get some loss of precision for int32 plus no updating on data changed
        this.app.addCircles({
            x: cx.datatype === "int32" ? new Float32Array(cx.data) : cx.data,
            y: cy.datatype === "int32" ? new Float32Array(cy.data) : cy.data,
            localFilter: this.dim.getLocalFilter(),
            globalFilter: this.dataStore.getFilter(),
            colorFunc: colorFunc,
        });

        this.defaultRadius = this._calculateRadius();

        c.radius = c.radius || 2;

        c.opacity = c.opacity == null ? 0.8 : c.opacity;

        this.app.setPointRadius(this.config.radius);
        this.app.setPointOpacity(this.config.opacity);

        this.onDataFiltered();
    }
    
    drawChart() {
        const { config } = this;
        const x_name = this.dataStore.getColumnName(config.param[0]);
        const y_name = this.dataStore.getColumnName(config.param[1]);
        config.axis.x.label = x_name;
        config.axis.y.label = y_name;

        this.x = this.config.param[0];
        this.y = this.config.param[1];
        this.dim = this.getDimension();
        const bf = config.background_filter;
        if (bf) {
            this.dim.setBackgroundFilter(bf.column, bf.category);
        }
        this.minMaxX = this.dataStore.getMinMaxForColumn(this.x);
        this.minMaxY = this.dataStore.getMinMaxForColumn(this.y);
        const colorFunc = this.afterAppCreation();
        const cx = this.dataStore.columnIndex[this.x];
        const cy = this.dataStore.columnIndex[this.y];
        //will get some loss of precision for int32 plus no updating on data changed
        this.app.addCircles({
            x: cx.datatype === "int32" ? new Float32Array(cx.data) : cx.data,
            y: cy.datatype === "int32" ? new Float32Array(cy.data) : cy.data,
            localFilter: this.dim.getLocalFilter(),
            globalFilter: this.dataStore.getFilter(),
            colorFunc: colorFunc,
        });


        this.onDataFiltered();
        this.updateAxis();
        super.drawChart();
    }

    /**
     * If called more than once, will result in more than one "RangeDimension" being created
     * this probably isn't a specific problem right now, but smells fishy to me.
     */
    getDimension() {
        return this.dataStore.getDimension("range_dimension");
    }

    handlePanZoom(range) {
        if (range.imageMoved) {
            const i = range.imageMoved;
            const c = this.config;
            const name = c.background_image.name;

            i.position[1] = -(i.position[1] + i.height);
            c.background_image.height = i.height;
            c.background_image.width = i.width;
            c.background_image.position = i.position;
            //obsolete
            if (c.image_choices) {
                const inf = c.image_choices.find((x) => x[0] === name);
                inf[2] = i.width;
                inf[3] = i.height;
                inf[4] = i.position[0];
                inf[5] = i.position[1];
            }
            if (this.dataStore.regions) {
                const im =
                    this.dataStore.regions.all_regions[c.region].images[
                        c.background_image.name
                    ];
                im.position = i.position;
                im.height = i.height;
                im.width = i.width;
                this.dataStore.dirtyMetadata.add("regions");
            }
        } else {
            this._updateScale(range);
            this.updateAxis();
        }
    }

    onDataAdded(newSize) {
        const config = this.getSetupConfig();
        config.x = this.dataStore.getRawColumn(this.x);
        config.y = this.dataStore.getRawColumn(this.y);
        this.app.updateSize(newSize, config);
        super.onDataAdded(newSize);
    }

    addBackgroundImage(ic, change = false) {
        const im = new Image();
        im.crossOrigin = "anonymous";
        const c = this.config;
        c.image_opacity = c.image_opacity || 1;
        im.src = ic.url
            ? ic.url
            : getProjectURL(this.dataStore.regions.base_url) + ic.file;
        im.onload = () => {
            if (change) {
                this.app.changeImage(im, ic, 0);
            } else {
                this.app.addImage(im, ic);
                this.app.changeImageOpacity(0, c.image_opacity);
            }
            this.app.refresh();
        };
        im.onerror = () => {
            this.config.background_image = undefined;
        };
    }

    linkToOtherChart(chart) {
        this.app.addHandler(
            "pan_or_zoom",
            (offset, x_scale, y_scale) => {
                chart.app.offset = [offset[0], offset[1]];
                chart.app.x_scale = x_scale;
                chart.app.y_scale = y_scale;
                chart.app.refresh();
                chart._updateScale(this.app.getRange());
                chart.updateAxis();
            },
            this.config.id,
        );
    }

    highlightDataItem(key) {
        this.app.highlightPoint(key);
    }

    getFilter() {
        if (!this.range || this.range === true) {
            return null;
        }
        const f = {};
        f[this.config.param[0]] = [this.range.x_min, this.range.x_max];
        f[this.config.param[1]] = [this.range.y_min, this.range.y_max];
        return f;
    }

    _updateScale(range) {
        this.x_scale.domain([range.x_range[0], range.x_range[1]]);
        this.y_scale.domain([-range.y_range[0], -range.y_range[1]]);
    }

    _calculateRadius() {
        // const width = this.width?this.width:1;
        let max_x =
            this.config.max_x || this.config.max_x === 0
                ? this.config.max_x
                : this.minMaxX[1];
        let min_x =
            this.config.min_x || this.config.min_x === 0
                ? this.config.min_x
                : this.minMaxX[0];
        if (this.config.axis.x_log_scale) {
            max_x = this._getLogValue(max_x);
            min_x = this._getLogValue(min_x);
        }

        const range = max_x - min_x;
        const pt_px = range ** 1.5 / (this.dataStore.size ** 0.8 + range);
        //pt_px = pt_px<0.001?0.001:pt_px;
        return pt_px * 5;
    }

    _createFilter(range) {
        this.range = range;
        if (range == null) {
            this.dim.removeFilter();
        } else {
            let y_max = range.y_max;
            let y_min = range.y_min;
            let x_max = range.x_max;
            let x_min = range.x_min;
            if (this.config.axis.y_log_scale) {
                y_max = this._getInverseLogValue(y_max);
                y_min = this._getInverseLogValue(y_min);
            }
            if (this.config.axis.x_log_scale) {
                x_max = this._getInverseLogValue(x_max);
                x_min = this._getInverseLogValue(x_min);
            }
            this.filter = [
                [x_min, x_max],
                [y_min, y_max],
            ];
            this.dim.filter("filterSquare", this.config.param, {
                range1: this.filter[0],
                range2: this.filter[1],
            });
        }
    }

    _createPolyFilter(vs) {
        this.range = true;
        for (const pt of vs) {
            pt[1] = -pt[1];
            if (this.config.axis.x_log_scale) {
                pt[0] = this._getInverseLogValue(pt[0]);
            }
            if (this.config.axis.y_log_scale) {
                pt[1] = this._getInverseLogValue(pt[1]);
            }
        }
        this.filter = vs;
        this.dim.filter("filterPoly", this.config.param, vs);
    }

    centerGraph() {
        const c = this.config;
        const roi = c.roi || {};

        let max_x = roi.max_x || this.minMaxX[1];
        let max_y = roi.max_y || this.minMaxY[1];
        let min_x = roi.min_x || this.minMaxX[0];
        let min_y = roi.min_y || this.minMaxY[0];
        if (c.center_using_quantiles) {
            const x_q = this.dataStore.columnIndex[c.param[0]].quantiles;
            const y_q = this.dataStore.columnIndex[c.param[1]].quantiles;
            const xmm = x_q[c.center_using_quantiles];
            const ymm = y_q[c.center_using_quantiles];
        }
        if (this.config.axis.y_log_scale) {
            max_y = this._getLogValue(max_y);
            min_y = this._getLogValue(min_y);
        }
        if (this.config.axis.x_log_scale) {
            max_x = this._getLogValue(max_x);
            min_x = this._getLogValue(min_x);
        }

        const x_margin = (max_x - min_x) / 20;
        const y_margin = (max_y - min_y) / 20;
        const x_range = max_x - min_x + 2 * x_margin;
        const y_range = max_y - min_y + 2 * y_margin;

        const dim = this._getContentDimensions();

        this.app.x_scale = dim.width / x_range;
        this.app.y_scale = dim.height / y_range;
        this.app.pointScale = (this.app.x_scale + this.app.y_scale) / 2;
        this.app.offset[0] = -(min_x - x_margin);
        this.app.offset[1] = max_y + y_margin;
        this._updateScale(this.app.getRange());
        this.updateAxis();
    }

    getContextMenu() {
        const m = super.getContextMenu();
        const c = this.config;
        const d = this.dataStore;
        //can we offset/rotate the x,y values
        if (d.offsets) {
            const o = d.offsets;
            const n = d.getColumnName(o.groups);
            //check to see if plot has correct x,y and background filter for altering offsets
            if (o.param[0] === c.param[0] && o.param[1] === c.param[1]) {
                if (
                    !o.background_filter ||
                    o.background_filter === c?.background_filter?.column
                ) {
                    m.push({
                        text: `Adjust ${n} position`,
                        icon: "fas fa-arrows-alt",
                        func: () => {
                            new OffsetDialog(this);
                        },
                    });
                }
            }
        }
        return m;
    }

    _addSliderToSettings(settings, im, index) {
        settings.splice(1, 0, {
            type: "slider",
            current_value: im.opacity,
            min: 0,
            max: 1,
            func: (x) => {
                im.opacity = x;
                this.app.changeImageOpacity(index, x);
                this.app.refresh();
            },
        });
    }

    getConfig() {
        const c = super.getConfig();
        if (c.offsets) {
            c.offsets.offsets = this.app.offsets.values.offsets;
        }
        return c;
    }

    getSettings() {
        const settings = super.getSettings({ pointMax: 30, pointMin: 0 });
        const c = this.config;
        //this is legacy and probably not used anymore
        if (c.image_choices) {
            const ic = c.image_choices.slice(0);
            ic.unshift(["None", "__none__"]);
            let cv = "__none__";
            if (c.background_image) {
                cv = c.background_image.url;
            }
            settings.splice(1, 0, {
                type: "dropdown",
                current_value: cv,
                values: [ic, 0, 1],
                label: "Change Image",
                func: (x) => {
                    if (x === "__none__") {
                        c.background_image = undefined;
                        this.app.removeImages();
                        this.app.refresh();
                        return;
                    }
                    let replace = true;
                    if (!c.background_image) {
                        c.background_image = {};
                        replace = false;
                    }

                    c.background_image.url = x;
                    const t = c.title.split("-")[0];
                    const f = c.image_choices.find((i) => i[1] === x);
                    c.background_image.name = f[0];
                    if (f[2]) {
                        c.background_image.width = f[2];
                        c.background_image.height = f[3];
                        c.background_image.position = [f[4], f[5]];
                    }
                    this.setTitle(`${t}-${f[0]}`);
                    this.addBackgroundImage(c.background_image, replace);
                },
            });
        }
        //the choice of images should be in the datastore 
        const rs = this.dataStore.regions;
        if (c.region && rs && !c.viv) {
            const ic = rs.all_regions[c.region].images;
            const vals = [{"name":"None", "field":"__none__"}];
            for (const n in ic) {
                vals.push({ name: n, field: n });
            }
            const cv = c.background_image?.name || "__none__";
            settings.splice(1, 0, {
                type: "dropdown",
                current_value: cv,
                values: [vals, "name", "field"],
                label: "Change Image",
                func: (x) => {
                    if (x === "__none__") {
                        c.background_image = undefined;
                        this.app.removeImages();
                        this.app.refresh();
                        return;
                    }
                    const replace = !!c.background_image;
                    c.background_image = ic[x];
                    this.setTitle(`${c.region}-${x}`);
                    this.addBackgroundImage(c.background_image, replace);
                },
            });
        }
        if (c.background_image) {
            settings.splice(2, 0, {
                label: "Image Opacity",
                type: "slider",
                current_value: c.image_opacity,
                max: 1,
                min: 0,
                func: (x) => {
                    c.image_opacity = x;
                    this.app.changeImageOpacity(0, x);
                    this.app.refresh();
                },
            });
        }
        settings.push({
            type: "radiobuttons",
            label: "Background Color",
            choices: [
                ["none", "none"],
                ["white", "white"],
                ["light gray", "lightgray"],
                ["gray", "gray"],
                ["black", "black"],
            ],
            current_value: c.background_color || "none",
            func: (x) => {
                if (x === "none") {
                    c.background_color = undefined;
                } else {
                    c.background_color = x;
                }
                this.app.setBackGroundColor(c.background_color);
            },
        });

        return settings.concat([
            {
                type: "radiobuttons",
                label: "Brush Type",
                choices: [
                    ["Free Draw", "poly"],
                    ["Rectangle", "default"],
                ],
                current_value: c.brush,
                func: (x) => {
                    this.app.clearBrush();
                    this.app.config.brush = x;
                    c.brush = x;
                },
            },
        ]);
    }
}

BaseChart.types["wgl_scatter_plot"] = {
    name: "2D Scatter Plot",
    class: WGLScatterPlot,
    params: [
        {
            type: "number",
            name: "X axis",
        },
        {
            type: "number",
            name: "Y axis",
        },
    ],
};

export default WGLScatterPlot;
