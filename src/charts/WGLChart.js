import SVGChart from "./SVGChart.js";
import { createEl } from "../utilities/Elements.js";
import { hexToRGB } from "../datastore/DataStore.js";
import "d3-transition";
// import { loadColumnData } from "@/datastore/decorateColumnMethod";

class WGLChart extends SVGChart {
    constructor(dataStore, div, config, axisTypes) {
        super(dataStore, div, config, axisTypes);
        const c = this.config;
        if (!c.tooltip) {
            c.tooltip = { show: false };
        }

        c.default_color = c.default_color || "#377eb8";

        const box = this._getContentDimensions();
        this.graphDiv = createEl(
            "div",
            {
                styles: {
                    position: "absolute",
                    left: `${box.left}px`,
                    top: `${box.top}px`,
                    width: `${box.width}px`,
                    height: `${box.height}px`,
                },
            },
            this.contentDiv,
        );
        this.tooltip = createEl(
            "div",
            {
                classes: ["ciview-tooltip"],
                styles: {
                    position: "absolute",
                    zIndex: 100,
                    display: "none",
                },
            },
            this.graphDiv,
        );
    }

    onDataHighlighted(data) {
        this.app.setHighlightPoints(data.indexes);
        this.app.refresh();
    }

    afterAppCreation() {
        const c = this.config;
        this.app.addHandler("object_over", (e, index) => {
            if (c.tooltip.show) {
                this.showTooltip(e, index);
            }
        });
        this.app.addHandler("object_out", (e, index) => {
            this.tooltip.style.display = "none";
        });

        c.default_color = c.default_color || "#377eb8";
        c.on_filter = c.on_filter == null ? "hide" : c.on_filter;
        this.app.setFilterAction(c.on_filter);

        const dColor = hexToRGB(c.default_color);
        let colorFunc = () => dColor;

        if (c.color_by) {
            const conf = {
                asArray: true,
                overideValues: {
                    colorLogScale: this.config.log_color_scale,
                },
            };
            this._addTrimmedColor(c.color_by, conf);
            colorFunc = this.dataStore.getColorFunction(c.color_by, conf);
            if (!c.color_legend) {
                c.color_legend = {
                    display: true,
                };
            }

            this.setColorLegend();
        }
        this.app.addHandler("object_clicked", (index) => {
            this.dataStore.dataHighlighted([index], this);
        });

        return colorFunc;
    }

    // @loadColumnData
    setToolTipColumn(column) {
        this.config.tooltip.column = column;
    }

    showTooltip(e, index) {
        const rect = this.graphDiv.getBoundingClientRect();
        const row = this.dataStore.getRowAsObject(index);
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;
        this.tooltip.style.display = "block";
        this.tooltip.style.left = `${x}px`;
        this.tooltip.style.top = `${y}px`;
        const c = this.config.tooltip.column;
        const column = this.dataStore.columnIndex[c];
        this.tooltip.textContent = `${column.name}: ${row[c]}`;
    }

    getSetupConfig() {
        const c = this.config;
        const dColor = hexToRGB(c.default_color);
        let colorFunc = () => dColor;

        if (c.color_by) {
            colorFunc = this.dataStore.getColorFunction(c.color_by, {
                asArray: true,
            });
        }
        return {
            localFilter: this.dim.getLocalFilter(),
            globalFilter: this.dataStore.getFilter(),
            colorFunc: colorFunc,
        };
    }

    setSize(x, y) {
        super.setSize(x, y);
        const dim = this._getContentDimensions();
        if (dim.width < 1 || dim.height < 1) {
            // todo we should think about how we get here, what we do, and how we avoid it
            dim.width = 50;
            dim.height = 50;
            console.error("WGLChart: bad _getContentDimensions()");
        }
        this.app.setSize(dim.width, dim.height);
        this.graphDiv.style.left = `${dim.left}px`;
        this.graphDiv.style.top = `${dim.top}px`;
        this.updateAxis();
    }

    onDataFiltered(dim) {
        if (dim === "all_removed") {
            this.app.clearBrush();
            this.app.setFilter(false);
            this.resetButton.style.display = "none";
        }
        this.app.setHighlightPoints(null);
        this.app.refresh();
    }

    remove(notify = true) {
        this.dim.destroy(notify);
        this.app.remove();
        super.remove();
    }

    removeFilter() {
        if (this.filter === null) {
            return;
        }
        this.dim.removeFilter();
        this.app.clearBrush();
        this.app.setFilter(false);
        this.app.refresh();
    }

    // @loadColumnData
    colorByColumn(column) {
        this.config.color_by = column;
        const colorFunc = this.getColorFunction(column, true);

        const t = performance.now();
        this.app.colorPoints(colorFunc);
        console.log(`color time:${performance.now() - t}`);
        setTimeout(() => {
            this.app.refresh();
        }, 50);
    }

    colorByDefault() {
        if (this.legend) {
            this.legend.remove();
        }
        const col = hexToRGB(this.config.default_color);
        this.app.colorPoints(() => col);
        setTimeout(() => {
            this.app.refresh();
        }, 50);
    }

    changeBaseDocument(doc) {
        this.app.clearBrush();
        this.app.__doc__ = doc;
        super.changeBaseDocument(doc);
    }

    getColorOptions() {
        return {
            colorby: "all",
            has_default_color: true,
        };
    }

    pinChart() {
        this.isPinned = true;
        //store a copy of the global filter and pass it
        //to the shader
        this.pinnedFilter = this.dataStore.getFilter().slice(0);
        this.app.setGlobalFilter(this.pinnedFilter);
    }

    unpinChart() {
        this.isPinned = false;
        this.pinnedFilter = null;
        this.app.setGlobalFilter(this.dataStore.getFilter());
        this.app.refresh();
    }

    addToImage() {
        this.app.refresh();
        let url = this.app.canvas.toDataURL("image/png");
        url = url.replace(/^data:image\/png/, "data:application/octet-stream");
        const box = this._getContentDimensions();

        this.svg
            .append("image")
            .attr("xlink:href", url)
            .attr("width", box.width)
            .attr("height", box.height)
            .attr("x", box.left)
            .attr("y", box.top);
    }
    removeFromImage() {
        this.svg.select("image").remove();
    }

    getSettings(conf) {
        if (!conf) {
            conf = {
                pointMax: 50,
                pointMin: 1,
            };
        }
        const settings = super.getSettings();
        const c = this.config;
        const cols = this.dataStore.getColumnList();
        settings.splice(
            2,
            0,
            {
                type: "slider",
                max: 1,
                min: 0,
                doc: this.__doc__,
                current_value: c.opacity,
                label: "Point Opacity",
                func: (x) => {
                    c.opacity = x;
                    this.app.setPointOpacity(x);
                    this.app.refresh();
                },
            },
            {
                type: "slider",
                max: 10,
                min: 0,

                doc: this.__doc__,
                current_value: c.radius,
                label: "Point Size",
                func: (x) => {
                    c.radius = x;
                    this.app.setPointRadius(x);
                    this.app.refresh();
                    console.log(x);
                },
            },
        );

        return settings.concat([
            {
                type: "check",
                label: "Show Tooltip",
                current_value: c.tooltip.show,
                func: (x) => {
                    c.tooltip.show = x;
                    // #297 degrades active link
                    const cl = c.tooltip.column || cols[0].field;
                    this.setToolTipColumn(cl);
                },
            },
            {
                type: "column",
                label: "Tooltip value",
                // #297 degrades active link
                current_value: c.tooltip.column || cols[0].field,
                // values: [cols, "name", "field"],
                func: (x) => {
                    // checkbox won't update as a result of this...
                    //todo consider how we have more state-full settings
                    c.tooltip.show = true;
                    this.setToolTipColumn(x);
                },
            },
            {
                type: "button",
                label: "Centre Plot",
                func: (x) => {
                    this.centerGraph();
                    this.app.refresh();
                },
            },
            {
                type: "radiobuttons",
                label: "Action on Filter",
                choices: [
                    ["Hide Points", "hide"],
                    ["Gray Out Points", "grey"],
                ],
                current_value: c.on_filter,
                func: (x) => {
                    c.on_filter = x;
                    this.app.setFilterAction(x);
                    this.app.refresh();
                },
            },
        ]);
    }
}

export default WGLChart;
