import BaseChart from "./BaseChart.js";
import { select } from "d3-selection";
import { scaleLinear, scaleBand } from "d3-scale";
import "d3-transition";
import { easeLinear } from "d3-ease";
import { axisLeft, axisBottom } from "d3-axis";
import { cluster } from "d3-hierarchy";

class SVGChart extends BaseChart {
    constructor(dataStore, div, config, axisTypes = {}) {
        super(dataStore, div, config);
        this.svg = select(this.contentDiv)
            .append("svg")
            .attr("width", this.width)
            .attr("height", this.height);

        this.config.axis = this.config.axis || {};
        const cax = this.config.axis;
        function setDefault(ob, param, val) {
            if (ob[param] == null) {
                ob[param] = val;
            }
        }

        //add default values for axis if not specified in config
        for (const at of [
            ["x", 25],
            ["y", 40],
            ["ry", 25],
            ["tx", 25],
        ]) {
            const ati = axisTypes[at[0]];
            if (ati) {
                cax[at[0]] = cax[at[0]] || {}; //textsize:13,label:"",size:ati.size||at[1]};
                const c = cax[at[0]];
                setDefault(c, "textSize", 13);
                setDefault(c, "label", "");
                setDefault(c, "size", ati.size || at[1]);
                setDefault(c, "tickfont", 10);
            }
        }
        this._setMargins();
        this.graph_area = this.svg
            .append("g")
            .attr(
                "transform",
                `translate(${this.margins.left},${this.margins.top})`,
            );
        this.ry_axis_svg = this.svg.append("g");
        this.tx_axis_svg = this.svg.append("g");
        this.x_axis_svg = this.svg.append("g");
        this.y_axis_svg = this.svg.append("g");

        // add the x Axis
        const xa = axisTypes["x"];
        if (xa) {
            if (!xa.custom) {
                this.x_scale = xa.type === "band" ? scaleBand() : scaleLinear();
                this.x_axis_call = axisBottom(this.x_scale);
            }
            this.x_axis_label = this.x_axis_svg
                .append("text")
                .style("text-anchor", "middle")
                .style("fill", "black")
                .attr("alignment-baseline", "auto");
        }

        const ya = axisTypes["y"];
        if (ya) {
            if (!ya.custom) {
                this.y_scale = ya.type === "band" ? scaleBand() : scaleLinear();
                this.y_axis_call = axisLeft(this.y_scale);
            }
            this.y_axis_label = this.y_axis_svg
                .append("text")
                .attr("text-anchor", "middle")
                .style("fill", "black")
                .attr("alignment-baseline", "auto");
        }
        const rya = axisTypes["ry"];
        if (rya?.label) {
            this.ry_axis_label = this.ry_axis_svg
                .append("text")
                .attr("text-anchor", "middle")
                .style("fill", "black")
                .attr("alignment-baseline", "text-top");
        }
    }

    _setMargins() {
        const ax = this.config.axis;
        this.margins = {
            top: ax.tx ? ax.tx.size : 10,
            right: ax.ry ? ax.ry.size : 10,
            bottom: ax.x ? ax.x.size : 10,
            left: ax.y ? ax.y.size : 10,
        };
    }

    _getContentDimensions() {
        //this can end up with -ve width/height e.g. in the process of popping out,
        //or when gridstack does something weird on window resize.
        //we avoid this in WGLChart, but may want a more general way of making things robust.
        return {
            top: this.margins.top,
            left: this.margins.left,
            height: this.height - this.margins.bottom - this.margins.top,
            width: this.width - this.margins.left - this.margins.right,
        };
    }

    setAxisSize(axis, size) {
        switch (axis) {
            case "tx":
                this.margins.top = size;
                break;

            case "ry":
                this.config.axis.ry.size = size;
                this.margins.right = size;
                break;

            case "y":
                this.config.axis.y.size = size;
                this.margins.left = size;
                break;

            case "x":
                this.config.axis.x.size = size;
                this.margins.bottom = size;
        }
        this.graph_area.attr(
            "transform",
            `translate(${this.margins.left},${this.margins.top})`,
        );
    }

    updateAxis() {
        const dim = this._getContentDimensions();
        const ax = this.config.axis;

        this.ry_axis_svg.attr(
            "transform",
            `translate(${dim.left + dim.width},${dim.top})`,
        );
        this.tx_axis_svg.attr("transform", `translate(${dim.left},0)`);
        this.x_axis_svg.attr(
            "transform",
            `translate(${dim.left},${this.height - this.margins.bottom})`,
        );
        this.y_axis_svg.attr("transform", `translate(${dim.left},${dim.top})`);
        if (ax.x) {
            if (this.x_scale) {
                //PJT- TODO fix over-crowded ticks.
                this.x_scale.range([0, dim.width]);
                this.x_axis_svg
                    .selectAll(".tick text")
                    .attr("font-family", "Helvetica");
                this.x_axis_svg.transition().call(this.x_axis_call);
                if (ax.x.rotate_labels) {
                    this.x_axis_svg
                        .selectAll(".tick text")
                        .style("text-anchor", "end")
                        .attr("dx", "-.8em")
                        .attr("dy", ".10em")
                        .attr("transform", "rotate(-45)");
                } else {
                    this.x_axis_svg
                        .selectAll(".tick text")
                        .style("text-anchor", null)
                        .attr("dx", null)
                        .attr("dy", ".71em")
                        .attr("transform", null);
                }
                this.x_axis_svg
                    .selectAll(".tick text")
                    .attr("font-size", ax.x.tickfont);
            }
            this.x_axis_label
                .attr("x", dim.width / 2)
                .attr("y", this.margins.bottom - 4)
                .text(this.config.axis.x.label)
                .attr("font-size", `${this.config.axis.x.textsize}px`)
                .attr("font-family", "Helvetica");
        }

        if (ax.y) {
            if (this.y_scale) {
                this.y_scale.range([0, dim.height]);
                const tinterval = dim.height / this.y_scale.domain().length;
                if (tinterval < 8) {
                    const ii = Math.round(
                        this.y_scale.domain().length / (dim.height / 10),
                    );
                    this.y_axis_call.tickFormat((d, i) => {
                        return i % ii === 0 ? d : "";
                    });
                } else {
                    this.y_axis_call.tickFormat(null);
                }

                this.y_axis_svg.transition().call(this.y_axis_call);
                if (ax.y.rotate_labels) {
                    this.y_axis_svg
                        .selectAll(".tick text")
                        .style("text-anchor", "end")
                        .attr("dy", "1em")
                        .attr("transform", "rotate(+45)");
                } else {
                    this.y_axis_svg
                        .selectAll(".tick text")
                        .style("text-anchor", null)
                        .attr("dx", null)
                        .attr("dy", "0.32em")
                        .attr("transform", null);
                }
                this.y_axis_svg
                    .selectAll(".tick text")
                    .attr("font-size", ax.y.tickfont)
                    .attr("font-family", "Helvetica");
            }

            this.y_axis_label
                .attr("x", -dim.height / 2)
                .attr("y", -dim.left + 10)
                .attr("transform", "rotate(-90)")
                .text(this.config.axis.y.label)
                .attr("font-size", `${this.config.axis.y.textsize}px`)
                .attr("font-family", "Helvetica");
        }

        if (this.ry_axis_label) {
            this.ry_axis_label
                .attr("x", dim.height / 2)
                .attr("y", -this.margins.right + 15)
                .attr("transform", "rotate(90)")
                .text(this.config.axis.ry.label)
                .attr("font-size", `${this.config.axis.ry.textsize}px`);
        }
    }

    setSize(x, y) {
        super.setSize(x, y);
        if (this.svg) {
            this.svg.attr("width", this.width).attr("height", this.height);
        }
    }

    drawXTree(nodeData, number, height, offset) {
        const trans = select(this.contentDiv)
            .transition()
            .duration(400)
            .ease(easeLinear);
        const b = this._getContentDimensions();
        const treeMap = cluster().size([b.width, height]);
        const nodes = treeMap(nodeData);
        const desc = nodes.descendants();
        const interval = b.width / number;
        for (const n of desc) {
            if (!n.children) {
                n.x = interval / 2 + interval * n.data.order;
            }
        }
        for (const n of desc) {
            if (!n.children) {
                let i = n.parent;
                while (i) {
                    i.x =
                        i.children[0].x +
                        (i.children[1].x - i.children[0].x) / 2;
                    i = i.parent;
                }
            }
        }

        const links = this.tx_axis_svg
            .selectAll(".tree-link")
            .data(desc.slice(1))
            .join("path")
            .attr("class", "tree-link")
            .style("stroke", "currentColor")
            .style("fill", "none")
            .style("stroke-width", "1px")
            .transition("trans")
            .attr("d", (d) => {
                let p = d.x;
                if (!d.children) {
                    p = interval / 2 + interval * d.data.order;
                }
                return `M ${p} ${d.y + offset} L${p} ${d.parent.y + offset} L${d.parent.x} ${d.parent.y + offset}`;
            });
    }

    drawYTree(nodeData, number) {
        // trans
        select(this.contentDiv)
            .transition()
            .duration(400)
            .ease(easeLinear);
        const b = this._getContentDimensions();
        const width = this.margins.right;
        const treeMap = cluster().size([b.height, width]);
        const nodes = treeMap(nodeData);
        const desc = nodes.descendants();
        const interval = b.height / number;
        for (const n of desc) {
            if (!n.children) {
                n.x = interval / 2 + interval * n.data.order;
            }
        }
        for (const n of desc) {
            if (!n.children) {
                let i = n.parent;
                while (i) {
                    i.x =
                        i.children[0].x +
                        (i.children[1].x - i.children[0].x) / 2;
                    i = i.parent;
                }
            }
        }

        // links
        this.ry_axis_svg
            .selectAll(".tree-link")
            .data(desc.slice(1))
            .join("path")
            .attr("class", "tree-link")
            .style("stroke", "currentColor")
            .style("fill", "none")
            .style("stroke-width", "1px")
            .transition("trans")
            .attr("d", (d) => {
                let p = d.x;
                if (!d.children) {
                    p = interval / 2 + interval * d.data.order;
                }
                return `M ${width - d.y} ${p} L${width - d.parent.y} ${p} L${width - d.parent.y} ${d.parent.x}`;
            });
    }

    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        let arr = [];
        if (this.config.axis.y) {
            arr = arr.concat([
                {
                    type: "check",
                    current_value: c.axis.y.rotate_labels,
                    label: "Rotate Y axis labels",
                    func: (x) => {
                        c.axis.y.rotate_labels = x;
                        this.setSize();
                    },
                },
                {
                    type: "slider",
                    max: 20,
                    min: 4,
                    step: 1,
                    doc: this.__doc__,
                    current_value: c.axis.y.tickfont,
                    label: "Y axis text size",
                    func: (x) => {
                        c.axis.y.tickfont = x;
                        this.setSize();
                    },
                },
                {
                    type: "slider",
                    max: 200,
                    min: 20,
                    step: 1,
                    doc: this.__doc__,
                    current_value: c.axis.y.size,
                    label: "Y axis width",
                    func: (x) => {
                        this.setAxisSize("y", x);
                        this.setSize();
                    },
                },
            ]);
        }
        if (this.config.axis.x) {
            arr = arr.concat([
                {
                    type: "check",
                    current_value: c.axis.x.rotate_labels,
                    label: "Rotate X axis labels",
                    func: (x) => {
                        c.axis.x.rotate_labels = x;
                        this.updateAxis();
                    },
                },
                {
                    type: "slider",
                    max: 20,
                    min: 4,
                    step: 1,
                    doc: this.__doc__,
                    current_value: c.axis.x.tickfont,
                    label: "X axis text size",
                    func: (x) => {
                        c.axis.x.tickfont = x;
                        this.setSize();
                    },
                },

                {
                    type: "slider",
                    max: 200,
                    min: 20,
                    step: 1,
                    doc: this.__doc__,
                    current_value: c.axis.x.size,
                    label: "X axis Height",
                    func: (x) => {
                        this.setAxisSize("x", x);
                        this.setSize();
                    },
                },
            ]);

            if (this.config.axis.ry) {
                arr = arr.concat([
                    {
                        type: "slider",
                        max: 200,
                        min: 20,
                        step: 1,
                        doc: this.__doc__,
                        current_value: c.axis.x.size,
                        label: "Right Y Axis size",
                        func: (x) => {
                            this.setAxisSize("ry", x);
                            this.setSize();
                        },
                    },
                ]);
            }
        }
        return settings.concat([{
            type: "folder",
            label: "Axis controls",
            current_value: arr,
            func: (x) => {},
        }]);
    }

    _addLegendToSVG(param) {
        const legend = this[param];
        if (!legend) {
            return;
        }
        const cl = [legend.offsetLeft, legend.offsetTop];
        legend._leg = legend.querySelector("g");
        legend._leg_g = select(legend._leg).attr(
            "transform",
            `translate(${cl[0]},${cl[1]})`,
        );
        this.svg.node().append(legend._leg_g.node());
    }

    _removeLegendFromSVG(param) {
        const legend = this[param];
        if (!legend) {
            return;
        }
        legend._leg_g.attr("transform", "translate(0,0)");
        this.legend.querySelector("svg").append(legend._leg);
        this.legend._leg = undefined;
        this.legend._leg_g = undefined;
    }

    getImage(callback, type) {
        if (this.addToImage) {
            this.addToImage();
        }

        //temporarily attach legend
        this._addLegendToSVG("legend");
        if (this.extra_legends) {
            this.extra_legends.forEach((x) => this._addLegendToSVG(x));
        }

        const svgAsXML = new XMLSerializer().serializeToString(this.svg.node());

        if (type === "png") {
            this.getImageFromSVG(this.svg.node(), (canvas) => {
                callback(canvas);
                if (this.removeFromImage) {
                    this.removeFromImage();
                }
                this._removeLegendFromSVG("legend");
                if (this.extra_legends) {
                    this.extra_legends.forEach((x) =>
                        this._removeLegendFromSVG(x),
                    );
                }
            });
        } else {
            callback(svgAsXML);
            if (this.removeFromImage) {
                this.removeFromImage();
            }
            this._removeLegendFromSVG("legend");
            if (this.extra_legends) {
                this.extra_legends.forEach((x) => this._removeLegendFromSVG(x));
            }
        }
    }
}

export default SVGChart;
