import SVGChart from "./SVGChart.js";
import BaseChart from "./BaseChart";
import { stratify as d3Stratify, scaleLinear, tree as d3Tree, arc } from "d3";
import { getColorLegend } from "../utilities/Color.js";
import { select } from "d3";

class CellRadialChart extends SVGChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config, {});
        const c = this.config;
        c.label_size = c.label_size || 12;

        this.stratify = d3Stratify()
            .id((d) => d[0])
            .parentId((d) => d[1]);
        //thickness
        this.linkThicknessScale = scaleLinear();
        c.scales = {};
        if (!c.link_thickness) {
            c.link_thickness = {
                type: "constant",
                width: 3,
            };
            const tq = this.dataStore.getColumnQuantile(c.param[4], 0.01);
            this._changeLinkThicknessScale(tq, [1, 10]);
        }
        c.state_legend = c.state_legend || {
            show: true,
            position: [10, 20],
        };
        c.outer_cell_color = c.outer_cell_color || "state";
        this.extra_legends = ["state_legend"];
        c.node_size = c.node_size || { inner: 5, outer: 5 };

        c.center_cells =
            c.center_cells ||
            this.dataStore.getColumnValues(c.param[1]).slice(0, 2);
        c.states =
            c.states || this.dataStore.getColumnValues(c.param[0]).slice(0, 2);

        //node size
        /*this.nodeScale= scaleSqrt();
        if (!c.node_size){
            c.node_size={};
            const tq= this.dataStore.getMinMaxForColumn(c.param[5]);
            this._changeNodeSizeScale(tq,[1,25]);
        }

        if (!c.link_color){
            c.link_color={};
            const tq= this.dataStore.getColumnQuantile(c.param[6],0.01);
            this._changeLinkColorScale(tq);

        }
        if (!c.color_legend){
            c.color_legend={display:true}
        }

        c.node_color = c.node_color || "cells"
        */

        c.label_size = c.label_size || 13;
        this.makeHierarchy();
        this.drawTree();
        this.setStateLegend();
    }

    onDataFiltered() {
        this.makeHierarchy();
        this.drawTree();
        this.setStateLegend();
    }

    getConfig() {
        const config = super.getConfig();
        const sl = this.state_legend;
        if (sl) {
            config.state_legend.position = [sl.offsetLeft, sl.offsetTop];
        }
        return config;
    }

    setStateLegend() {
        const c = this.config;
        const sl = this.config.state_legend;
        if (this.state_legend) {
            this.state_legend.remove();
        }
        if (!sl.show) {
            return;
        }
        const col = this.dataStore.columnIndex[c.param[0]];
        let names = col.values;

        let colors = this.dataStore.getColumnColors(c.param[0]);

        if (c.states) {
            colors = colors.filter((x, i) => c.states.indexOf(names[i]) !== -1);
            names = names.filter((x) => c.states.indexOf(x) !== -1);
        }
        this.state_legend = getColorLegend(colors, names, { label: col.name });
        this.contentDiv.append(this.state_legend);
        this.state_legend.style.top = `${sl.position[1]}px`;
        this.state_legend.style.left = `${sl.position[0]}px`;
    }

    _changeNodeSize() {
        const c = this.config;
        this.nodes.attr("r", (d) => {
            if (d.depth === 1) {
                return c.node_size.inner;
            }
            if (d.depth === 3) {
                return c.node_size.outer;
            }
        });
    }

    _changeLinkThicknessScale(domain, range, update) {
        const c = this.config;
        const tdata = this.dataStore.getRawColumn(c.param[4]);
        if (range) {
            c.link_thickness.range = range;
        }
        if (domain) {
            c.link_thickness.domain = domain;
        }

        this.linkThicknessScale
            .domain(c.link_thickness.domain)
            .range(c.link_thickness.range)
            .clamp(true);
        const cwidth = c.link_thickness.type === "constant";
        const width = `${c.link_thickness.width}px`;
        if (update) {
            this.links.attr("stroke-width", (d) => {
                if (d.depth === 3) {
                    return cwidth
                        ? width
                        : this.linkThicknessScale(tdata[d.id]);
                }
                return width;
            });
        }
    }

    setLinkThicknessColumn(col) {
        this.config.param[4] = col;
        const tq = this.dataStore.getColumnQuantile(col, 0.01);
        this._changeLinkThicknessScale(tq, null, true);
    }

    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        const c1 = this.dataStore.getColumnValues(c.param[1]);
        const st = this.dataStore.getColumnValues(c.param[0]);
        const tCols = this.dataStore.getColumnList("number");
        const tQuant = this.dataStore.getColumnQuantile(c.param[4], 0.01);

        const s = settings.concat([
            {
                type: "slider",
                max: 20,
                min: 5,
                doc: this.__doc__,
                current_value: c.label_size,
                label: "Label Size",
                func: (x) => {
                    c.label_size = x;
                    this.graph_area.selectAll("text").attr("font-size", x);
                },
            },
            {
                type: "multidropdown",
                current_value: c.center_cells,
                values: [c1.map((x) => [x]), 0, 0],
                label: "Central Cells",
                func: (x) => {
                    c.center_cells = x;
                    this.makeHierarchy();
                    this.drawTree();
                },
            },
            {
                type: "multidropdown",
                current_value: c.states,
                values: [st.map((x) => [x]), 0, 0],
                label: "Cell states",
                func: (x) => {
                    c.states = x;
                    this.makeHierarchy();
                    this.drawTree();
                    this.setStateLegend();
                },
            },

            {
                type: "slider",
                label: "Inner Node Size",
                max: 20,
                min: 1,
                doc: this.__doc__,
                current_value: c.node_size.inner,
                func: (x) => {
                    c.node_size.inner = x;
                    this.drawTree();
                },
            },
            {
                type: "check",
                label: "Show Inner Links",

                current_value: c.central_links,
                func: (x) => {
                    c.central_links = x;
                    this.drawTree();
                },
            },
            {
                type: "slider",
                label: "Outer Node Size",
                max: 20,
                min: 1,
                doc: this.__doc__,
                current_value: c.node_size.outer,
                func: (x) => {
                    c.node_size.outer = x;
                    this.drawTree();
                },
            },
            {
                type: "radiobuttons",
                label: "Outer Node Color",
                current_value: c.outer_cell_color,
                choices: [
                    ["state", "state"],
                    ["Cell color", "cell_color"],
                ],
                func: (x) => {
                    c.outer_cell_color = x;
                    this.drawTree();
                },
            },
            {
                type: "radiobuttons",
                label: "Link Thickness Type",
                current_value: c.link_thickness.type,
                choices: [
                    ["constant", "constant"],
                    ["based on column", "column"],
                ],
                func: (x) => {
                    c.link_thickness.type = x;
                    this._changeLinkThicknessScale(null, null, true);
                },
            },
            {
                type: "slider",
                label: "Link Thickness",
                max: 10,
                min: 1,
                doc: this.__doc__,
                current_value: this.config.link_thickness.width,
                func: (x) => {
                    c.link_thickness.width = x;
                    this._changeLinkThicknessScale(null, null, true);
                },
            },
            {
                label: "Link Thickness Column",
                type: "dropdown",
                values: [tCols, "name", "field"],
                current_value: c.param[4],
                func: (v, controls) => {
                    const sl = controls["Link Thickness Domain"].noUiSlider;
                    const tq = this.dataStore.getColumnQuantile(v, 0.01);
                    sl.updateOptions({ range: { max: tq[1], min: tq[0] } });
                    sl.set(tq);
                    this.setLinkThicknessColumn(v);
                },
            },
            {
                type: "doubleslider",
                max: 20,
                min: 1,
                doc: this.__doc__,
                current_value: this.config.link_thickness.range,
                label: "Link Thickness Range",
                func: (x, y) => {
                    this._changeLinkThicknessScale(null, [x, y], true);
                },
            },
            {
                type: "doubleslider",
                max: tQuant[1],
                min: tQuant[0],
                doc: this.__doc__,
                current_value: this.config.link_thickness.domain,
                label: "Link Thickness Domain",
                func: (x, y) => {
                    this._changeLinkThicknessScale([x, y], null, true);
                },
            },
        ]);
        return s;
    }

    drawTree() {
        const b = this._getContentDimensions();
        const c = this.config;
        const thick = this.dataStore.getRawColumn(this.config.param[4]);
        const size_factor = Math.min(b.width, b.height);
        const radius = size_factor / 2.8;
        this.graph_area.attr(
            "transform",
            `translate(${b.left + b.width / 2},${b.top + b.height / 2})`,
        );
        const wedgeColors = this.dataStore.getColumnColors(c.param[0]);
        const c1Colors = this.dataStore.getColumnColors(c.param[1]);
        const c2Colors = this.dataStore.getColumnColors(c.param[2]);
        const stData = this.dataStore.getRawColumn(c.param[0]);

        const tree = d3Tree()

            .size([360, radius])
            .separation((a, b) => (a.parent === b.parent ? 1 : 2) / a.depth);

        const root = tree(this.root);

        this.graph_area.selectAll(".link").remove();
        this.graph_area.selectAll(".node").remove();
        this.graph_area.selectAll(".wedge").remove();
        this.graph_area.selectAll(".ring").remove();
        //this.g.selectAll(".level").remove();

        const level_radii = {};
        const data = root.descendants().slice(1);
        const inc = radius / 10;
        data.forEach((d) => {
            if (d.depth === 1) {
                d.y += inc;
            }
        });
        const wedgeData = {};
        for (const d of data) {
            if (d.depth !== 3) {
                continue;
            }
            const id = d.data[1];
            if (!wedgeData[id]) {
                wedgeData[id] = [d.x, d.x, d.data[3]];
            } else {
                wedgeData[id][0] = Math.min(wedgeData[id][0], d.x);
                wedgeData[id][1] = Math.max(wedgeData[id][1], d.x);
            }
        }
        const wData = [];
        for (const wd in wedgeData) {
            const item = wedgeData[wd];
            const info = wd.split("|");
            wData.push({
                startAngle: ((item[0] - 2) * Math.PI) / 180,
                endAngle: ((item[1] + 2) * Math.PI) / 180,
                cell: info[0],
                state: info[1],
                sid: item[2],
            });
        }

        const arcGen = arc()
            .innerRadius(radius * (2 / 3))
            .outerRadius(radius + 5);
        this.graph_area
            .selectAll(".wedge")
            .data(wData)
            .enter()
            .append("path")
            .attr("class", "wedge")
            .attr("d", arcGen)
            .attr("opacity", 0.8)
            .attr("fill", (d) => wedgeColors[d.sid])
            .attr("stroke", "gray")
            .attr("stroke-width", 1);

        const ringGen = arc()
            .startAngle(0)
            .endAngle(2 * Math.PI);
        const ringData = [
            { outerRadius: radius + 5, innerRadius: radius - 5, color: "blue" },
            {
                outerRadius: (radius * 1) / 3 + inc + 5,
                innerRadius: (radius * 1) / 3 + inc - 5,
                color: "green",
            },
        ];

        this.graph_area
            .selectAll(".ring")
            .data(ringData)
            .enter()
            .append("path")
            .attr("class", "ring")
            .attr("d", ringGen)
            .attr("opacity", 0.4)
            .attr("fill", (d) => d.color)
            .attr("stroke", "gray")
            .attr("stroke-width", 1);

        const cwidth = c.link_thickness.type === "constant";

        this.links = this.graph_area
            .selectAll(".link")
            .data(data)
            .enter()
            .filter((d) => {
                return c.central_links ? true : d.depth > 1;
            })
            .append("path")
            .attr("class", "link")
            .attr("fill", "none")
            .attr("stroke", "currentColor")
            .attr("stroke-width", (d) => {
                if (d.depth === 3) {
                    return cwidth
                        ? c.link_thickness.type
                        : this.linkThicknessScale(thick[d.id]);
                }
                return c.link_thickness.type;
            })
            .attr("d", (d) => {
                level_radii[d.depth] = d.y;
                if (d.depth > 1) {
                    return (
                        // biome-ignore lint/style/useTemplate: clearer to have this on several lines
                        "M" +
                        project(d.x, d.y) +
                        "C" +
                        project(d.x, (d.y + d.parent.y) / 2) +
                        " " +
                        project(d.parent.x, (d.y + d.parent.y) / 2) +
                        " " +
                        project(d.parent.x, d.parent.y)
                    );
                }

                return `M${project(d.x, d.y)}L${project(d.parent.x, d.parent.y)}`;
            });

        const node = this.graph_area
            .selectAll(".node")
            .data(data)
            .enter()
            .append("g")
            .attr(
                "class",
                (d) => `node${d.children ? " node--internal" : " node--leaf"}`,
            )
            .attr("transform", (d) => `translate(${project(d.x, d.y)})`)
            .on("click", (d) => {
                this.nodeClicked(d);
            });

        this.nodes = node
            .filter((d) => d.depth !== 2)
            .append("circle")
            .style("pointer-events", "visible")
            .attr("fill", (d) => {
                if (d.depth === 1) {
                    return c1Colors[d.data[3]];
                }

                return c.outer_cell_color === "state"
                    ? wedgeColors[stData[d.data[0]]]
                    : c2Colors[d.data[4]];
            })
            .attr("stroke", "black")
            .attr("stroke-width", 1)
            .attr("r", (d) => {
                if (d.depth === 1) {
                    return c.node_size.inner;
                }
                if (d.depth === 3) {
                    return c.node_size.outer;
                }
            });

        const fs = `${this.config.label_size}px`;
        node.append("text")
            .attr("dy", ".31em")
            .attr("font-size", fs)
            .style("fill", "currentColor")
            .attr("font-family", "Helvetica")
            .attr("x", (d) => {
                const x_off =
                    d.depth === 1
                        ? c.node_size.inner + 1
                        : c.node_size.outer + 1;
                return d.x < 180 === !d.children ? x_off : -x_off;
            })
            .attr("alignment-baseline", (d) =>
                d.depth === 1 && c.central_links ? "hanging" : "baseline",
            )
            .style("text-anchor", (d) => {
                if (!d.children) {
                    return d.x < 180 ? "start" : "end";
                }

                return d.x > 180 ? "start" : "end";
            })
            .attr("transform", (d) => {
                if (d.depth === 1 && !c.central_links) {
                    return "rotate(0)";
                }
                return `rotate(${d.x < 180 ? d.x - 90 : d.x + 90})`;
            })
            .text((d) => (d.depth === 2 ? "" : d.data[2]));
    }

    setSize(x, y) {
        super.setSize(x, y);
        this.drawTree();
    }

    makeHierarchy() {
        //get cell type columns
        const c = this.config;
        const p = c.param;
        //cell types
        const ct1 = this.dataStore.columnIndex[p[1]];
        const ct2 = this.dataStore.columnIndex[p[2]];
        //state
        const st = this.dataStore.columnIndex[p[0]];

        //data structure [child,parent,other_data....]

        //add the center cells-first ring
        let data = c.center_cells.map((x) => [
            x,
            "root",
            x,
            ct1.values.indexOf(x),
        ]);
        //add the root
        data.unshift(["root", null]);
        //add the second (inner) ring celltype|state
        const innerRingEntries = new Set();
        for (const state of c.states) {
            for (const cell of c.center_cells) {
                const name = `${state}|${cell}`;
                data.push([name, cell, state]);
                innerRingEntries.add(name);
            }
        }
        //get entries in inner ring
        const st_inter = c.center_cells.length + 1;
        const en_inter = st_inter + innerRingEntries.size;

        //to hold which inner ring entries have links
        const hasLinks = new Set();
        const f = this.dataStore.filterArray;

        //go through the data looking for links
        for (let n = 0; n < this.dataStore.size; n++) {
            if (f[n] > 0) {
                continue;
            }
            const c1 = ct1.values[ct1.data[n]];
            const c2 = ct2.values[ct2.data[n]];
            const s = st.values[st.data[n]];
            const name = `${s}|${c1}`;
            //exclude certain interactions no way odf user adding this
            let exclude = false;
            if (c.exclude_interactions) {
                for (const e of c.exclude_interactions) {
                    if (e.indexOf(c1) !== -1) {
                        if (e.indexOf(c2) !== -1) {
                            exclude = true;
                            break;
                        }
                    }
                }
            }
            if (exclude) {
                continue;
            }
            if (innerRingEntries.has(name)) {
                hasLinks.add(name);
                //child is index of link, parern is inner ring id
                data.push([n, name, c2, st.data[n], ct2.data[n]]);
            }
        }
        //filter out cell/state (inner ring) with no interactions
        data = data.filter((x, i) => {
            //mot in inner ring
            if (i < st_inter || i >= en_inter) {
                return true;
            }
            return hasLinks.has(x[0]);
        });
        this.root = this.stratify(data);
    }
}

function project(x, y) {
    const angle = ((x - 90) / 180) * Math.PI;
    const radius = y;
    return [radius * Math.cos(angle), radius * Math.sin(angle)];
}

BaseChart.types["cell_radial_chart"] = {
    name: "Radial Connectivty Map",
    required: ["interactions"],
    methodsUsingColumns: ["setLinkThicknessColumn"],
    init: (config, dataSource, extraControls) => {
        const i = dataSource.interactions;
        const param = dataSource.interactions.cell_radial_chart;
        config.param = [
            i.pivot_column,
            i.interaction_columns[0],
            i.interaction_columns[1],
            i.pivot_column, //stat cutoff - not used
            param.link_thickness,
        ];
        Object.assign(config, param.default_parameters);
    },
    class: CellRadialChart,
};
