import {
    createEl,
    createSVGEl,
    makeResizable,
    makeDraggable,
} from "./Elements.js";
import { select } from "d3-selection";
import { scaleLinear } from "d3-scale";
import { axisBottom } from "d3-axis";
import { getRandomString } from "./Utilities";

function getColorLegendCustom(scale, config = {}) {
    const ticks = config.tickValues || scale.ticks(config.ticks || 4);
    const widths = ticks.map((x) => scale(x));
    return getColorLegend(widths, ticks, config);
}

function getColorLegend(colors, names, config = {}) {
    const len = colors.length;
    const type = config.type || "color";
    const h_fac = type === "circle" ? (config.itemGap || 24) : 12;
    const itemTop = config.itemTop || 0;
    const dynamicCircleSpacing = type === "circle" && config.dynamicCircleSpacing;
    const circleGap = config.circleGap ?? 4;
    const baseCy = itemTop + 6;
    const itemY = [];
    let contentBottom = 0;
    if (dynamicCircleSpacing) {
        let currentCy = baseCy;
        for (let i = 0; i < len; i++) {
            const radius = Number(colors[i]) || 0;
            if (i > 0) {
                const prevRadius = Number(colors[i - 1]) || 0;
                currentCy += prevRadius + radius + circleGap;
            }
            itemY.push(currentCy);
            contentBottom = Math.max(contentBottom, currentCy + radius);
        }
    } else {
        for (let i = 0; i < len; i++) {
            const y = itemTop + (i + 1) * 2 + i * h_fac;
            itemY.push(y + 4);
            contentBottom = Math.max(contentBottom, y + h_fac);
        }
    }
    const height = Math.max(len * h_fac + (len + 2) * 2 + itemTop, contentBottom + 8);
    const container = createEl("div", {
        styles: {
            width: `${config.width || 120}px`,
            height: `${(height > 215 ? 215 : height) + 35}px`,
            position: "absolute",
            border: "0.5px solid black",
        },
        classes: ["legend-container"],
    });

    createEl(
        "div",
        {
            styles: {
                height: "20px",
                whiteSpace: "nowrap",
            },
            classes: ["legend-title"],
            text: config.label,
        },
        container,
    );
    const body = createEl(
        "div",
        {
            styles: {
                overflowY: "auto",
                overflowX: "hidden",
                height: "calc(100% - 10px)",
                width: "100%",
            },
            classes: ["legend-body"],
        },
        container,
    );

    const legend = createSVGEl(
        "svg",
        {
            height: height,
            width: 180,
            styles: {
                position: "relative",
            },
        },
        body,
    );

    const legendg = createSVGEl("g", {}, legend);
    const defaultTextOffset = type === "color" ? 14 : type === "line" ? 20 : 30;
    const t_offset = config.textOffset || defaultTextOffset;
    const circleX = config.circleX || 15;

    for (let i = 0; i < len; i++) {
        if (type === "color") {
            const y = itemTop + (i + 1) * 2 + i * h_fac;
            createSVGEl(
                "rect",
                {
                    y,
                    x: 2,
                    height: "10",
                    width: "10",
                    styles: {
                        fill: colors[i],
                    },
                },
                legendg,
            );
        } else if (type === "line") {
            const y = itemTop + (i + 1) * 2 + i * h_fac;
            createSVGEl(
                "rect",
                {
                    y,
                    x: 2,
                    height: "10",
                    width: colors[i],
                    styles: {
                        fill: "gray",
                    },
                },
                legendg,
            );
        } else if (type === "circle") {
            const cy = itemY[i];
            createSVGEl(
                "circle",
                {
                    cy,
                    cx: circleX,
                    r: colors[i],
                    styles: {
                        fill: "gray",
                    },
                },
                legendg,
            );
        }

        const t = createSVGEl(
            "text",
            {
                y: type === "circle" ? itemY[i] + 5 : itemTop + (i + 1) * 2 + i * h_fac + 9,
                x: t_offset,
                //Firefox was wrong. Changing to "center" made Chrome wrong in same way
                // "alignment-baseline":"middle",
                //so... since change to "center" makes them consistent, adjusting "+6" offset above by "+3" to compensate
                "alignment-baseline": "center",
                styles: {
                    "font-size": "12px",
                    fill: "currentcolor",
                },
                text: names[i] === "" ? "none" : names[i],
            },
            legendg,
        );
        select(t).style("fill", "currentcolor"); //seems redundant?
    }

    makeDraggable(container, { handle: ".legend-body" });
    makeResizable(container);
    return container;
}

function getColorBar(colors, config = {}) {
    const c = config;
    const len = colors.length;
    const colorPct = colors.map((c, i) => [
        `${Math.floor((i / len) * 100)}%`,
        c,
    ]);
    const width = c.width || 120;
    const height = c.height || 45;
    const svg = createSVGEl("svg", {
        height: height,
        width: width,
        styles: {
            position: "relative",
        },
    });
    const g = select(svg).append("g");
    const id = getRandomString();
    if (c.label) {
        g.append("text")
            .text(c.label)
            .attr("x", 10)
            .style("fill", "currentColor")
            .attr("alignment-baseline", "hanging")
            .style("font-size", "12px");
    }

    const grad = g
        .append("defs")
        .append("linearGradient")
        .attr("id", id)
        .attr("x1", "0%") // bottom
        .attr("y1", "0%")
        .attr("x2", "100%") // to top
        .attr("y2", "0%")
        .attr("spreadMethod", "pad");
    let bH = height;
    if (c.label) {
        bH -= 10;
    }
    if (c.range) {
        bH -= 25;
    }

    const bar = g
        .append("rect")
        .attr("x", c.range ? 10 : 0)
        .attr("y", c.label ? 10 : 0)
        .attr("width", c.range ? width - 20 : 20)
        .attr("height", bH)
        .style("fill", `url(#${id})`);

    colorPct.forEach((d) => {
        grad.append("stop")
            .attr("offset", d[0])
            .attr("stop-color", d[1])
            .attr("stop-opacity", 1);
    });

    //add bottom axis
    if (c.range) {
        const scale = scaleLinear()
            .domain([c.range[0], c.range[1]])
            .range([0, width - 20]);
        const axis = axisBottom(scale)
            .tickFormat((v, i) =>
                v >= 10000 ? Number.parseFloat(v).toPrecision(2) : v,
            )
            .ticks(Math.ceil(width / 20));

        const axisg = g
            .append("g")
            .call(axis)
            .attr("transform", `translate(10,${height - 25})`)
            .selectAll("text")
            .style("text-anchor", "end")
            .attr("dx", "-.4em")
            .attr("dy", ".4em")
            .attr("transform", "rotate(-45)");
    }
    const container = createEl("div", {
        height: height,
        width: width,
        styles: {
            position: "absolute",
        },
        classes: ["legend-container"],
    });
    container.append(svg);
    makeDraggable(container);
    return container;
}

export { getColorLegend, getColorBar, getColorLegendCustom };
