import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import { addLegendsToSvg } from "@/charts/legendImage";

const SVG_NS = "http://www.w3.org/2000/svg";

/** A legend as `LegendWrapper` draws it: a host holding a frame, which holds the legend's SVG. */
function wrappedLegend(left: number, top: number, minimized = false) {
    const host = document.createElement("div");
    const frame = document.createElement("div");
    frame.className = "legend-container";
    frame.style.left = `${left}px`;
    frame.style.top = `${top}px`;
    host.append(frame);
    if (!minimized) {
        const svg = document.createElementNS(SVG_NS, "svg");
        const g = document.createElementNS(SVG_NS, "g");
        g.setAttribute("transform", "translate(2,3)");
        const after = document.createElementNS(SVG_NS, "text");
        svg.append(g, after);
        frame.append(svg);
    }
    document.body.append(host);
    return host;
}

function chartSvg() {
    const svg = document.createElementNS(SVG_NS, "svg");
    document.body.append(svg);
    return svg;
}

beforeEach(() => {
    // jsdom has no layout, so offsets come from inline styles
    vi.spyOn(HTMLElement.prototype, "offsetLeft", "get").mockImplementation(function (this: HTMLElement) {
        return Number.parseFloat(this.style.left) || 0;
    });
    vi.spyOn(HTMLElement.prototype, "offsetTop", "get").mockImplementation(function (this: HTMLElement) {
        return Number.parseFloat(this.style.top) || 0;
    });
});

afterEach(() => {
    document.body.innerHTML = "";
    vi.restoreAllMocks();
});

describe("addLegendsToSvg", () => {
    test("puts a legend's drawing where its frame sits over the plot", () => {
        const legend = wrappedLegend(40, 25);
        const svg = chartSvg();

        addLegendsToSvg(svg, [legend]);

        const drawing = svg.querySelector("g");
        expect(drawing?.getAttribute("transform")).toBe("translate(40,25)");
        expect(legend.querySelector("g")).toBeNull();
    });

    test("puts the drawing back where it was, as it was", () => {
        const legend = wrappedLegend(40, 25);
        const svg = chartSvg();

        const restore = addLegendsToSvg(svg, [legend]);
        restore();

        const drawing = legend.querySelector("g");
        expect(drawing?.getAttribute("transform")).toBe("translate(2,3)");
        // before the text that followed it
        expect(drawing?.nextElementSibling?.tagName).toBe("text");
        expect(svg.childElementCount).toBe(0);
    });

    test("leaves out minimised and missing legends", () => {
        const svg = chartSvg();
        addLegendsToSvg(svg, [wrappedLegend(10, 10, true), undefined, null]);
        expect(svg.childElementCount).toBe(0);
    });

    test("uses an older legend's own position when it has no frame", () => {
        const legend = document.createElement("div");
        legend.style.left = "7px";
        legend.style.top = "9px";
        const inner = document.createElementNS(SVG_NS, "svg");
        inner.append(document.createElementNS(SVG_NS, "g"));
        legend.append(inner);
        document.body.append(legend);
        const svg = chartSvg();

        addLegendsToSvg(svg, [legend]);

        expect(svg.querySelector("g")?.getAttribute("transform")).toBe("translate(7,9)");
    });
});
