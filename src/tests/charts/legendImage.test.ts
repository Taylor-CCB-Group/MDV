import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import { addLegendsToSvg } from "@/charts/legendImage";

const SVG_NS = "http://www.w3.org/2000/svg";

/** Where elements sit on screen; jsdom has no layout. */
const boxes = new Map<Element, { left: number; top: number; width: number; height: number }>();

function placeAt(element: Element, left: number, top: number, width = 100, height = 50) {
    boxes.set(element, { left, top, width, height });
}

/** A legend as `LegendWrapper` draws it: a host holding a frame, which holds the legend's SVG. */
function wrappedLegend(left: number, top: number, minimized = false) {
    const host = document.createElement("div");
    const frame = document.createElement("div");
    frame.className = "legend-container";
    host.append(frame);
    if (!minimized) {
        const svg = document.createElementNS(SVG_NS, "svg");
        const g = document.createElementNS(SVG_NS, "g");
        g.setAttribute("transform", "translate(2,3)");
        const after = document.createElementNS(SVG_NS, "text");
        svg.append(g, after);
        frame.append(svg);
        placeAt(svg, left, top);
    }
    document.body.append(host);
    return host;
}

/** The chart's SVG, at (10, 100) on screen. */
function chartSvg() {
    const svg = document.createElementNS(SVG_NS, "svg");
    svg.setAttribute("width", "300");
    svg.setAttribute("height", "200");
    document.body.append(svg);
    placeAt(svg, 10, 100, 300, 200);
    return svg;
}

beforeEach(() => {
    vi.spyOn(Element.prototype, "getBoundingClientRect").mockImplementation(function (this: Element) {
        const box = boxes.get(this) ?? { left: 0, top: 0, width: 0, height: 0 };
        return DOMRect.fromRect({ x: box.left, y: box.top, width: box.width, height: box.height });
    });
});

afterEach(() => {
    document.body.innerHTML = "";
    boxes.clear();
    vi.restoreAllMocks();
});

describe("addLegendsToSvg", () => {
    test("puts a floating legend's drawing where it sits over the plot", () => {
        const legend = wrappedLegend(50, 125);
        const svg = chartSvg();

        addLegendsToSvg(svg, [legend]);

        const drawing = svg.querySelector("g");
        expect(drawing?.getAttribute("transform")).toBe("translate(40,25) translate(2,3)");
        expect(legend.querySelector("g")).toBeNull();
        // nothing above the plot, so the image keeps its size
        expect(svg.getAttribute("height")).toBe("200");
        expect(svg.hasAttribute("viewBox")).toBe(false);
    });

    test("draws legends in the strip above the plot, making the image taller to fit them", () => {
        const legend = wrappedLegend(10, 60);
        const svg = chartSvg();

        addLegendsToSvg(svg, [legend]);

        expect(svg.querySelector("g")?.getAttribute("transform")).toBe("translate(0,-40) translate(2,3)");
        expect(svg.getAttribute("height")).toBe("240");
        expect(svg.getAttribute("viewBox")).toBe("0 -40 300 240");
    });

    test("puts everything back as it was", () => {
        const floating = wrappedLegend(50, 125);
        const above = wrappedLegend(10, 60);
        const svg = chartSvg();

        const restore = addLegendsToSvg(svg, [floating, above]);
        restore();

        for (const legend of [floating, above]) {
            const drawing = legend.querySelector("g");
            expect(drawing?.getAttribute("transform")).toBe("translate(2,3)");
            // before the text that followed it
            expect(drawing?.nextElementSibling?.tagName).toBe("text");
        }
        expect(svg.childElementCount).toBe(0);
        expect(svg.getAttribute("height")).toBe("200");
        expect(svg.hasAttribute("viewBox")).toBe(false);
    });

    test("leaves out minimised and missing legends", () => {
        const svg = chartSvg();
        addLegendsToSvg(svg, [wrappedLegend(10, 10, true), undefined, null]);
        expect(svg.childElementCount).toBe(0);
    });

    test("places an older legend, whose drawing has no transform, by its SVG", () => {
        const legend = document.createElement("div");
        const inner = document.createElementNS(SVG_NS, "svg");
        inner.append(document.createElementNS(SVG_NS, "g"));
        legend.append(inner);
        document.body.append(legend);
        placeAt(inner, 17, 109);
        const svg = chartSvg();

        addLegendsToSvg(svg, [legend]);

        expect(svg.querySelector("g")?.getAttribute("transform")).toBe("translate(7,9)");
    });
});
