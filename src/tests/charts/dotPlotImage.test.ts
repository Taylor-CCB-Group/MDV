import { afterEach, describe, expect, it, vi } from "vitest";
import { buildChartSvg } from "@/react/components/dotplot/dotPlotImage";

/** jsdom has no layout, so give elements a box by their `data-box` attribute. */
function mockBoxes() {
    vi.spyOn(Element.prototype, "getBoundingClientRect").mockImplementation(function (this: Element) {
        const [x, y, width, height] = (this.getAttribute("data-box") ?? "0 0 0 0").split(" ").map(Number);
        return DOMRect.fromRect({ x, y, width, height });
    });
}

function makeChart() {
    const container = document.createElement("div");
    container.setAttribute("data-box", "100 50 400 300");
    container.innerHTML = `
        <svg data-box="100 80 400 270" class="absolute" style="fill: red">
            <g><circle r="5" class="fill-current"></circle></g>
            <svg data-box="120 300 20 10"><text>CD4</text></svg>
        </svg>
        <div class="legend-container">
            <div class="legend-title" data-box="310 60 80 20">% expressing</div>
            <div class="legend-body"><svg data-box="310 80 80 30"><circle r="3"></circle></svg></div>
        </div>
        <svg data-box="0 0 0 0"></svg>`;
    document.body.append(container);
    return container;
}

describe("buildChartSvg", () => {
    afterEach(() => {
        vi.restoreAllMocks();
        document.body.innerHTML = "";
    });

    it("nests each top-level SVG at its place in the chart", () => {
        mockBoxes();
        const out = buildChartSvg(makeChart());
        expect(out.getAttribute("width")).toBe("400");
        expect(out.getAttribute("height")).toBe("300");
        // the plot and the legend body; the label's inner SVG comes with the plot, the empty one is left out
        const nested = [...out.children].filter((el) => el.tagName === "svg");
        expect(nested.map((el) => [el.getAttribute("x"), el.getAttribute("y")])).toEqual([
            ["0", "30"],
            ["210", "30"],
        ]);
        expect(nested[0].querySelectorAll("svg")).toHaveLength(1);
    });

    it("inlines computed styles and drops classes", () => {
        mockBoxes();
        const out = buildChartSvg(makeChart());
        const plot = out.querySelector("svg");
        expect(plot?.hasAttribute("class")).toBe(false);
        expect(plot?.getAttribute("style")).toContain("fill:");
        expect(out.querySelector("circle")?.hasAttribute("class")).toBe(false);
    });

    it("draws HTML legend titles as SVG text", () => {
        mockBoxes();
        const out = buildChartSvg(makeChart());
        const title = [...out.children].find((el) => el.tagName === "text");
        expect(title?.textContent).toBe("% expressing");
        expect(title?.getAttribute("x")).toBe("210");
        expect(title?.getAttribute("y")).toBe("20");
    });

    it("takes in the strip of legends above the plot, growing the image to cover both", () => {
        mockBoxes();
        const chart = makeChart();
        const strip = document.createElement("div");
        strip.setAttribute("data-box", "100 10 400 40");
        strip.innerHTML = `<div class="legend-title" data-box="104 12 80 16">Cell type</div>
            <svg data-box="104 28 200 20"><rect width="10" height="10"></rect></svg>`;
        document.body.append(strip);

        const out = buildChartSvg(chart, strip);

        expect(out.getAttribute("height")).toBe("340");
        const nested = [...out.children].filter((el) => el.tagName === "svg");
        expect(nested.map((el) => el.getAttribute("y"))).toEqual(["70", "70", "18"]);
        const titles = [...out.children].filter((el) => el.tagName === "text").map((el) => el.textContent);
        expect(titles).toEqual(["% expressing", "Cell type"]);
    });

    it("cuts a legend scrolled in the strip to the part that can be seen", () => {
        mockBoxes();
        const chart = makeChart();
        const strip = document.createElement("div");
        strip.setAttribute("data-box", "100 10 400 40");
        // a 300px wide legend in a body that shows 200px of it, scrolled right by 50px
        strip.innerHTML = `<div class="legend-body" style="overflow-x: auto" data-box="104 28 200 20">
            <svg data-box="54 28 300 20"><rect width="10" height="10"></rect></svg></div>`;
        document.body.append(strip);

        const out = buildChartSvg(chart, strip);

        const frame = [...out.children].filter((el) => el.tagName === "svg")[2];
        expect(["x", "y", "width", "height"].map((name) => frame.getAttribute(name))).toEqual(["4", "18", "200", "20"]);
        const legend = frame.querySelector("svg");
        expect(["x", "y", "width", "height"].map((name) => legend?.getAttribute(name))).toEqual([
            "-50",
            "0",
            "300",
            "20",
        ]);
    });
});
