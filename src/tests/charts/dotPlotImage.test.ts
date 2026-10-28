import { buildChartSvg } from "@/react/components/dotplot/dotPlotImage";
import { afterEach, describe, expect, it, vi } from "vitest";

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
});
