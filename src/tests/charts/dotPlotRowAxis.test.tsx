import { render } from "@testing-library/react";
import { scaleBand } from "@visx/scale";
import { describe, expect, test } from "vitest";
import DotPlotRowAxis from "@/react/components/dotplot/DotPlotRowAxis";

/** Renders the row axis for `count` rows in a plot `height` pixels tall; returns its labels. */
function renderRowAxis(count: number, height: number, tickfont = 10) {
    const rows = Array.from({ length: count }, (_, i) => i);
    const labels = rows.map((i) => `cell state ${i}`);
    const scale = scaleBand({ domain: rows, range: [0, height] });
    const { container } = render(
        <svg>
            <DotPlotRowAxis rows={rows} labels={labels} scale={scale} left={100} top={0} tickfont={tickfont} />
        </svg>,
    );
    return [...container.querySelectorAll(".visx-axis-tick text")];
}

describe("DotPlotRowAxis", () => {
    test("labels every row, not just about 10 of them", () => {
        // the Myeloid cell states: 22 rows used to get only 11 labels
        const labels = renderRowAxis(22, 440);
        expect(labels.map((t) => t.textContent)).toEqual(Array.from({ length: 22 }, (_, i) => `cell state ${i}`));
    });

    test("labels every row of a long list", () => {
        expect(renderRowAxis(67, 1340)).toHaveLength(67);
    });

    test("keeps the chosen font size when rows are tall enough", () => {
        const [label] = renderRowAxis(22, 440, 10); // 20px rows
        expect(label.getAttribute("font-size")).toBe("10");
    });

    test("shrinks labels to the row height so they don't overlap", () => {
        const [label] = renderRowAxis(22, 176, 10); // 8px rows
        expect(label.getAttribute("font-size")).toBe("8");
    });
});
