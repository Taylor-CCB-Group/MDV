import { describe, expect, test, vi } from "vitest";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import ColorLegend from "@/react/components/legend/ColorLegend";

describe("ColorLegend", () => {
    test("renders categorical legend title and labels", () => {
        const { container } = render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "categorical",
                        label: "Cell type",
                        column: "cell_type",
                        items: [
                            {
                                color: "#ff0000",
                                name: "T cell",
                                value: "T cell",
                            },
                            {
                                color: "#00ff00",
                                name: "B cell",
                                value: "B cell",
                            },
                        ],
                    }}
                />
            </div>,
        );
        expect(screen.getByText("Cell type")).toBeTruthy();
        expect(screen.getByText("T cell")).toBeTruthy();
        expect(screen.getByText("B cell")).toBeTruthy();
        expect(container.querySelector(".legend-drag-handle")).toBeTruthy();
        expect(container.querySelector(".legend-title")?.className).toContain(
            "text-ellipsis",
        );
        const svg = container.querySelector("svg");
        expect(svg?.classList).toContain("overflow-hidden");
        expect(svg?.getAttribute("width")).toBe("100%");
    });

    test("marks hovered categorical row as focused", () => {
        render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "categorical",
                        label: "Cell type",
                        column: "cell_type",
                        items: [
                            {
                                color: "#ff0000",
                                name: "T cell",
                                value: "T cell",
                            },
                            {
                                color: "#00ff00",
                                name: "B cell",
                                value: "B cell",
                            },
                        ],
                    }}
                />
            </div>,
        );

        const hoveredRow = screen.getByText("T cell").closest("g");
        if (!hoveredRow) {
            throw new Error("Expected categorical legend rows to render");
        }
        fireEvent.mouseEnter(hoveredRow);

        expect(
            hoveredRow?.querySelector("rect[fill]")?.getAttribute("stroke"),
        ).toBe("currentColor");
    });

    test("calls categorical click handler with item value", () => {
        const onCategoricalItemClick = vi.fn();
        render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "categorical",
                        label: "Cell type",
                        column: "cell_type",
                        items: [
                            {
                                color: "#ff0000",
                                name: "T cell",
                                value: "T cell",
                            },
                        ],
                    }}
                    onCategoricalItemClick={onCategoricalItemClick}
                />
            </div>,
        );

        const row = screen.getByText("T cell").closest("g");
        if (!row) {
            throw new Error("Expected categorical legend row to render");
        }
        fireEvent.click(row);

        expect(onCategoricalItemClick).toHaveBeenCalledWith("T cell");
    });

    test("adds SVG title tooltip for truncated categorical row labels", () => {
        const longLabel =
            "CD8-positive cytotoxic T lymphocyte population cluster";
        const { container } = render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "categorical",
                        label: "Cell type",
                        column: "cell_type",
                        items: [
                            {
                                color: "#ff0000",
                                name: longLabel,
                                value: longLabel,
                            },
                        ],
                    }}
                />
            </div>,
        );

        const titles = container.querySelectorAll("title");
        expect(titles.length).toBeGreaterThan(0);
        expect(Array.from(titles).some((t) => t.textContent === longLabel)).toBe(
            true,
        );
        const rowText = container.querySelector(".legend-body svg text");
        expect(rowText?.textContent?.endsWith("…")).toBe(true);
    });

    test("marks active categorical row as selected", () => {
        render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "categorical",
                        label: "Cell type",
                        column: "cell_type",
                        items: [
                            {
                                color: "#ff0000",
                                name: "T cell",
                                value: "T cell",
                            },
                            {
                                color: "#00ff00",
                                name: "B cell",
                                value: "B cell",
                            },
                        ],
                    }}
                    activeCategoricalValue="T cell"
                />
            </div>,
        );

        const row = screen.getByText("T cell").closest("g");
        if (!row) {
            throw new Error("Expected categorical legend row to render");
        }
        expect(row?.querySelector("rect[fill]")?.getAttribute("stroke")).toBe(
            "currentColor",
        );
    });

    test("keeps active categorical row selected while another row is hovered", () => {
        render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "categorical",
                        label: "Cell type",
                        column: "cell_type",
                        items: [
                            {
                                color: "#ff0000",
                                name: "T cell",
                                value: "T cell",
                            },
                            {
                                color: "#00ff00",
                                name: "B cell",
                                value: "B cell",
                            },
                            {
                                color: "#0000ff",
                                name: "NK cell",
                                value: "NK cell",
                            },
                        ],
                    }}
                    activeCategoricalValue="T cell"
                />
            </div>,
        );

        const activeRow = screen.getByText("T cell").closest("g");
        const hoveredRow = screen.getByText("B cell").closest("g");
        if (!activeRow || !hoveredRow) {
            throw new Error("Expected categorical legend rows to render");
        }

        fireEvent.mouseEnter(hoveredRow);

        expect(activeRow.querySelector("rect[fill]")?.getAttribute("stroke")).toBe(
            "currentColor",
        );
        expect(hoveredRow.querySelector("rect[fill]")?.getAttribute("stroke")).toBe(
            "currentColor",
        );
    });

    test("renders continuous legend label and gradient bar", () => {
        const { container } = render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "continuous",
                        label: "Expression",
                        column: "expression",
                        colors: ["#000000", "#ffffff"],
                        range: [0, 1],
                    }}
                />
            </div>,
        );
        expect(screen.getByText("Expression")).toBeTruthy();
        expect(container.querySelector("linearGradient")).toBeTruthy();
        expect(container.querySelector("rect[fill^='url(#']")).toBeTruthy();
        expect(container.querySelector("svg")?.getAttribute("width")).toBe(
            "100%",
        );
        expect(container.querySelector("svg")?.getAttribute("height")).toBe(
            "100%",
        );
    });

    test("renders continuous brush over the gradient bar when interactive", async () => {
        const onContinuousRangeChange = vi.fn();
        const { container } = render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "continuous",
                        label: "Expression",
                        column: "expression",
                        colors: ["#000000", "#ffffff"],
                        range: [0, 100],
                    }}
                    onContinuousRangeChange={onContinuousRangeChange}
                />
            </div>,
        );

        await waitFor(() => {
            expect(container.querySelector(".brush .overlay")).toBeTruthy();
        });
        const overlay = container.querySelector(".brush .overlay");
        expect(overlay?.getAttribute("x")).toBe("28");
        expect(overlay?.getAttribute("width")).toBe("124");
    });

    test("renders active continuous range selection", async () => {
        const { container } = render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "continuous",
                        label: "Expression",
                        column: "expression",
                        colors: ["#000000", "#ffffff"],
                        range: [0, 100],
                    }}
                    activeContinuousRange={[25, 75]}
                />
            </div>,
        );

        await waitFor(() => {
            expect(container.querySelector(".brush .selection")).toBeTruthy();
        });
        expect(container.querySelector(".brush .selection")?.getAttribute("x")).toBe(
            "59",
        );
        expect(
            container.querySelector(".brush .selection")?.getAttribute("width"),
        ).toBe("62");
    });
});
