import { describe, expect, test, vi } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import ColorLegend from "./ColorLegend";

describe("ColorLegend", () => {
    test("renders categorical legend title and labels", () => {
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
        expect(screen.getByText("Cell type")).toBeTruthy();
        expect(screen.getByText("T cell")).toBeTruthy();
        expect(screen.getByText("B cell")).toBeTruthy();
    });

    test("applies local hover styling for categorical rows", () => {
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
        const dimmedRow = screen.getByText("B cell").closest("g");
        if (!hoveredRow || !dimmedRow) {
            throw new Error("Expected categorical legend rows to render");
        }
        fireEvent.mouseEnter(hoveredRow);

        expect(hoveredRow?.querySelector("rect")?.getAttribute("stroke")).toBe(
            "currentColor",
        );
        expect(dimmedRow?.getAttribute("style")).toContain("opacity: 0.4");
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

    test("applies active categorical styling", () => {
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
        const dimmedRow = screen.getByText("B cell").closest("g");
        if (!row) {
            throw new Error("Expected categorical legend row to render");
        }
        expect(row?.querySelector("rect")?.getAttribute("stroke")).toBe(
            "currentColor",
        );
        expect(screen.getByText("T cell").getAttribute("style")).toContain(
            "font-weight: 700",
        );
        expect(dimmedRow?.getAttribute("style")).toContain("opacity: 0.4");
    });

    test("renders continuous legend label and gradient bar", () => {
        const { container } = render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "continuous",
                        label: "Expression",
                        colors: ["#000000", "#ffffff"],
                        range: [0, 1],
                    }}
                />
            </div>,
        );
        expect(screen.getByText("Expression")).toBeTruthy();
        expect(container.querySelector("linearGradient")).toBeTruthy();
        expect(container.querySelector("rect[fill^='url(#']")).toBeTruthy();
    });
});
