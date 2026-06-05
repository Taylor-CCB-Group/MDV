import { describe, expect, test } from "vitest";
import { render, screen } from "@testing-library/react";
import ColorLegend from "./ColorLegend";

describe("ColorLegend", () => {
    test("renders categorical legend title and labels", () => {
        render(
            <div className="legend-container">
                <ColorLegend
                    spec={{
                        kind: "categorical",
                        label: "Cell type",
                        items: [
                            { color: "#ff0000", name: "T cell" },
                            { color: "#00ff00", name: "B cell" },
                        ],
                    }}
                />
            </div>,
        );
        expect(screen.getByText("Cell type")).toBeTruthy();
        expect(screen.getByText("T cell")).toBeTruthy();
        expect(screen.getByText("B cell")).toBeTruthy();
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
