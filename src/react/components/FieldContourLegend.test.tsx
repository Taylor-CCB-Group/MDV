import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, test, vi } from "vitest";
import FieldContourLegend, { type CategoryLegendItem, type FieldLegendItem } from "./FieldContourLegend";

const fields: FieldLegendItem[] = [{ name: "CD3", field: "cd3", color: [10, 20, 30] }];
const categories: CategoryLegendItem[] = [
    { id: "contour1", name: "T cell", colors: [[255, 0, 0]] },
    { id: "contour2", name: "B cell, NK", colors: [[0, 0, 0], [255, 255, 255]] },
];

describe("FieldContourLegend", () => {
    test("renders nothing without fields or categories", () => {
        const { container } = render(<FieldContourLegend fields={[]} categories={[]} />);
        expect(container.firstChild).toBeNull();
    });

    test("keeps the fields-only legend unchanged", () => {
        render(<FieldContourLegend fields={fields} />);
        expect(screen.getByText("Density Fields")).toBeTruthy();
        expect(screen.queryByText("Fields")).toBeNull();
    });

    test("groups categorical contours under the contour parameter name", () => {
        const { container } = render(
            <FieldContourLegend fields={fields} categories={categories} categoryLabel="Cell type" />,
        );
        expect(screen.getByText("Density")).toBeTruthy();
        expect(screen.getByText("Cell type")).toBeTruthy();
        expect(screen.getByText("Fields")).toBeTruthy();
        expect(screen.getByText("T cell", { selector: "text" })).toBeTruthy();
        expect(screen.getByText("CD3")).toBeTruthy();

        const swatches = [...container.querySelectorAll("rect")].map((rect) => rect.getAttribute("fill"));
        expect(swatches[0]).toBe("rgb(255, 0, 0)");
        // a contour over several categories shows its colour range as a gradient
        const gradient = container.querySelector("linearGradient");
        expect(gradient?.querySelectorAll("stop")).toHaveLength(2);
        expect(swatches[1]).toBe(`url(#${gradient?.id})`);
    });

    test("reports hovered categorical contours and fields", () => {
        const onContourHover = vi.fn();
        render(<FieldContourLegend fields={fields} categories={categories} onContourHover={onContourHover} />);

        fireEvent.mouseEnter(screen.getByText("T cell", { selector: "text" }));
        expect(onContourHover).toHaveBeenLastCalledWith({ type: "category", id: "contour1" });
        fireEvent.mouseLeave(screen.getByText("T cell", { selector: "text" }));
        expect(onContourHover).toHaveBeenLastCalledWith(null);

        fireEvent.mouseEnter(screen.getByText("CD3"));
        expect(onContourHover).toHaveBeenLastCalledWith({ type: "field", field: "cd3" });
    });
});
