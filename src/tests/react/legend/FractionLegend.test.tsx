import { render, screen, waitFor } from "@testing-library/react";
import { describe, expect, test, vi } from "vitest";
import FractionLegend from "@/react/components/legend/FractionLegend";

function getWrapper(container: HTMLElement) {
    const wrapper = container.firstElementChild;
    if (!(wrapper instanceof HTMLElement)) {
        throw new Error("Expected fraction legend wrapper element");
    }
    return wrapper;
}

describe("FractionLegend", () => {
    test("renders title, circles, labels, and sizes wrapper", async () => {
        const onLayoutReady = vi.fn();
        const { container } = render(
            <div className="legend-container">
                <FractionLegend
                    spec={{
                        label: "fraction",
                        width: 100,
                        items: [
                            { key: "20", label: "20", radius: 8 },
                            { key: "40", label: "40", radius: 12 },
                        ],
                    }}
                    onLayoutReady={onLayoutReady}
                />
            </div>,
        );

        expect(screen.getByText("fraction")).toBeTruthy();
        expect(screen.getByText("20")).toBeTruthy();
        expect(screen.getByText("40")).toBeTruthy();
        expect(container.querySelectorAll("circle")).toHaveLength(2);
        expect(container.querySelector(".legend-drag-handle")).toBeTruthy();

        await waitFor(() => expect(onLayoutReady).toHaveBeenCalled());
        expect(getWrapper(container).style.width).toBe("100px");
    });

    test("auto-sizes wrapper to tall content without applying the old scroll cap", async () => {
        const { container } = render(
            <div className="legend-container">
                <FractionLegend
                    spec={{
                        label: "fraction",
                        width: 100,
                        items: [
                            { key: "20", label: "20", radius: 18 },
                            { key: "40", label: "40", radius: 21 },
                            { key: "60", label: "60", radius: 23 },
                            { key: "80", label: "80", radius: 24 },
                            { key: "100", label: "100", radius: 25 },
                        ],
                    }}
                />
            </div>,
        );

        const wrapper = getWrapper(container);
        await waitFor(() => expect(wrapper.style.height).toBe("285px"));
    });

    test("caps auto height when maxHeight is provided", async () => {
        const { container } = render(
            <div className="legend-container">
                <FractionLegend
                    spec={{
                        label: "fraction",
                        width: 120,
                        maxHeight: 80,
                        items: [
                            { key: "20", label: "20", radius: 18 },
                            { key: "40", label: "40", radius: 21 },
                        ],
                    }}
                />
            </div>,
        );

        const wrapper = getWrapper(container);
        await waitFor(() => {
            expect(wrapper.style.width).toBe("120px");
            expect(wrapper.style.height).toBe("80px");
        });
    });
});
