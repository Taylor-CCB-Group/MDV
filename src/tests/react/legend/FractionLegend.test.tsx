import { render, screen, waitFor } from "@testing-library/react";
import type { ReactNode } from "react";
import { describe, expect, test } from "vitest";
import FractionLegend from "@/react/components/legend/FractionLegend";
import LegendFrame from "@/react/components/legend/LegendFrame";

/** The legend in the frame the chart draws it in, floating or in the strip above the plot. */
function Frame({ children, inBand = false }: { children: ReactNode; inBand?: boolean }) {
    return (
        <LegendFrame
            inBand={inBand}
            container={document.createElement("div")}
            start={{ left: 0, top: 0 }}
            minimized={false}
            summary={{ label: "fraction" }}
        >
            {children}
        </LegendFrame>
    );
}

function getWrapper(container: HTMLElement) {
    const wrapper = container.querySelector(".legend-container");
    if (!(wrapper instanceof HTMLElement)) {
        throw new Error("Expected fraction legend frame element");
    }
    return wrapper;
}

describe("FractionLegend", () => {
    test("renders title, circles, labels, and sizes wrapper", async () => {
        const { container } = render(
            <Frame>
                <FractionLegend
                    spec={{
                        label: "fraction",
                        width: 100,
                        items: [
                            { key: "20", label: "20", radius: 8 },
                            { key: "40", label: "40", radius: 12 },
                        ],
                    }}
                />
            </Frame>,
        );

        expect(screen.getByText("fraction")).toBeTruthy();
        expect(screen.getByText("20")).toBeTruthy();
        expect(screen.getByText("40")).toBeTruthy();
        expect(container.querySelectorAll("circle")).toHaveLength(2);
        expect(container.querySelector(".legend-drag-handle")).toBeTruthy();

        await waitFor(() => expect(getWrapper(container).style.width).toBe("100px"));
    });

    test("auto-sizes wrapper to tall content without applying the old scroll cap", async () => {
        const { container } = render(
            <Frame>
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
            </Frame>,
        );

        const wrapper = getWrapper(container);
        await waitFor(() => expect(wrapper.style.height).toBe("285px"));
    });

    test("caps auto height when maxHeight is provided", async () => {
        const { container } = render(
            <Frame>
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
            </Frame>,
        );

        const wrapper = getWrapper(container);
        await waitFor(() => {
            expect(wrapper.style.width).toBe("120px");
            expect(wrapper.style.height).toBe("80px");
        });
    });

    test("lays the circles out in a row in the strip, leaving the width to the strip", async () => {
        const { container } = render(
            <Frame inBand>
                <FractionLegend
                    spec={{
                        label: "fraction",
                        width: 100,
                        items: [
                            { key: "20", label: "20", radius: 8 },
                            { key: "40", label: "40", radius: 12 },
                        ],
                    }}
                />
            </Frame>,
        );

        const circles = [...container.querySelectorAll("circle")];
        expect(new Set(circles.map((circle) => circle.getAttribute("cy"))).size).toBe(1);
        const wrapper = getWrapper(container);
        await waitFor(() => expect(wrapper.style.height).not.toBe(""));
        expect(wrapper.style.width).toBe("");
    });
});
