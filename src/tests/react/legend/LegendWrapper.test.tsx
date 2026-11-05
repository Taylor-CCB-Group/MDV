import { act, type ReactElement, StrictMode, useLayoutEffect, useRef } from "react";
import { createRoot } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import LegendWrapper, { type LegendWrapperComponentProps } from "@/react/components/legend/LegendWrapper";

// The real portal adds app providers the wrapper doesn't need; keep StrictMode,
// since it runs layout effects twice, which is what the saved size has to survive.
vi.mock("@/react/react_utils", () => ({
    createMdvPortal: (component: ReactElement, container: HTMLElement) => {
        const root = createRoot(container);
        act(() => root.render(<StrictMode>{component}</StrictMode>));
        return root;
    },
}));

const DEFAULT_SIZE = { width: 180, height: 250 };

/** Stands in for ColorLegend: sets the wrapper's default size, then reports layout ready. */
function FakeLegend({ onLayoutReady }: LegendWrapperComponentProps<string>) {
    const ref = useRef<HTMLDivElement>(null);
    useLayoutEffect(() => {
        const el = ref.current?.parentElement;
        if (!el) return;
        el.style.width = `${DEFAULT_SIZE.width}px`;
        el.style.height = `${DEFAULT_SIZE.height}px`;
        onLayoutReady?.();
    }, [onLayoutReady]);
    return (
        <div ref={ref}>
            <div className="legend-drag-handle">Legend</div>
        </div>
    );
}

/** jsdom has no layout, so read sizes and positions from inline styles and data attributes. */
function px(value: string) {
    return Number.parseFloat(value) || 0;
}

function createChart(content: { width: number; height: number }) {
    const contentDiv = document.createElement("div");
    contentDiv.dataset.testWidth = String(content.width);
    contentDiv.dataset.testHeight = String(content.height);
    document.body.append(contentDiv);
    const chart = { contentDiv, __doc__: document };
    return {
        contentDiv,
        // LegendWrapper only reads contentDiv and __doc__ from the chart. A real
        // BaseChart would pull in the whole chart code, which takes minutes to load.
        wrapper: new LegendWrapper<string>(chart as unknown as BaseChart<BaseConfig>),
    };
}

let wrappers: LegendWrapper<string>[] = [];

function renderLegend(
    content: { width: number; height: number },
    position: { left: number; top: number; size?: [number, number] },
    resizable = true,
) {
    const { wrapper, contentDiv } = createChart(content);
    wrappers.push(wrapper);
    act(() => {
        wrapper.render("spec", position, FakeLegend, {
            dragHandle: ".legend-drag-handle",
            resizable,
        });
    });
    const el = wrapper.getWrapperElement();
    if (!el) throw new Error("legend wrapper was not created");
    return { wrapper, el, contentDiv };
}

beforeEach(() => {
    vi.spyOn(HTMLElement.prototype, "offsetLeft", "get").mockImplementation(function (this: HTMLElement) {
        return px(this.style.left);
    });
    vi.spyOn(HTMLElement.prototype, "offsetTop", "get").mockImplementation(function (this: HTMLElement) {
        return px(this.style.top);
    });
    vi.spyOn(HTMLElement.prototype, "offsetWidth", "get").mockImplementation(function (this: HTMLElement) {
        return px(this.style.width);
    });
    vi.spyOn(HTMLElement.prototype, "offsetHeight", "get").mockImplementation(function (this: HTMLElement) {
        return px(this.style.height);
    });
    vi.spyOn(Element.prototype, "clientWidth", "get").mockImplementation(function (this: Element) {
        return this instanceof HTMLElement ? Number(this.dataset.testWidth ?? 0) : 0;
    });
    vi.spyOn(Element.prototype, "clientHeight", "get").mockImplementation(function (this: Element) {
        return this instanceof HTMLElement ? Number(this.dataset.testHeight ?? 0) : 0;
    });
});

afterEach(() => {
    act(() => {
        for (const wrapper of wrappers) wrapper.unmount();
    });
    wrappers = [];
    document.body.innerHTML = "";
    vi.restoreAllMocks();
});

describe("LegendWrapper saved size", () => {
    test("applies a saved size after the legend sets its default size", () => {
        const { el, wrapper } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10, size: [120, 90] });
        expect([el.style.width, el.style.height]).toEqual(["120px", "90px"]);
        expect(wrapper.getUserSize()).toEqual([120, 90]);
    });

    test("keeps the default size when there is no saved size", () => {
        const { el, wrapper } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 });
        expect([el.style.width, el.style.height]).toEqual(["180px", "250px"]);
        expect(wrapper.getUserSize()).toBeNull();
    });

    test("ignores a saved size when the legend can't be resized", () => {
        const { el } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10, size: [120, 90] }, false);
        expect([el.style.width, el.style.height]).toEqual(["180px", "250px"]);
    });

    test("records the size the user resizes the legend to", () => {
        const { el, wrapper } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 });
        const corner = el.querySelector(".resizer-se");
        if (!corner) throw new Error("expected a resize handle");

        corner.dispatchEvent(new MouseEvent("mousedown", { bubbles: true, clientX: 100, clientY: 100 }));
        document.documentElement.dispatchEvent(new MouseEvent("mousemove", { clientX: 60, clientY: 40 }));
        document.documentElement.dispatchEvent(new MouseEvent("mouseup"));

        expect(wrapper.getUserSize()).toEqual([140, 190]);
    });
});

describe("LegendWrapper keeps the legend inside the chart", () => {
    test("moves a legend that would hang over the bottom-right edge back inside", () => {
        // 180x250 legend saved at the bottom right of a chart it no longer fits at
        const { el } = renderLegend({ width: 400, height: 300 }, { left: 350, top: 280 });
        expect([el.style.left, el.style.top]).toEqual(["220px", "50px"]);
    });

    test("uses the saved size when deciding whether it fits", () => {
        const { el } = renderLegend({ width: 400, height: 300 }, { left: 250, top: 200, size: [120, 90] });
        // 250 + 120 <= 400 and 200 + 90 <= 300, so it stays where it was saved
        expect([el.style.left, el.style.top]).toEqual(["250px", "200px"]);
    });

    test("pins a legend larger than the chart to the top left", () => {
        const { el } = renderLegend({ width: 100, height: 100 }, { left: 40, top: 40 });
        expect([el.style.left, el.style.top]).toEqual(["0px", "0px"]);
    });

    test("leaves the position alone while the chart has no size yet", () => {
        const { el } = renderLegend({ width: 0, height: 0 }, { left: 350, top: 280 });
        expect([el.style.left, el.style.top]).toEqual(["350px", "280px"]);
    });
});
