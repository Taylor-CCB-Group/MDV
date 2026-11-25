import { act, type ReactElement, StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import { useLegendPreferredSize } from "@/react/components/legend/LegendFrame";
import LegendWrapper, {
    type LegendPlacement,
    type LegendWrapperInteractions,
} from "@/react/components/legend/LegendWrapper";

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

/** Stands in for ColorLegend: asks the frame for its default size. */
function FakeLegend() {
    useLegendPreferredSize(DEFAULT_SIZE);
    return (
        <div>
            <div className="legend-drag-handle">Legend</div>
            <div className="legend-body">Items</div>
        </div>
    );
}

/** jsdom has no layout, so read sizes and positions from inline styles and data attributes. */
function px(value: string) {
    return Number.parseFloat(value) || 0;
}

function createChart(content: { width: number; height: number }, placement?: Partial<LegendPlacement>) {
    const contentDiv = document.createElement("div");
    contentDiv.dataset.testWidth = String(content.width);
    contentDiv.dataset.testHeight = String(content.height);
    document.body.append(contentDiv);
    const band = document.createElement("div");
    document.body.append(band);
    const chart = {
        contentDiv,
        __doc__: document,
        ...(placement && {
            getLegendPlacement: (): LegendPlacement => ({
                onTop: false,
                band,
                setLayout: () => {},
                onBandChanged: () => {},
                ...placement,
            }),
        }),
    };
    return {
        contentDiv,
        band,
        // LegendWrapper only reads contentDiv and getLegendPlacement from the chart. A real
        // BaseChart would pull in the whole chart code, which takes minutes to load.
        wrapper: new LegendWrapper<string>(chart as unknown as BaseChart<BaseConfig>),
    };
}

let wrappers: LegendWrapper<string>[] = [];

function renderLegend(
    content: { width: number; height: number },
    position: { left: number; top?: number; bottom?: number; size?: [number, number] },
    resizable = true,
    interactions: LegendWrapperInteractions = {},
    placement?: Partial<LegendPlacement>,
) {
    const { wrapper, contentDiv, band } = createChart(content, placement);
    wrappers.push(wrapper);
    act(() => {
        wrapper.render("spec", position, FakeLegend, {
            dragHandle: ".legend-drag-handle",
            resizable,
            ...interactions,
        });
    });
    const el = wrapper.getWrapperElement()?.querySelector(".legend-container");
    if (!(el instanceof HTMLElement)) throw new Error("legend frame was not drawn");
    return { wrapper, el, contentDiv, band };
}

/** Presses on `target`, moves the pointer by (dx, dy) and lets go. */
function drag(target: Element, dx: number, dy: number) {
    act(() => {
        target.dispatchEvent(new MouseEvent("pointerdown", { bubbles: true, button: 0, clientX: 100, clientY: 100 }));
        document.dispatchEvent(new MouseEvent("pointermove", { clientX: 100 + dx, clientY: 100 + dy }));
        document.dispatchEvent(new MouseEvent("pointerup"));
    });
}

/** A pointer event from a given pointer; jsdom's mouse events have no pointerId. */
function pointer(type: string, pointerId: number, init: MouseEventInit = {}) {
    const event = new MouseEvent(type, { bubbles: true, button: 0, ...init });
    Object.defineProperty(event, "pointerId", { value: pointerId });
    return event;
}

/** jsdom has no ResizeObserver; this one lets a test say when the chart has resized. */
class FakeResizeObserver {
    static instances: FakeResizeObserver[] = [];
    targets: Element[] = [];
    disconnected = false;
    callback: ResizeObserverCallback;
    constructor(callback: ResizeObserverCallback) {
        this.callback = callback;
        FakeResizeObserver.instances.push(this);
    }
    observe(target: Element) {
        this.targets.push(target);
    }
    unobserve() {}
    disconnect() {
        this.disconnected = true;
    }
    trigger() {
        this.callback([], this);
    }
}

/** Changes the chart's size and tells the legend's observer about it. */
function resizeChart(contentDiv: HTMLElement, width: number, height: number) {
    contentDiv.dataset.testWidth = String(width);
    contentDiv.dataset.testHeight = String(height);
    act(() => {
        for (const observer of FakeResizeObserver.instances) {
            if (!observer.disconnected && observer.targets.includes(contentDiv)) observer.trigger();
        }
    });
}

beforeEach(() => {
    FakeResizeObserver.instances = [];
    vi.stubGlobal("ResizeObserver", FakeResizeObserver);
    vi.spyOn(HTMLElement.prototype, "offsetLeft", "get").mockImplementation(function (this: HTMLElement) {
        return px(this.style.left);
    });
    vi.spyOn(HTMLElement.prototype, "offsetTop", "get").mockImplementation(function (this: HTMLElement) {
        if (this.style.top || !this.style.bottom) return px(this.style.top);
        // placed from the bottom of the chart it's in
        const chart = this.closest<HTMLElement>("[data-test-height]");
        return Number(chart?.dataset.testHeight ?? 0) - px(this.style.bottom) - px(this.style.height);
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
    vi.unstubAllGlobals();
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

        drag(corner, -40, -60);

        expect(wrapper.getUserSize()).toEqual([140, 190]);
        expect([el.style.width, el.style.height]).toEqual(["140px", "190px"]);
    });
});

describe("LegendWrapper resized while redrawn", () => {
    test("ignores a resize that finishes on a legend that's since been redrawn", () => {
        const { el, wrapper } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 });
        const corner = el.querySelector(".resizer-se");
        if (!corner) throw new Error("expected a resize handle");

        act(() => {
            corner.dispatchEvent(pointer("pointerdown", 1, { clientX: 100, clientY: 100 }));
            document.dispatchEvent(pointer("pointermove", 1, { clientX: 60, clientY: 40 }));
            // redrawn mid-resize; the old legend is only taken down after this
            wrapper.render("spec", { left: 10, top: 10 }, FakeLegend, {
                dragHandle: ".legend-drag-handle",
                resizable: true,
            });
            document.dispatchEvent(pointer("pointerup", 1));
        });

        expect(wrapper.getUserSize()).toBeNull();
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

    test("fits the legend once the chart gets a size", () => {
        const { el, contentDiv } = renderLegend({ width: 0, height: 0 }, { left: 350, top: 280 });
        resizeChart(contentDiv, 400, 300);
        expect([el.style.left, el.style.top]).toEqual(["220px", "50px"]);
    });

    test("moves the legend back inside when the chart shrinks", () => {
        // fits at first: 350 + 180 <= 600 and 280 + 250 <= 600
        const { el, contentDiv } = renderLegend({ width: 600, height: 600 }, { left: 350, top: 280 });
        expect([el.style.left, el.style.top]).toEqual(["350px", "280px"]);
        resizeChart(contentDiv, 400, 300);
        expect([el.style.left, el.style.top]).toEqual(["220px", "50px"]);
    });

    test("starts a legend placed from the bottom at that distance from the chart's bottom", () => {
        const { el, wrapper } = renderLegend({ width: 400, height: 300 }, { left: 10, bottom: 10 });
        expect([el.style.left, el.style.bottom]).toEqual(["10px", "10px"]);
        // 300 - 250 - 10
        expect(wrapper.getFloatingPosition()).toEqual([10, 40]);
    });

    test("moves a legend placed from the bottom down from the top when it's too tall to fit", () => {
        const { el } = renderLegend({ width: 400, height: 200 }, { left: 10, bottom: 10 });
        expect([el.style.left, el.style.top]).toEqual(["10px", "0px"]);
    });

    test("stops watching the chart when the legend is removed", async () => {
        const { wrapper, contentDiv } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 });
        const watching = () =>
            FakeResizeObserver.instances.filter((o) => !o.disconnected && o.targets.includes(contentDiv));
        expect(watching()).toHaveLength(1);
        await act(async () => wrapper.unmount());
        expect(watching()).toHaveLength(0);
    });
});

describe("LegendWrapper dragging", () => {
    test("moves the legend by its handle and reports where it was left", () => {
        const onMoved = vi.fn();
        const { el, wrapper } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 }, true, { onMoved });
        const handle = el.querySelector(".legend-drag-handle");
        if (!handle) throw new Error("expected a drag handle");

        drag(handle, 50, 30);

        expect([el.style.left, el.style.top]).toEqual(["60px", "40px"]);
        expect(onMoved).toHaveBeenCalledWith([60, 40]);
        expect(wrapper.getFloatingPosition()).toEqual([60, 40]);
    });

    test("doesn't move when pressed outside the handle", () => {
        const { el } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 });
        const body = el.querySelector(".legend-body");
        if (!body) throw new Error("expected the legend body");

        drag(body, 50, 30);

        expect([el.style.left, el.style.top]).toEqual(["10px", "10px"]);
    });

    test("follows only the pointer that started the drag", () => {
        const onMoved = vi.fn();
        const { el } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 }, true, { onMoved });
        const handle = el.querySelector(".legend-drag-handle");
        if (!handle) throw new Error("expected a drag handle");

        act(() => {
            handle.dispatchEvent(pointer("pointerdown", 1, { clientX: 100, clientY: 100 }));
            document.dispatchEvent(pointer("pointermove", 2, { clientX: 300, clientY: 300 }));
            document.dispatchEvent(pointer("pointerup", 2));
        });
        expect([el.style.left, el.style.top]).toEqual(["10px", "10px"]);
        expect(onMoved).not.toHaveBeenCalled();

        act(() => {
            document.dispatchEvent(pointer("pointermove", 1, { clientX: 150, clientY: 130 }));
            document.dispatchEvent(pointer("pointerup", 1));
        });
        expect([el.style.left, el.style.top]).toEqual(["60px", "40px"]);
        expect(onMoved).toHaveBeenCalledWith([60, 40]);
    });

    test("keeps following the first pointer when a second presses the handle before it's released", () => {
        const onMoved = vi.fn();
        const { el } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 }, true, { onMoved });
        const handle = el.querySelector(".legend-drag-handle");
        if (!handle) throw new Error("expected a drag handle");

        act(() => {
            handle.dispatchEvent(pointer("pointerdown", 1, { clientX: 100, clientY: 100 }));
            handle.dispatchEvent(pointer("pointerdown", 2, { clientX: 200, clientY: 200 }));
            document.dispatchEvent(pointer("pointermove", 2, { clientX: 400, clientY: 400 }));
            document.dispatchEvent(pointer("pointermove", 1, { clientX: 150, clientY: 130 }));
        });
        expect([el.style.left, el.style.top]).toEqual(["60px", "40px"]);

        act(() => {
            document.dispatchEvent(pointer("pointerup", 2));
            document.dispatchEvent(pointer("pointermove", 1, { clientX: 160, clientY: 140 }));
            document.dispatchEvent(pointer("pointerup", 1));
        });
        expect([el.style.left, el.style.top]).toEqual(["70px", "50px"]);
        expect(onMoved).toHaveBeenCalledTimes(1);
        expect(onMoved).toHaveBeenCalledWith([70, 50]);
    });

    test("reports where the legend is shown when the chart shrinks during the drag", () => {
        const onMoved = vi.fn();
        const { el, wrapper, contentDiv } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 }, true, {
            onMoved,
        });
        const handle = el.querySelector(".legend-drag-handle");
        if (!handle) throw new Error("expected a drag handle");

        act(() => {
            handle.dispatchEvent(pointer("pointerdown", 1, { clientX: 100, clientY: 100 }));
            document.dispatchEvent(pointer("pointermove", 1, { clientX: 300, clientY: 300 }));
        });
        expect([el.style.left, el.style.top]).toEqual(["210px", "210px"]);
        // the legend (180 x 250) is moved back inside a chart now 400 x 400
        resizeChart(contentDiv, 400, 400);
        act(() => document.dispatchEvent(pointer("pointerup", 1)));

        expect([el.style.left, el.style.top]).toEqual(["210px", "150px"]);
        expect(onMoved).toHaveBeenCalledWith([210, 150]);
        expect(wrapper.getFloatingPosition()).toEqual([210, 150]);
    });

    test("drops a drag still going when the legend is redrawn", async () => {
        const onMoved = vi.fn();
        const { el, wrapper } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 }, true, { onMoved });
        const handle = el.querySelector(".legend-drag-handle");
        if (!handle) throw new Error("expected a drag handle");
        act(() => {
            handle.dispatchEvent(pointer("pointerdown", 1, { clientX: 100, clientY: 100 }));
            document.dispatchEvent(pointer("pointermove", 1, { clientX: 150, clientY: 130 }));
        });

        await act(async () => wrapper.unmount());
        act(() => document.dispatchEvent(pointer("pointerup", 1)));

        expect(onMoved).not.toHaveBeenCalled();
    });

    test("can be dragged in a chart popped out into another window", () => {
        const frame = document.createElement("iframe");
        document.body.append(frame);
        const popout = frame.contentDocument;
        if (!popout) throw new Error("expected the other window's document");
        const { wrapper, contentDiv } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 });
        popout.body.append(contentDiv);
        // redrawn there, so it's made of the other window's elements
        // (at 0, 0: the stand-in offsets above are this window's, so they don't apply to them)
        act(() => wrapper.render("spec", { left: 0, top: 0 }, FakeLegend, { dragHandle: ".legend-drag-handle" }));
        const el = wrapper.getWrapperElement()?.querySelector<HTMLElement>(".legend-container");
        const handle = el?.querySelector(".legend-drag-handle");
        if (!el || !handle) throw new Error("expected a drag handle");
        expect(handle instanceof Element).toBe(false);

        act(() => {
            handle.dispatchEvent(pointer("pointerdown", 1, { clientX: 100, clientY: 100 }));
            popout.dispatchEvent(pointer("pointermove", 1, { clientX: 150, clientY: 130 }));
            popout.dispatchEvent(pointer("pointerup", 1));
        });

        expect([el.style.left, el.style.top]).toEqual(["50px", "30px"]);
    });

    test("keeps a dragged legend inside the chart", () => {
        const { el } = renderLegend({ width: 400, height: 300 }, { left: 10, top: 10 });
        const handle = el.querySelector(".legend-drag-handle");
        if (!handle) throw new Error("expected a drag handle");

        drag(handle, 1000, 1000);

        expect([el.style.left, el.style.top]).toEqual(["220px", "50px"]);
    });
});

describe("LegendWrapper minimised", () => {
    test("opens only from its + button, and can be dragged by the rest", () => {
        const onMinimizedChange = vi.fn();
        const { el } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 }, true, {
            minimized: true,
            onMinimizedChange,
            summary: { label: "Cell type", swatches: ["#f00"] },
        });
        expect(el.classList).toContain("legend-minimized");
        expect(el.textContent).toContain("Cell type");

        const pill = el.querySelector(".legend-summary");
        if (!pill) throw new Error("expected the minimised legend");
        act(() => pill.dispatchEvent(new MouseEvent("click", { bubbles: true })));
        expect(onMinimizedChange).not.toHaveBeenCalled();

        drag(pill, 40, 20);
        expect([el.style.left, el.style.top]).toEqual(["50px", "30px"]);

        const open = el.querySelector<HTMLButtonElement>(".legend-summary-open");
        act(() => open?.click());
        expect(onMinimizedChange).toHaveBeenCalledWith(false);
    });

    test("minimises from the legend's own button", () => {
        const onMinimizedChange = vi.fn();
        const { el } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 }, true, { onMinimizedChange });
        const minimise = el.querySelector<HTMLButtonElement>('[aria-label="Minimise legend"]');
        act(() => minimise?.click());
        expect(onMinimizedChange).toHaveBeenCalledWith(true);
    });
});

describe("LegendWrapper keyboard focus", () => {
    /** A legend whose minimise and + buttons redraw it, as the charts' do. */
    function renderRedrawnLegend() {
        let redraw = (_minimized: boolean) => {};
        const interactions: LegendWrapperInteractions = {
            summary: { label: "Cell type" },
            onMinimizedChange: (minimized) => redraw(minimized),
        };
        const { wrapper } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 }, true, interactions);
        redraw = (minimized) =>
            wrapper.render("spec", { left: 10, top: 10 }, FakeLegend, { ...interactions, minimized });
        const find = (selector: string) => wrapper.getWrapperElement()?.querySelector<HTMLButtonElement>(selector);
        return { wrapper, redraw, find };
    }

    test("stays on the legend's button as it's minimised and opened from the keyboard", () => {
        const { find } = renderRedrawnLegend();

        find('[aria-label="Minimise legend"]')?.focus();
        act(() => find('[aria-label="Minimise legend"]')?.click());
        expect(document.activeElement).toBe(find(".legend-summary-open"));

        act(() => find(".legend-summary-open")?.click());
        expect(document.activeElement).toBe(find('[aria-label="Minimise legend"]'));
    });

    test("is left alone when the legend didn't have it, or is redrawn for another reason", () => {
        const { find, redraw } = renderRedrawnLegend();
        const elsewhere = document.createElement("button");
        document.body.append(elsewhere);
        elsewhere.focus();

        act(() => find('[aria-label="Minimise legend"]')?.click());
        expect(document.activeElement).toBe(elsewhere);

        act(() => redraw(false));
        expect(document.activeElement).toBe(elsewhere);
    });
});

describe("LegendWrapper in the strip above the plot", () => {
    test("draws the legend in the strip, sized by the strip, with no floating position", () => {
        const { el, band, wrapper } = renderLegend(
            { width: 600, height: 600 },
            { left: 10, top: 10 },
            true,
            {},
            {
                onTop: true,
            },
        );
        expect(band.contains(el)).toBe(true);
        expect(el.classList).toContain("legend-in-band");
        expect(el.style.width).toBe("");
        expect(el.style.height).toBe("96px");
        expect(el.querySelector(".resizer-se")).toBeNull();
        expect(el.querySelector('[aria-label="Put legends above the plot"]')).toBeNull();
        expect(wrapper.getFloatingPosition()).toBeNull();
    });

    test("uses the legend's own height in the strip, not the size it was resized to while floating", () => {
        const { el } = renderLegend(
            { width: 600, height: 600 },
            { left: 10, top: 10, size: [120, 40] },
            true,
            {},
            {
                onTop: true,
            },
        );
        // 250 asked for, capped at the strip's 96
        expect(el.style.height).toBe("96px");
    });

    test("offers to put the legends above the plot while floating", () => {
        const setLayout = vi.fn();
        const { el } = renderLegend({ width: 600, height: 600 }, { left: 10, top: 10 }, true, {}, { setLayout });
        const button = el.querySelector<HTMLButtonElement>('[aria-label="Put legends above the plot"]');
        act(() => button?.click());
        expect(setLayout).toHaveBeenCalledWith({ position: "top" });
    });
});
