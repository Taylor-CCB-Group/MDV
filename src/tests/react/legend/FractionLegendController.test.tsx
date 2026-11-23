import { act, type ReactElement, StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import FractionLegendController, {
    type FractionLegendSettings,
} from "@/react/components/legend/FractionLegendController";
import type { LegendPlacement } from "@/react/components/legend/LegendWrapper";
import type { FractionLegendSpec } from "@/react/legend/fraction_legend/types";

// The real portal adds app providers the legend doesn't need.
vi.mock("@/react/react_utils", () => ({
    createMdvPortal: (component: ReactElement, container: HTMLElement) => {
        const root = createRoot(container);
        act(() => root.render(<StrictMode>{component}</StrictMode>));
        return root;
    },
}));

const SPEC: FractionLegendSpec = {
    label: "% expressing",
    items: [
        { key: "20", label: "20", radius: 4 },
        { key: "40", label: "40", radius: 6 },
    ],
};

/** The parts of a chart the controller uses; a real BaseChart takes minutes to load in tests. */
function createChart() {
    const contentDiv = document.createElement("div");
    const band = document.createElement("div");
    document.body.append(contentDiv, band);
    const listeners = new Set<() => void>();
    let onTop = false;
    const chart = {
        contentDiv,
        getLegendPlacement: (): LegendPlacement => ({
            onTop,
            band,
            setLayout: () => {},
            onBandChanged: () => {},
        }),
        subscribeLegendLayout: (listener: () => void) => {
            listeners.add(listener);
            return () => listeners.delete(listener);
        },
    };
    return {
        chart: chart as unknown as BaseChart<BaseConfig>,
        contentDiv,
        band,
        listeners,
        /** What `BaseChart.setLegendLayout` does, as far as the controller can tell. */
        setOnTop(value: boolean) {
            onTop = value;
            act(() => {
                for (const listener of listeners) listener();
            });
        },
    };
}

let controllers: FractionLegendController[] = [];

function setUp(settings: FractionLegendSettings) {
    const parts = createChart();
    const controller = new FractionLegendController(parts.chart, () => settings);
    controllers.push(controller);
    return { ...parts, controller, settings };
}

function frameIn(container: HTMLElement) {
    return container.querySelector<HTMLElement>(".legend-container");
}

beforeEach(() => {
    vi.spyOn(HTMLElement.prototype, "offsetLeft", "get").mockImplementation(function (this: HTMLElement) {
        return Number.parseFloat(this.style.left) || 0;
    });
    vi.spyOn(HTMLElement.prototype, "offsetTop", "get").mockImplementation(function (this: HTMLElement) {
        return Number.parseFloat(this.style.top) || 0;
    });
});

afterEach(async () => {
    await act(async () => {
        for (const controller of controllers) controller.remove();
    });
    controllers = [];
    document.body.innerHTML = "";
    vi.restoreAllMocks();
});

describe("FractionLegendController", () => {
    test("starts at the default position and saves it, so it stays put", () => {
        const { controller, contentDiv, settings } = setUp({ display: true });
        act(() => controller.show(() => SPEC, [12, 8]));

        const frame = frameIn(contentDiv);
        expect([frame?.style.left, frame?.style.top]).toEqual(["12px", "8px"]);
        expect(settings.position).toEqual([12, 8]);
    });

    test("prefers the saved position, then an older view's pos, over the default", () => {
        const saved = setUp({ display: true, position: [30, 40], pos: [1, 1] });
        act(() => saved.controller.show(() => SPEC, [12, 8]));
        expect(frameIn(saved.contentDiv)?.style.left).toBe("30px");

        const older = setUp({ display: true, pos: [50, 60] });
        act(() => older.controller.show(() => SPEC, [12, 8]));
        expect(frameIn(older.contentDiv)?.style.left).toBe("50px");
    });

    test("gives the builder the position it will start at", () => {
        const { controller } = setUp({ display: true, position: [5, 70] });
        const build = vi.fn(() => SPEC);
        act(() => controller.show(build, [0, 0]));
        expect(build).toHaveBeenCalledWith([5, 70]);
    });

    test("hides when switched off or when there's nothing to show", async () => {
        const { controller, contentDiv, settings } = setUp({ display: true });
        act(() => controller.show(() => SPEC, [0, 0]));
        expect(frameIn(contentDiv)).not.toBeNull();

        await act(async () => controller.show(null, [0, 0]));
        expect(frameIn(contentDiv)).toBeNull();

        settings.display = false;
        await act(async () => controller.show(() => SPEC, [0, 0]));
        expect(frameIn(contentDiv)).toBeNull();
    });

    test("saves where it was dragged", () => {
        const { controller, contentDiv, settings } = setUp({ display: true });
        act(() => controller.show(() => SPEC, [0, 0]));
        const frame = frameIn(contentDiv);
        if (!frame) throw new Error("expected the legend");
        frame.style.left = "90px";
        frame.style.top = "45px";

        controller.savePosition();

        expect(settings.position).toEqual([90, 45]);
    });

    test("minimises on its own button and saves that", () => {
        const { controller, contentDiv, settings } = setUp({ display: true });
        act(() => controller.show(() => SPEC, [0, 0]));

        const minimise = contentDiv.querySelector<HTMLButtonElement>('[aria-label="Minimise legend"]');
        act(() => minimise?.click());

        expect(settings.minimized).toBe(true);
        expect(frameIn(contentDiv)?.classList).toContain("legend-minimized");
        expect(contentDiv.textContent).toContain("% expressing");
    });

    test("moves into the strip above the plot and back with the chart's other legends", () => {
        const { controller, contentDiv, band, setOnTop } = setUp({ display: true, position: [30, 40] });
        act(() => controller.show(() => SPEC, [0, 0]));

        setOnTop(true);
        expect(frameIn(band)?.classList).toContain("legend-in-band");
        expect(frameIn(contentDiv)).toBeNull();

        setOnTop(false);
        expect(frameIn(contentDiv)?.style.left).toBe("30px");
    });

    test("stops following the chart's layout once removed", async () => {
        const { controller, listeners } = setUp({ display: true });
        expect(listeners.size).toBe(1);
        await act(async () => controller.remove());
        expect(listeners.size).toBe(0);
    });
});
