import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import DotPlot from "@/charts/DotPlot.js";
import type CatColDimension from "@/datastore/CatColDimension";
import type DataStore from "@/datastore/DataStore";
import { createStore } from "../datastore/dimensionHarness";

// These tests only call DotPlot's own data methods. Stubbing the chart base classes, legends and
// loaders it imports keeps the file from compiling most of the app (minutes, not milliseconds).
vi.mock("@/charts/BaseChart", () => ({ default: { types: {} } }));
vi.mock("@/charts/SVGChart.js", () => ({ default: class {} }));
vi.mock("@/datastore/decorateColumnMethod", () => ({ loadColumnData: (method: unknown) => method }));
vi.mock("@/dataloaders/DataLoaderUtil", () => ({ loadColumn: vi.fn() }));
vi.mock("@/react/legend/color_legend/buildColorLegendSpec", () => ({ buildColorLegendSpec: vi.fn() }));
vi.mock("@/react/components/legend/FractionLegend", () => ({ default: () => null }));
vi.mock("@/react/components/legend/LegendWrapper", () => ({ default: class {} }));

type Averages = { mean_range: [number, number]; frac_range: [number, number] };

/**
 * Stands in for catColWorker.js: holds requests until the test answers them. Like the real
 * worker, answers come back in the order the requests were sent.
 */
class QueuedWorker {
    onmessage: ((event: { data: Averages }) => void) | null = null;
    onerror: ((event: { message: string }) => void) | null = null;
    requests: unknown[] = [];
    postMessage(message: unknown) {
        this.requests.push(message);
    }
    answer(data: Averages) {
        this.onmessage?.({ data });
    }
    fail() {
        this.onerror?.({ message: "worker error" });
    }
    terminate() {}
}

let worker: QueuedWorker;
let ds: DataStore;
beforeEach(() => {
    vi.stubGlobal(
        "Worker",
        vi.fn().mockImplementation(function MockWorker() {
            worker = new QueuedWorker();
            return worker;
        }),
    );
    ds = createStore(4, [
        { field: "cell_type", datatype: "text", data: ["a", "b", "a", "b"] },
        { field: "default_gene", datatype: "double", data: [0, 0.001, 0, 0.003] },
        { field: "selected_gene", datatype: "double", data: [2, 8.5, 0, 3] },
    ]);
});
afterEach(() => {
    vi.unstubAllGlobals();
});

// A dot plot with only what onDataFiltered and the colour range need, on a real CatColDimension.
function makeDotPlot() {
    const chart = Object.create(DotPlot.prototype);
    Object.assign(chart, {
        config: { param: ["cell_type"], threshold_value: 0 },
        dataStore: ds,
        dim: ds.getDimension("catcol_dimension"),
        x_scale: { domain: () => {} },
        fractionScale: { domain: () => {} },
        clusterRows: () => {},
        clusterColumns: () => {},
        // keep the real range caching, without building colours or legends
        setColorFunction: () => chart.getColorScaleRange(),
        drawChart: () => {},
    });
    return chart;
}

const defaultGenes: Averages = { mean_range: [0, 0.003], frac_range: [0, 50] };
const selectedGenes: Averages = { mean_range: [0, 8.5], frac_range: [0, 100] };

describe("DotPlot averages from the worker", () => {
    it("ignores an answer for fields that have since changed", () => {
        const chart = makeDotPlot();
        // the link's default genes, then the selected genes before the first answer arrives
        chart._setFields(["default_gene"]);
        chart._setFields(["selected_gene"]);
        worker.answer(defaultGenes);
        worker.answer(selectedGenes);

        expect(chart.data).toBe(selectedGenes);
        expect(chart.colorLegendRange).toEqual([0, 8.5]);
    });

    it("keeps the colour range when only the filter changes", () => {
        const chart = makeDotPlot();
        chart._setFields(["selected_gene"]);
        worker.answer(selectedGenes);
        chart.onDataFiltered();
        worker.answer({ mean_range: [0, 2], frac_range: [0, 50] });

        expect(chart.data.mean_range).toEqual([0, 2]);
        expect(chart.colorLegendRange).toEqual([0, 8.5]);
    });
});

describe("CatColDimension.getAverages", () => {
    it("gives each answer to the request that asked for it", () => {
        const dim = ds.getDimension("catcol_dimension") as CatColDimension;
        const first = vi.fn();
        const second = vi.fn();
        const third = vi.fn();
        dim.getAverages(first, ["cell_type", "default_gene"]);
        dim.getAverages(second, ["cell_type", "selected_gene"]);
        dim.getAverages(third, ["cell_type", "selected_gene"]);
        worker.answer(defaultGenes);
        // the second request failed in the worker, so the next answer is the third's
        worker.fail();
        worker.answer(selectedGenes);

        expect(first).toHaveBeenCalledExactlyOnceWith(defaultGenes);
        expect(second).not.toHaveBeenCalled();
        expect(third).toHaveBeenCalledExactlyOnceWith(selectedGenes);
    });
});
