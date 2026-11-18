import { afterEach, describe, expect, test, vi } from "vitest";
import type BaseChart from "@/charts/BaseChart";
import ViewManager from "@/charts/ViewManager";
import { convertLegacyChartsSetting } from "@/lib/browserSettings";

afterEach(() => {
    convertLegacyChartsSetting.set("off");
});

function viewManagerWithCharts(types: string[]) {
    const charts = Object.fromEntries(
        types.map((type, i) => [`c${i}`, { chart: { config: { id: `c${i}`, type }, dataStore: {} } }]),
    );
    (window as unknown as { mdv: unknown }).mdv = { chartManager: { charts } };
    const viewManager = new ViewManager();
    const convert = vi.spyOn(viewManager, "convertLegacyChart").mockResolvedValue(true);
    const report = vi.spyOn(viewManager, "reportLegacyConversion").mockImplementation(() => {});
    return { convert, report };
}

describe("turning on Convert legacy charts", () => {
    test("converts the open charts that have an equivalent, sharing one alert", async () => {
        const { convert, report } = viewManagerWithCharts(["table_chart", "DeckContourScatter", "wgl_scatter_plot"]);
        convertLegacyChartsSetting.set("on");
        await vi.waitFor(() => expect(report).toHaveBeenCalledTimes(2));
        const converted = convert.mock.calls.map(([chart]) => (chart as BaseChart<any>).config.id);
        expect(converted).toEqual(["c0", "c2"]);
        expect(convert.mock.calls.every(([chart, options]) => options?.announce === false && !!chart)).toBe(true);
    });

    test("does nothing when turned off", async () => {
        convertLegacyChartsSetting.set("on");
        const { convert } = viewManagerWithCharts(["table_chart"]);
        convertLegacyChartsSetting.set("off");
        await new Promise((resolve) => setTimeout(resolve, 0));
        expect(convert).not.toHaveBeenCalled();
    });

    test("counts only the charts that were converted", async () => {
        const { convert, report } = viewManagerWithCharts(["table_chart", "table_chart"]);
        convert.mockResolvedValueOnce(false);
        convertLegacyChartsSetting.set("on");
        await vi.waitFor(() => expect(convert).toHaveBeenCalledTimes(2));
        // the second conversion's report would follow its own await
        await new Promise((resolve) => setTimeout(resolve, 0));
        expect(report).toHaveBeenCalledTimes(1);
    });
});
