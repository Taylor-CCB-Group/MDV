import { describe, expect, test, vi } from "vitest";
import WGLScatterPlot from "@/charts/WGLScatterPlot.js";

function makeStub(overrides: {
    x: string;
    y: string;
    columns: Record<string, { is_date?: boolean; date_unit?: string }>;
    axis: {
        x?: { size?: number };
        y?: { size?: number };
        x_log_scale?: boolean;
        y_log_scale?: boolean;
    };
}) {
    const stub = {
        x: overrides.x,
        y: overrides.y,
        dataStore: { columnIndex: overrides.columns },
        config: { axis: overrides.axis },
        _dateAxisCache: undefined as
            | {
                  x?: { log_scale?: boolean; size?: number } | null;
                  y?: { log_scale?: boolean; size?: number } | null;
              }
            | undefined,
        setAxisSize: vi.fn(function (
            this: { config: { axis: Record<string, { size?: number }> } },
            axis: "x" | "y",
            size: number,
        ) {
            this.config.axis[axis].size = size;
        }),
    };
    return stub;
}

describe("_applyDateAxisDefaults date↔numeric restore", () => {
    test("date-to-numeric switch restores log scale and original axis size", () => {
        const stub = makeStub({
            x: "qc_date",
            y: "n_genes",
            columns: {
                qc_date: { is_date: true },
                n_genes: {},
                score: {},
            },
            axis: {
                x: { size: 30 },
                y: { size: 45 },
                x_log_scale: true,
                y_log_scale: false,
            },
        });

        WGLScatterPlot.prototype._applyDateAxisDefaults.call(stub);

        expect(stub.config.axis.x_log_scale).toBe(false);
        expect(stub.config.axis.x?.size).toBe(40);
        expect(stub.setAxisSize).toHaveBeenCalledWith("x", 40);

        stub.x = "score";
        WGLScatterPlot.prototype._applyDateAxisDefaults.call(stub);

        expect(stub.config.axis.x_log_scale).toBe(true);
        expect(stub.config.axis.x?.size).toBe(30);
        expect(stub.setAxisSize).toHaveBeenCalledWith("x", 30);
    });
});
