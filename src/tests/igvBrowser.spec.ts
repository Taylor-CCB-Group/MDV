import { describe, expect, it } from "vitest";
import BaseChart from "@/charts/BaseChart";
import "@/charts/GenomeBrowsers/IGVBrowser/IGVBrowser";
import {
    applyViewMargins,
    getLocationFieldsFromGenome,
    locationFromFieldValues,
    getStructuralVariantStyle,
} from "@/charts/GenomeBrowsers/IGVBrowser/igvUtils";

describe("igv utils", () => {
    it("gets location fields from genomic_location and svs", () => {
        expect(
            getLocationFieldsFromGenome({ genomic_location: { columns: { chr: "chr", start: "start", end: "end" } } }),
        ).toEqual(["chr", "start", "end"]);

        expect(
            getLocationFieldsFromGenome({ svs: { sv_columns: { chr1: "c1", pos1: "p1", pos2: "p2", chr2: "c2" } } }),
        ).toEqual(["c1", "p1", "p2", "c2"]);
    });

    it("applies margins consistently", () => {
        expect(applyViewMargins({ chr: "chr1", start: 100, end: 200 }, { type: "absolute", value: 25 })).toEqual({
            chr: "chr1",
            start: 75,
            end: 225,
        });

        expect(applyViewMargins({ chr: "chr1", start: 100, end: 200 }, { type: "percentage", value: 10 })).toEqual({
            chr: "chr1",
            start: 90,
            end: 210,
        });

        expect(applyViewMargins({ chr: "chr1", start: 100, end: 200 }, { type: "fixed_length", value: 40 })).toEqual({
            chr: "chr1",
            start: 60,
            end: 240,
        });
    });

    it("returns paired loci for cross-chromosome svs", () => {
        expect(
            locationFromFieldValues({ chr: "chr1", start: 1000, end: 500000, chr2: "chr2" }, true),
        ).toEqual([
            { chr: "chr1", start: 900, end: 1100 },
            { chr: "chr2", start: 499900, end: 500100 },
        ]);
    });

    it("maps sv types to distinct glyph styles", () => {
        expect(getStructuralVariantStyle("DEL")).toEqual({ fillStyle: "#d94841", strokeStyle: "#922b21", glyph: "capped_line" });
        expect(getStructuralVariantStyle("INS")).toEqual({ fillStyle: "#8e44ad", strokeStyle: "#5b2c6f", glyph: "triangle" });
        expect(getStructuralVariantStyle("INV")).toEqual({ fillStyle: "#2ecc71", strokeStyle: "#1e8449", glyph: "inversion" });
        expect(getStructuralVariantStyle("DUP")).toEqual({ fillStyle: "#2980b9", strokeStyle: "#1f618d", glyph: "double_bar" });
        expect(getStructuralVariantStyle("BND")).toEqual({ fillStyle: "#e67e22", strokeStyle: "#af601a", glyph: "breakend" });
    });
});

describe("igv_browser chart type", () => {
    it("registers required/init for genome contracts", () => {
        const chartType: any = BaseChart.types["igv_browser"];
        expect(chartType).toBeDefined();

        const dsGenomic = { genome: { genomic_location: { columns: { chr: "chr", start: "start", end: "end" } } } } as any;
        const dsSvs = { genome: { svs: { sv_columns: { chr1: "chrA", pos1: "s", pos2: "e", chr2: "chrB" } } } } as any;
        const dsNone = { genome: {} } as any;

        expect(chartType.required(dsGenomic)).toBeTruthy();
        expect(chartType.required(dsSvs)).toBeTruthy();
        expect(chartType.required(dsNone)).toBeFalsy();

        const confA: any = {};
        chartType.init(confA, dsGenomic, undefined);
        expect(confA.param).toEqual(["chr", "start", "end"]);

        const confB: any = {};
        chartType.init(confB, dsSvs, undefined);
        expect(confB.param).toEqual(["chrA", "s", "e", "chrB"]);
    });
});
