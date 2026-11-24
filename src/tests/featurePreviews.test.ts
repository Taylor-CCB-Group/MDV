import { afterEach, describe, expect, test } from "vitest";
import BaseChart from "@/charts/BaseChart";
import "@/spatial_stats/CrossPcfChart";
import {
    convertLegacyChartsSetting,
    featurePreviewsSetting,
    isSettingListed,
    showFeaturePreviews,
} from "@/lib/browserSettings";

describe("feature previews", () => {
    afterEach(() => {
        featurePreviewsSetting.set("hide");
        convertLegacyChartsSetting.set("off");
    });

    test("are hidden by default", () => {
        expect(showFeaturePreviews()).toBe(false);
    });

    test("the cross-PCF chart is a preview, named without a qualifier", () => {
        const type = BaseChart.types.cross_pcf;
        expect(type.preview).toBe(true);
        expect(type.name).toBe("Cross-PCF matrix");
    });

    test("a preview setting is listed only while previews show, unless it has been changed", () => {
        expect(isSettingListed(featurePreviewsSetting)).toBe(true);
        expect(isSettingListed(convertLegacyChartsSetting)).toBe(false);
        featurePreviewsSetting.set("show");
        expect(isSettingListed(convertLegacyChartsSetting)).toBe(true);
        convertLegacyChartsSetting.set("on");
        featurePreviewsSetting.set("hide");
        // still listed, so it can be turned back off
        expect(isSettingListed(convertLegacyChartsSetting)).toBe(true);
    });
});
