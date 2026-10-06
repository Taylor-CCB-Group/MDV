import { describe, expect, it } from "vitest";
import {
    layoutXLabels,
    type MeasureText,
    stripSharedSubgroupSuffix,
    truncateMiddle,
    X_MARGIN_MAX_PX,
    X_MIN_FONT_SIZE,
    type XLabelLayoutInput,
} from "@/charts/dotplot/xAxisLabels";

// monospace stand-in for canvas measurement: every character is 0.6em wide
const measure: MeasureText = (text, fontSize) => text.length * fontSize * 0.6;

function subgroupColumn(gene: string, subgroup: string, index = 0) {
    return {
        field: `${subgroup}|${gene}(${subgroup})|${index}`,
        name: `${gene}(${subgroup})`,
        inSubgroup: true,
    };
}

function plainColumn(name: string) {
    return { field: name, name, inSubgroup: false };
}

describe("stripSharedSubgroupSuffix", () => {
    it("drops the suffix when every column shares one subgroup", () => {
        const columns = ["CD4", "CD8A", "GZMB"].map((g, i) =>
            subgroupColumn(g, "rna_expr", i),
        );
        expect(stripSharedSubgroupSuffix(columns)).toEqual(["CD4", "CD8A", "GZMB"]);
    });

    it("leaves plain columns unchanged", () => {
        const columns = [plainColumn("age(years)"), plainColumn("score")];
        expect(stripSharedSubgroupSuffix(columns)).toEqual(["age(years)", "score"]);
    });

    it("keeps suffixes when columns come from different subgroups", () => {
        const columns = [subgroupColumn("CD4", "rna_expr"), subgroupColumn("CD4", "gs")];
        expect(stripSharedSubgroupSuffix(columns)).toEqual(["CD4(rna_expr)", "CD4(gs)"]);
    });

    it("keeps suffixes when subgroup and plain columns are mixed", () => {
        const columns = [subgroupColumn("CD4", "rna_expr"), plainColumn("score")];
        expect(stripSharedSubgroupSuffix(columns)).toEqual(["CD4(rna_expr)", "score"]);
    });

    it("keeps a renamed column that no longer ends in the suffix", () => {
        const columns = [
            subgroupColumn("CD4", "rna_expr"),
            { ...subgroupColumn("CD8A", "rna_expr", 1), name: "my favourite gene" },
        ];
        expect(stripSharedSubgroupSuffix(columns)).toEqual(["CD4", "my favourite gene"]);
    });

    it("does not reduce a name to nothing", () => {
        const columns = [{ ...subgroupColumn("CD4", "rna_expr"), name: "(rna_expr)" }];
        expect(stripSharedSubgroupSuffix(columns)).toEqual(["(rna_expr)"]);
    });
});

describe("truncateMiddle", () => {
    it("returns text that fits unchanged", () => {
        expect(truncateMiddle("CD4", 100, 10, measure)).toBe("CD4");
    });

    it("cuts from the middle, keeping as many characters as fit", () => {
        // 60px at 6px per char = 10 chars: 5 head + "…" + 4 tail
        expect(truncateMiddle("ENSG00000228253-AS1", 60, 10, measure)).toBe("ENSG0…-AS1");
    });

    it("keeps labels with a shared prefix distinguishable", () => {
        const a = truncateMiddle("ENSG00000228253-AS1", 60, 10, measure);
        const b = truncateMiddle("ENSG00000254876.2", 60, 10, measure);
        expect(a).not.toBe(b);
    });

    it("falls back to an ellipsis when nothing else fits", () => {
        expect(truncateMiddle("ENSG00000228253-AS1", 3, 10, measure)).toBe("…");
    });
});

describe("layoutXLabels", () => {
    const base: XLabelLayoutInput = {
        labels: [],
        plotLeft: 110,
        plotWidth: 560,
        chartHeight: 400,
        fontSize: 10,
        titleSize: 0,
        truncate: true,
        thin: false,
    };
    const longLabels = Array.from({ length: 10 }, (_, i) => `LONG_GENE_NAME_${i}_abcd`);

    it("lays labels flat at the user's font size when they fit", () => {
        const layout = layoutXLabels(
            { ...base, labels: ["CD4", "CD8A", "GZMB", "NKG7"] },
            measure,
        );
        expect(layout).toEqual({
            angle: 0,
            fontSize: 10,
            margin: 23, // tick offset 9 + font 10 + gap 4
            labels: ["CD4", "CD8A", "GZMB", "NKG7"],
        });
    });

    it("adds room for the axis title", () => {
        const layout = layoutXLabels({ ...base, labels: ["CD4"], titleSize: 19 }, measure);
        expect(layout.margin).toBe(42);
    });

    it("tilts labels that do not fit flat", () => {
        const layout = layoutXLabels({ ...base, labels: longLabels }, measure);
        expect(layout.angle).toBe(45);
        expect(layout.fontSize).toBe(10);
    });

    it("caps the margin and shortens labels when truncating", () => {
        const layout = layoutXLabels({ ...base, labels: longLabels }, measure);
        expect(layout.margin).toBeLessThanOrEqual(X_MARGIN_MAX_PX);
        for (const label of layout.labels) {
            expect(label).toContain("…");
        }
    });

    it("caps the margin at 30% of a short chart's height", () => {
        const layout = layoutXLabels(
            { ...base, labels: longLabels, chartHeight: 200 },
            measure,
        );
        expect(layout.margin).toBeLessThanOrEqual(60);
    });

    it("shortens the first labels more when they would reach past the left edge", () => {
        const layout = layoutXLabels(
            { ...base, labels: longLabels, plotLeft: 20 },
            measure,
        );
        expect(layout.labels[0].length).toBeLessThan(layout.labels[1].length);
    });

    it("shows full labels and grows the margin when not truncating", () => {
        const layout = layoutXLabels(
            { ...base, labels: longLabels, truncate: false },
            measure,
        );
        expect(layout.labels).toEqual(longLabels);
        expect(layout.margin).toBeGreaterThan(X_MARGIN_MAX_PX);
    });

    it("leaves room for the dots when full labels would not fit", () => {
        const layout = layoutXLabels(
            { ...base, labels: longLabels, truncate: false, chartHeight: 100 },
            measure,
        );
        expect(layout.margin).toBe(60);
    });

    it("shrinks the font to fit narrow columns", () => {
        // 20 columns of 13px: 13·sin45 / 1.1 ≈ 8.4px
        const layout = layoutXLabels(
            { ...base, labels: Array(20).fill("LONGERNAME"), plotWidth: 260 },
            measure,
        );
        expect(layout.fontSize).toBeCloseTo((13 * Math.SQRT1_2) / 1.1);
        expect(layout.fontSize).toBeLessThan(10);
    });

    it("never enlarges labels past a user size below the minimum", () => {
        const layout = layoutXLabels(
            { ...base, labels: Array(60).fill("GENE00"), plotWidth: 315, fontSize: 5 },
            measure,
        );
        expect(layout.fontSize).toBe(5);
    });

    describe("crowded columns", () => {
        // 60 columns in 315px: 5.25px each, too narrow even at the minimum font
        const crowded = { ...base, plotWidth: 315, labels: Array(60).fill("GENE00") };

        it("labels every column by default, at the minimum font", () => {
            const layout = layoutXLabels(crowded, measure);
            expect(layout.fontSize).toBe(X_MIN_FONT_SIZE);
            expect(layout.labels.every((l) => l !== "")).toBe(true);
        });

        it("labels only every nth column when thinning", () => {
            const layout = layoutXLabels({ ...crowded, thin: true }, measure);
            // 7px · 1.1 / (5.25 · sin45) → every 3rd column
            const labelled = layout.labels
                .map((l, i) => (l === "" ? null : i))
                .filter((i) => i !== null);
            expect(labelled).toHaveLength(20);
            expect(labelled.slice(0, 3)).toEqual([0, 3, 6]);
        });
    });
});
