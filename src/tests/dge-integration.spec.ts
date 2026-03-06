import { describe, expect, test } from "vitest";
import {
	computeEffectSize,
	log2FoldChange,
	computeGeneStats,
} from "@/datastore/dgeStats";

// ── Data type detection logic ───────────────────────────────────────────────
// Mirrors the multi-gene probe loop in dgeIntegration.ts: runDGEOnDataStore.
// The real code scans values across multiple genes, tracking globalMaxVal and
// globalHasNegative, and only makes a decision after examining all probes (or
// early-exits on a definitive signal).

function resolveDataTypeMultiGene(
	geneArrays: (Float32Array | null)[],
): "log1p" | "linear" | "zscored" {
	let globalMaxVal = -Infinity;
	let globalHasNegative = false;
	let anyFinite = false;

	for (const arr of geneArrays) {
		if (!arr) continue;
		for (let j = 0; j < arr.length; j++) {
			const v = arr[j];
			if (Number.isNaN(v)) continue;
			anyFinite = true;
			if (v < 0) { globalHasNegative = true; break; }
			if (v > globalMaxVal) globalMaxVal = v;
		}
		if (globalHasNegative || globalMaxVal > 20) break;
	}

	if (!anyFinite) return "log1p";
	if (globalHasNegative) return "zscored";
	return globalMaxVal > 20 ? "linear" : "log1p";
}

describe("multi-gene data type detection", () => {
	test("detects linear when LATER genes have high values (raw count data)", () => {
		const gene0 = new Float32Array([NaN, NaN, NaN, NaN]);
		const gene1 = new Float32Array([NaN, 1, NaN, 2]);
		const gene2 = new Float32Array([NaN, NaN, 3, NaN]);
		const gene3 = new Float32Array([NaN, 50, NaN, 116]);
		expect(resolveDataTypeMultiGene([gene0, gene1, gene2, gene3])).toBe("linear");
	});

	test("first gene with low values does NOT prematurely decide log1p", () => {
		const gene0 = new Float32Array([NaN, 2.0, NaN, NaN]);
		const gene1 = new Float32Array([NaN, NaN, NaN, NaN]);
		const gene2 = new Float32Array([NaN, 100, NaN, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1, gene2])).toBe("linear");
	});

	test("detects log1p when ALL genes have small values", () => {
		const gene0 = new Float32Array([NaN, 0.5, NaN, 1.2]);
		const gene1 = new Float32Array([NaN, NaN, 3.1, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1])).toBe("log1p");
	});

	test("detects zscored immediately on negative value", () => {
		const gene0 = new Float32Array([NaN, -0.5, NaN, 2.3]);
		const gene1 = new Float32Array([NaN, 50, NaN, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1])).toBe("zscored");
	});

	test("all-NaN across all genes defaults to log1p", () => {
		const gene0 = new Float32Array([NaN, NaN, NaN]);
		const gene1 = new Float32Array([NaN, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1])).toBe("log1p");
	});

	test("null buffers are skipped gracefully", () => {
		const gene0 = null;
		const gene1 = new Float32Array([NaN, 25, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1])).toBe("linear");
	});
});

describe("computeEffectSize", () => {
	test("uses log2FC formula for log1p data", () => {
		const result = computeEffectSize(2.0, 1.0, "log1p");
		const expected = log2FoldChange(2.0, 1.0, true);
		expect(result).toBe(expected);
	});

	test("uses log2FC formula for linear data (without expm1)", () => {
		const result = computeEffectSize(2.0, 1.0, "linear");
		const expected = log2FoldChange(2.0, 1.0, false);
		expect(result).toBe(expected);
	});

	test("uses mean difference for z-scored data", () => {
		const result = computeEffectSize(2.0, 1.0, "zscored");
		expect(result).toBe(1.0);
	});

	test("mean difference works with negative values (z-scored)", () => {
		const result = computeEffectSize(-0.5, 0.3, "zscored");
		expect(result).toBeCloseTo(-0.8, 10);
	});

	test("mean difference is zero when means are equal", () => {
		expect(computeEffectSize(3.5, 3.5, "zscored")).toBe(0);
	});

	test("log2FC is zero when means are equal", () => {
		expect(computeEffectSize(1.5, 1.5, "log1p")).toBeCloseTo(0, 6);
	});

	test("log2FC is positive when target > reference (log1p data)", () => {
		expect(computeEffectSize(3.0, 1.0, "log1p")).toBeGreaterThan(0);
	});

	test("mean diff is positive when target > reference (z-scored data)", () => {
		expect(computeEffectSize(3.0, 1.0, "zscored")).toBeGreaterThan(0);
	});
});

describe("computeGeneStats with dataType parameter", () => {
	function makeArrays(values: number[], groups: number[]) {
		return {
			v: new Float32Array(values),
			g: new Uint8Array(groups),
			f: new Uint8Array(values.length),
		};
	}

	test("uses log2FC for log1p data type", () => {
		const values = [2, 3, 4, 5, 6, 10, 11, 12, 13, 14];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("Gene1", v, g, f, 0, -1, "log1p");
		const expectedEffectSize = log2FoldChange(4, 12, true);
		expect(result.effectSize).toBeCloseTo(expectedEffectSize, 3);
	});

	test("uses log2FC for linear data type (without expm1)", () => {
		const values = [2, 3, 4, 5, 6, 10, 11, 12, 13, 14];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("Gene1", v, g, f, 0, -1, "linear");
		const expectedEffectSize = log2FoldChange(4, 12, false);
		expect(result.effectSize).toBeCloseTo(expectedEffectSize, 3);
	});

	test("uses mean difference for zscored data type", () => {
		const values = [-1, -0.5, 0, 0.5, 1, 2, 2.5, 3, 3.5, 4];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("Gene1", v, g, f, 0, -1, "zscored");
		expect(result.effectSize).toBeCloseTo(0 - 3, 3);
	});

	test("effect size sign is consistent across data types", () => {
		const values = [1, 1, 1, 1, 1, 5, 5, 5, 5, 5];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const resultLog1p = computeGeneStats("Gene1", v, g, f, 0, -1, "log1p");
		const resultLinear = computeGeneStats("Gene1", v, g, f, 0, -1, "linear");
		const resultZscored = computeGeneStats("Gene1", v, g, f, 0, -1, "zscored");

		expect(resultLog1p.effectSize).toBeLessThan(0);
		expect(resultLinear.effectSize).toBeLessThan(0);
		expect(resultZscored.effectSize).toBeLessThan(0);
	});
});

describe("gene ordering alignment", () => {
	test("maps DGE results to correct gene row indices", () => {
		const valueToRowIndex = new Map<string, number>([
			["GeneA", 0],
			["GeneB", 1],
			["GeneC", 2],
			["GeneD", 3],
			["GeneE", 4],
		]);
		const n = 5;

		const dgeResults = [
			{ gene: "GeneD", effectSize: 2.5, pval: 0.001 },
			{ gene: "GeneA", effectSize: -1.0, pval: 0.01 },
			{ gene: "GeneC", effectSize: 0.5, pval: 0.05 },
		];

		const effectSizeColumn = new Float32Array(n).fill(NaN);
		const pvalColumn = new Float32Array(n).fill(NaN);

		for (const result of dgeResults) {
			const idx = valueToRowIndex.get(result.gene);
			if (idx !== undefined) {
				effectSizeColumn[idx] = result.effectSize;
				pvalColumn[idx] = result.pval;
			}
		}

		expect(effectSizeColumn[0]).toBeCloseTo(-1.0, 5);
		expect(pvalColumn[0]).toBeCloseTo(0.01, 5);
		expect(effectSizeColumn[1]).toBeNaN();
		expect(pvalColumn[1]).toBeNaN();
		expect(effectSizeColumn[2]).toBeCloseTo(0.5, 5);
		expect(effectSizeColumn[3]).toBeCloseTo(2.5, 5);
		expect(effectSizeColumn[4]).toBeNaN();
	});

	test("handles duplicate gene names by using last occurrence", () => {
		const valueToRowIndex = new Map<string, number>([
			["Gene1", 0],
			["Gene2", 1],
		]);

		const dgeResults = [
			{ gene: "Gene1", effectSize: 1.0 },
			{ gene: "Gene1", effectSize: 2.0 },
		];

		const column = new Float32Array(2).fill(NaN);
		for (const result of dgeResults) {
			const idx = valueToRowIndex.get(result.gene);
			if (idx !== undefined) {
				column[idx] = result.effectSize;
			}
		}

		expect(column[0]).toBeCloseTo(2.0, 5);
	});

	test("handles genes not in the DataStore gracefully", () => {
		const valueToRowIndex = new Map<string, number>([
			["Gene1", 0],
		]);

		const column = new Float32Array(1).fill(NaN);
		const idx = valueToRowIndex.get("UnknownGene");
		expect(idx).toBeUndefined();
		expect(column[0]).toBeNaN();
	});
});
