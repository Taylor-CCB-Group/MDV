import { describe, expect, test } from "vitest";
import {
	computeEffectSize,
	detectDataIsLog1p,
	log2FoldChange,
	computeGeneStats,
} from "@/datastore/dgeStats";

// ── sgtype-driven detection logic ───────────────────────────────────────────
// Mirrors the decision tree in dgeIntegration.ts: runDGEOnDataStore

function resolveDataIsLog1p(
	sgType: string | undefined,
	probeValues: Float32Array | null,
): { dataIsLog1p: boolean; source: "metadata" | "probe" | "default" } {
	if (sgType === "sparse") {
		return { dataIsLog1p: true, source: "metadata" };
	}
	if (probeValues) {
		const detection = detectDataIsLog1p(probeValues);
		if (detection !== null) {
			return { dataIsLog1p: detection, source: "probe" };
		}
	}
	return { dataIsLog1p: true, source: "default" };
}

describe("sgtype-driven data detection", () => {
	test("sparse sgtype -> dataIsLog1p=true without probing", () => {
		const result = resolveDataIsLog1p("sparse", null);
		expect(result.dataIsLog1p).toBe(true);
		expect(result.source).toBe("metadata");
	});

	test("sparse sgtype takes precedence even if probe data has negatives", () => {
		const zScoredData = new Float32Array([-1.5, 0.2, -0.8, 1.1]);
		const result = resolveDataIsLog1p("sparse", zScoredData);
		expect(result.dataIsLog1p).toBe(true);
		expect(result.source).toBe("metadata");
	});

	test("dense sgtype -> falls back to probing", () => {
		const log1pData = new Float32Array([0, 0.5, 1.2, 0, 3.1]);
		const result = resolveDataIsLog1p("dense", log1pData);
		expect(result.dataIsLog1p).toBe(true);
		expect(result.source).toBe("probe");
	});

	test("dense sgtype with z-scored data -> detects z-scored via probe", () => {
		const zScoredData = new Float32Array([-1.5, 0.2, -0.8, 1.1]);
		const result = resolveDataIsLog1p("dense", zScoredData);
		expect(result.dataIsLog1p).toBe(false);
		expect(result.source).toBe("probe");
	});

	test("undefined sgtype -> falls back to probing", () => {
		const log1pData = new Float32Array([0, 0.5, 1.2, 0, 3.1]);
		const result = resolveDataIsLog1p(undefined, log1pData);
		expect(result.dataIsLog1p).toBe(true);
		expect(result.source).toBe("probe");
	});

	test("dense sgtype with all-NaN probe data -> defaults to log1p", () => {
		const allNaN = new Float32Array([NaN, NaN, NaN, NaN]);
		const result = resolveDataIsLog1p("dense", allNaN);
		expect(result.dataIsLog1p).toBe(true);
		expect(result.source).toBe("default");
	});

	test("undefined sgtype with no probe data -> defaults to log1p", () => {
		const result = resolveDataIsLog1p(undefined, null);
		expect(result.dataIsLog1p).toBe(true);
		expect(result.source).toBe("default");
	});
});

describe("computeEffectSize", () => {
	test("uses log2FC formula for log1p data", () => {
		const result = computeEffectSize(2.0, 1.0, true);
		const expected = log2FoldChange(2.0, 1.0);
		expect(result).toBe(expected);
	});

	test("uses mean difference for z-scored data", () => {
		const result = computeEffectSize(2.0, 1.0, false);
		expect(result).toBe(1.0);
	});

	test("mean difference works with negative values (z-scored)", () => {
		const result = computeEffectSize(-0.5, 0.3, false);
		expect(result).toBeCloseTo(-0.8, 10);
	});

	test("mean difference is zero when means are equal", () => {
		expect(computeEffectSize(3.5, 3.5, false)).toBe(0);
	});

	test("log2FC is zero when means are equal", () => {
		expect(computeEffectSize(1.5, 1.5, true)).toBeCloseTo(0, 6);
	});

	test("log2FC is positive when target > reference (log1p data)", () => {
		expect(computeEffectSize(3.0, 1.0, true)).toBeGreaterThan(0);
	});

	test("mean diff is positive when target > reference (z-scored data)", () => {
		expect(computeEffectSize(3.0, 1.0, false)).toBeGreaterThan(0);
	});
});

describe("detectDataIsLog1p", () => {
	test("returns true for all non-negative values", () => {
		const data = new Float32Array([0, 0.5, 1.0, 2.3, 0, 4.5]);
		expect(detectDataIsLog1p(data)).toBe(true);
	});

	test("returns false when any value is negative", () => {
		const data = new Float32Array([0, 0.5, -0.1, 2.3, 0, 4.5]);
		expect(detectDataIsLog1p(data)).toBe(false);
	});

	test("returns true for all zeros", () => {
		const data = new Float32Array([0, 0, 0, 0]);
		expect(detectDataIsLog1p(data)).toBe(true);
	});

	test("returns null for empty array (inconclusive)", () => {
		const data = new Float32Array([]);
		expect(detectDataIsLog1p(data)).toBe(null);
	});

	test("returns null for all-NaN array (inconclusive)", () => {
		const data = new Float32Array([NaN, NaN, NaN, NaN]);
		expect(detectDataIsLog1p(data)).toBe(null);
	});

	test("returns true for sparse data with NaN and non-negative values", () => {
		const data = new Float32Array([NaN, 0.5, NaN, NaN, 2.3, NaN]);
		expect(detectDataIsLog1p(data)).toBe(true);
	});

	test("returns false for sparse data with NaN and negative values", () => {
		const data = new Float32Array([NaN, -0.5, NaN, NaN, 2.3, NaN]);
		expect(detectDataIsLog1p(data)).toBe(false);
	});

	test("returns false for typical z-scored data with negatives", () => {
		const data = new Float32Array([-1.5, 0.2, -0.8, 1.1, -0.3, 0.9]);
		expect(detectDataIsLog1p(data)).toBe(false);
	});
});

describe("computeGeneStats with dataIsLog1p flag", () => {
	function makeArrays(values: number[], groups: number[]) {
		return {
			v: new Float32Array(values),
			g: new Uint8Array(groups),
			f: new Uint8Array(values.length),
		};
	}

	test("uses log2FC by default (dataIsLog1p=true)", () => {
		const values = [2, 3, 4, 5, 6, 10, 11, 12, 13, 14];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("Gene1", v, g, f, 0, -1, true);
		const expectedEffectSize = log2FoldChange(4, 12);
		expect(result.effectSize).toBeCloseTo(expectedEffectSize, 3);
	});

	test("uses mean difference when dataIsLog1p=false", () => {
		const values = [-1, -0.5, 0, 0.5, 1, 2, 2.5, 3, 3.5, 4];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("Gene1", v, g, f, 0, -1, false);
		expect(result.effectSize).toBeCloseTo(0 - 3, 3); // mean(group0) - mean(group1)
	});

	test("effect size sign is consistent regardless of flag", () => {
		const values = [1, 1, 1, 1, 1, 5, 5, 5, 5, 5];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const resultLog1p = computeGeneStats("Gene1", v, g, f, 0, -1, true);
		const resultZscored = computeGeneStats("Gene1", v, g, f, 0, -1, false);

		expect(Math.sign(resultLog1p.effectSize)).toBe(Math.sign(resultZscored.effectSize));
		expect(resultLog1p.effectSize).toBeLessThan(0); // target (1) < reference (5)
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

		// GeneA at index 0
		expect(effectSizeColumn[0]).toBeCloseTo(-1.0, 5);
		expect(pvalColumn[0]).toBeCloseTo(0.01, 5);

		// GeneB at index 1 should be NaN (not in results)
		expect(effectSizeColumn[1]).toBeNaN();
		expect(pvalColumn[1]).toBeNaN();

		// GeneC at index 2
		expect(effectSizeColumn[2]).toBeCloseTo(0.5, 5);

		// GeneD at index 3
		expect(effectSizeColumn[3]).toBeCloseTo(2.5, 5);

		// GeneE at index 4 should be NaN (not in results)
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
