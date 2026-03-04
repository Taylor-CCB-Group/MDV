import { describe, expect, test } from "vitest";
import {
	welfordCreate,
	welfordAccumulate,
	welfordFinalize,
	welchTTest,
	log2FoldChange,
	clampEffectSize,
	benjaminiHochberg,
	tDistPValue,
	computeGeneStats,
	type GroupStats,
} from "@/datastore/dgeStats";

// ── Welford's Online Algorithm ──────────────────────────────────────────────

describe("Welford online mean/variance", () => {
	test("computes correct mean and variance for a known sequence", () => {
		const values = [2.3, 3.1, 2.8, 3.5, 2.9, 3.2, 2.7, 3.0, 2.6, 3.3];
		const acc = welfordCreate();
		for (const v of values) welfordAccumulate(acc, v);
		const stats = welfordFinalize(acc);

		expect(stats.n).toBe(10);
		expect(stats.mean).toBeCloseTo(2.94, 10);
		expect(stats.variance).toBeCloseTo(0.12711111111111112, 10);
	});

	test("handles a single value (variance = 0)", () => {
		const acc = welfordCreate();
		welfordAccumulate(acc, 42.0);
		const stats = welfordFinalize(acc);

		expect(stats.n).toBe(1);
		expect(stats.mean).toBe(42.0);
		expect(stats.variance).toBe(0);
	});

	test("handles all identical values", () => {
		const acc = welfordCreate();
		for (let i = 0; i < 100; i++) welfordAccumulate(acc, 7.0);
		const stats = welfordFinalize(acc);

		expect(stats.n).toBe(100);
		expect(stats.mean).toBe(7.0);
		expect(stats.variance).toBeCloseTo(0, 12);
	});

	test("handles empty accumulator", () => {
		const stats = welfordFinalize(welfordCreate());
		expect(stats.n).toBe(0);
		expect(stats.mean).toBe(0);
		expect(stats.variance).toBe(0);
	});

	test("handles negative values", () => {
		const acc = welfordCreate();
		for (const v of [-3, -1, 0, 1, 3]) welfordAccumulate(acc, v);
		const stats = welfordFinalize(acc);

		expect(stats.mean).toBeCloseTo(0, 10);
		expect(stats.variance).toBeCloseTo(5, 10);
	});

	test("handles large values without catastrophic cancellation", () => {
		const acc = welfordCreate();
		const base = 1e8;
		const values = [base + 1, base + 2, base + 3, base + 4, base + 5];
		for (const v of values) welfordAccumulate(acc, v);
		const stats = welfordFinalize(acc);

		expect(stats.mean).toBeCloseTo(base + 3, 5);
		expect(stats.variance).toBeCloseTo(2.5, 5);
	});
});

// ── t-Distribution CDF ─────────────────────────────────────────────────────

describe("tDistPValue", () => {
	test("known reference values from scipy.stats.t", () => {
		const cases: [number, number, number][] = [
			// [t, df, expected_two_sided_pvalue]
			[2.0, 10, 0.07338803477074039],
			[3.0, 5, 0.03009924789746257],
			[1.96, 100, 0.052778901366229654],
			[0.5, 30, 0.6207230048851273],
			[10.0, 3, 0.0021283990584141494],
		];

		for (const [t, df, expected] of cases) {
			const result = tDistPValue(t, df);
			expect(result).toBeCloseTo(expected, 5);
		}
	});

	test("symmetric for positive and negative t", () => {
		expect(tDistPValue(2.5, 15)).toBeCloseTo(tDistPValue(-2.5, 15), 10);
	});

	test("t=0 gives p=1", () => {
		expect(tDistPValue(0, 10)).toBeCloseTo(1, 10);
	});

	test("very large t gives p near 0", () => {
		const p = tDistPValue(50, 100);
		expect(p).toBeLessThan(1e-10);
	});

	test("returns NaN for invalid inputs", () => {
		expect(tDistPValue(NaN, 10)).toBeNaN();
		expect(tDistPValue(2, NaN)).toBeNaN();
		expect(tDistPValue(2, -1)).toBeNaN();
	});
});

// ── Welch's t-test ──────────────────────────────────────────────────────────

describe("welchTTest", () => {
	test("matches scipy reference for equal-size groups", () => {
		const g1: GroupStats = { mean: 2.94, variance: 0.12711111111111112, n: 10 };
		const g2: GroupStats = { mean: 1.45, variance: 0.09166666666666665, n: 10 };
		const result = welchTTest(g1, g2);

		expect(result.tStatistic).toBeCloseTo(10.073599143070807, 4);
		expect(result.pValue).toBeCloseTo(1.0273445355608373e-8, 3);
	});

	test("matches scipy reference for unequal-size groups", () => {
		const g1: GroupStats = { mean: 7.0, variance: 2.5, n: 5 };
		const g2: GroupStats = { mean: 4.5, variance: 6.0, n: 8 };
		const result = welchTTest(g1, g2);

		expect(result.tStatistic).toBeCloseTo(2.23606797749979, 4);
		expect(result.pValue).toBeCloseTo(0.047155428922714536, 4);
	});

	test("returns NaN for groups with < 2 observations", () => {
		const g1: GroupStats = { mean: 1.0, variance: 0, n: 1 };
		const g2: GroupStats = { mean: 2.0, variance: 0.5, n: 10 };
		const result = welchTTest(g1, g2);

		expect(result.pValue).toBeNaN();
	});

	test("returns p=1 for identical groups (zero variance)", () => {
		const g1: GroupStats = { mean: 5.0, variance: 0, n: 100 };
		const g2: GroupStats = { mean: 5.0, variance: 0, n: 100 };
		const result = welchTTest(g1, g2);

		expect(result.tStatistic).toBe(0);
		expect(result.pValue).toBe(1);
	});

	test("nearly identical groups give high p-value", () => {
		const g1: GroupStats = { mean: 1.003, variance: 2.5e-6, n: 5 };
		const g2: GroupStats = { mean: 1.002, variance: 2.5e-6, n: 5 };
		const result = welchTTest(g1, g2);

		expect(result.pValue).toBeCloseTo(0.3465935070873849, 1);
	});
});

// ── Log2 Fold Change ────────────────────────────────────────────────────────

describe("log2FoldChange", () => {
	test("matches Scanpy formula: log2((expm1(mean_g) + 1e-9) / (expm1(mean_r) + 1e-9))", () => {
		const cases: [number, number, number][] = [
			[2.0, 1.0, 1.894636123358204],
			[0.0, 0.0, 0.0],
			[3.0, 0.5, 4.8787372176909996],
			[0.1, 0.1, 0.0],
		];

		for (const [meanTarget, meanRef, expected] of cases) {
			expect(log2FoldChange(meanTarget, meanRef)).toBeCloseTo(expected, 6);
		}
	});

	test("positive when target > reference", () => {
		expect(log2FoldChange(3.0, 1.0)).toBeGreaterThan(0);
	});

	test("negative when target < reference", () => {
		expect(log2FoldChange(1.0, 3.0)).toBeLessThan(0);
	});

	test("handles zero means (both zero gives 0)", () => {
		expect(log2FoldChange(0, 0)).toBeCloseTo(0, 10);
	});
});

// ── clampEffectSize ─────────────────────────────────────────────────────────

describe("clampEffectSize", () => {
	test("clamps extreme positive values to 10", () => {
		expect(clampEffectSize(23.5)).toBe(10);
	});

	test("clamps extreme negative values to -10", () => {
		expect(clampEffectSize(-25.0)).toBe(-10);
	});

	test("passes through moderate values unchanged", () => {
		expect(clampEffectSize(3.5)).toBe(3.5);
		expect(clampEffectSize(-2.1)).toBe(-2.1);
		expect(clampEffectSize(0)).toBe(0);
	});

	test("boundary values", () => {
		expect(clampEffectSize(10)).toBe(10);
		expect(clampEffectSize(-10)).toBe(-10);
	});
});

// ── Benjamini-Hochberg Correction ───────────────────────────────────────────

describe("benjaminiHochberg", () => {
	test("matches statsmodels reference for ordered p-values", () => {
		const pvals = [0.001, 0.01, 0.03, 0.05, 0.1, 0.5, 0.8];
		const expected = [0.007, 0.035, 0.07, 0.0875, 0.14, 0.583333333, 0.8];
		const result = benjaminiHochberg(pvals);

		for (let i = 0; i < pvals.length; i++) {
			expect(result[i]).toBeCloseTo(expected[i], 6);
		}
	});

	test("preserves original index ordering", () => {
		const pvals = [0.5, 0.001, 0.1, 0.01];
		const result = benjaminiHochberg(pvals);

		// Smallest p-value at index 1 should have smallest adjusted p-value
		expect(result[1]).toBeLessThan(result[3]);
		expect(result[3]).toBeLessThan(result[2]);
		expect(result[2]).toBeLessThan(result[0]);
	});

	test("all p-values of 1 remain 1", () => {
		const result = benjaminiHochberg([1, 1, 1]);
		expect(result).toEqual([1, 1, 1]);
	});

	test("single p-value is unchanged", () => {
		expect(benjaminiHochberg([0.05])).toEqual([0.05]);
	});

	test("empty array returns empty", () => {
		expect(benjaminiHochberg([])).toEqual([]);
	});

	test("adjusted p-values never exceed 1", () => {
		const pvals = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
		const result = benjaminiHochberg(pvals);
		for (const p of result) {
			expect(p).toBeLessThanOrEqual(1);
		}
	});

	test("adjusted p-values are monotonically non-decreasing when input is sorted", () => {
		const pvals = [0.001, 0.005, 0.01, 0.02, 0.05];
		const result = benjaminiHochberg(pvals);
		for (let i = 1; i < result.length; i++) {
			expect(result[i]).toBeGreaterThanOrEqual(result[i - 1]);
		}
	});
});

// ── computeGeneStats (end-to-end single gene) ──────────────────────────────

describe("computeGeneStats", () => {
	function makeArrays(
		values: number[],
		groups: number[],
		filter?: number[],
	) {
		const v = new Float32Array(values);
		const g = new Uint8Array(groups);
		const f = filter ? new Uint8Array(filter) : new Uint8Array(values.length);
		return { v, g, f };
	}

	test("computes correct stats for two clear groups", () => {
		// Group 0: [2, 3, 4, 5, 6] mean=4, Group 1: [10, 11, 12, 13, 14] mean=12
		const values = [2, 3, 4, 5, 6, 10, 11, 12, 13, 14];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("TestGene", v, g, f, 0, -1);

		expect(result.gene).toBe("TestGene");
		expect(result.meanTarget).toBeCloseTo(4, 5);
		expect(result.meanReference).toBeCloseTo(12, 5);
		expect(result.effectSize).toBeLessThan(0); // target < reference
		expect(result.pval).toBeLessThan(0.001);
	});

	test("respects filterArray (skips filtered cells)", () => {
		const values = [100, 2, 3, 4, 5, 10, 11, 12, 13, 14];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		// Filter out the outlier at index 0
		const filter = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
		const { v, g, f } = makeArrays(values, groups, filter);

		const result = computeGeneStats("TestGene", v, g, f, 0, -1);

		// Without filter: mean would be much higher
		// With filter: group 0 is [2,3,4,5] mean=3.5
		expect(result.meanTarget).toBeCloseTo(3.5, 4);
	});

	test("handles specific reference group (not rest)", () => {
		const values = [1, 2, 3, 10, 11, 12, 50, 51, 52];
		const groups = [0, 0, 0, 1, 1, 1, 2, 2, 2];
		const { v, g, f } = makeArrays(values, groups);

		// Target=0, Reference=1 (ignore group 2)
		const result = computeGeneStats("TestGene", v, g, f, 0, 1);

		expect(result.meanTarget).toBeCloseTo(2, 4);
		expect(result.meanReference).toBeCloseTo(11, 4);
	});

	test("treats NaN values as zero (sparse data convention)", () => {
		const values = [1, NaN, 3, 10, NaN, 12];
		const groups = [0, 0, 0, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("TestGene", v, g, f, 0, -1);

		// NaN treated as 0: group 0 = [1, 0, 3] mean=4/3, group 1 = [10, 0, 12] mean=22/3
		expect(result.meanTarget).toBeCloseTo(4 / 3, 4);
		expect(result.meanReference).toBeCloseTo(22 / 3, 4);
	});
});
