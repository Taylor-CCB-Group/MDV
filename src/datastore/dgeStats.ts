/**
 * mdv-dge: Pure statistical functions for browser-based differential gene expression.
 *
 * All functions are pure (no side effects, no DOM, no SharedArrayBuffer dependency)
 * so they can be used in Web Workers, tests, or the main thread.
 */

// ── Welford's Online Algorithm ──────────────────────────────────────────────

export interface WelfordAccumulator {
	n: number;
	mean: number;
	m2: number;
}

export interface GroupStats {
	mean: number;
	variance: number;
	n: number;
}

export function welfordCreate(): WelfordAccumulator {
	return { n: 0, mean: 0, m2: 0 };
}

export function welfordAccumulate(acc: WelfordAccumulator, value: number): void {
	acc.n++;
	const delta = value - acc.mean;
	acc.mean += delta / acc.n;
	const delta2 = value - acc.mean;
	acc.m2 += delta * delta2;
}

export function welfordFinalize(acc: WelfordAccumulator): GroupStats {
	if (acc.n < 2) {
		return { mean: acc.mean, variance: 0, n: acc.n };
	}
	return { mean: acc.mean, variance: acc.m2 / (acc.n - 1), n: acc.n };
}

// ── t-Distribution CDF ─────────────────────────────────────────────────────

/**
 * Regularized incomplete beta function I_x(a, b) via continued fraction
 * (Lentz's algorithm). Used for the t-distribution CDF.
 */
function betaIncomplete(x: number, a: number, b: number): number {
	if (x <= 0) return 0;
	if (x >= 1) return 1;

	// Use the symmetry relation when x > (a+1)/(a+b+2) for convergence
	if (x > (a + 1) / (a + b + 2)) {
		return 1 - betaIncomplete(1 - x, b, a);
	}

	const lnBeta = lnGamma(a) + lnGamma(b) - lnGamma(a + b);
	const front = Math.exp(Math.log(x) * a + Math.log(1 - x) * b - lnBeta) / a;

	// Lentz's continued fraction
	const TINY = 1e-30;
	const EPS = 1e-14;
	let f = TINY;
	let c = 1;
	let d = 1 - (a + b) * x / (a + 1);
	if (Math.abs(d) < TINY) d = TINY;
	d = 1 / d;
	f = d;

	for (let m = 1; m <= 200; m++) {
		// Even step: a_2m
		let numerator = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m));
		d = 1 + numerator * d;
		if (Math.abs(d) < TINY) d = TINY;
		c = 1 + numerator / c;
		if (Math.abs(c) < TINY) c = TINY;
		d = 1 / d;
		f *= c * d;

		// Odd step: a_2m+1
		numerator = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1));
		d = 1 + numerator * d;
		if (Math.abs(d) < TINY) d = TINY;
		c = 1 + numerator / c;
		if (Math.abs(c) < TINY) c = TINY;
		d = 1 / d;
		const delta = c * d;
		f *= delta;

		if (Math.abs(delta - 1) < EPS) break;
	}

	return front * f;
}

/**
 * Log-gamma function using Lanczos approximation (g=7, n=9).
 */
function lnGamma(z: number): number {
	if (z < 0.5) {
		return Math.log(Math.PI / Math.sin(Math.PI * z)) - lnGamma(1 - z);
	}
	z -= 1;
	const c = [
		0.99999999999980993,
		676.5203681218851,
		-1259.1392167224028,
		771.32342877765313,
		-176.61502916214059,
		12.507343278686905,
		-0.13857109526572012,
		9.9843695780195716e-6,
		1.5056327351493116e-7,
	];
	let x = c[0];
	for (let i = 1; i < 9; i++) {
		x += c[i] / (z + i);
	}
	const t = z + 7.5;
	return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}

/**
 * Two-sided p-value from Student's t-distribution.
 */
export function tDistPValue(t: number, df: number): number {
	if (df <= 0 || !Number.isFinite(t) || !Number.isFinite(df)) return NaN;
	if (df === Number.POSITIVE_INFINITY) {
		// Normal distribution limit
		return 2 * (1 - normalCDF(Math.abs(t)));
	}
	const x = df / (df + t * t);
	const ibeta = betaIncomplete(x, df / 2, 0.5);
	return ibeta;
}

/** Standard normal CDF (for the df→∞ case). */
function normalCDF(x: number): number {
	const a1 = 0.254829592;
	const a2 = -0.284496736;
	const a3 = 1.421413741;
	const a4 = -1.453152027;
	const a5 = 1.061405429;
	const p = 0.3275911;
	const sign = x < 0 ? -1 : 1;
	x = Math.abs(x) / Math.SQRT2;
	const t = 1 / (1 + p * x);
	const y = 1 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
	return 0.5 * (1 + sign * y);
}

// ── Welch's t-Test ──────────────────────────────────────────────────────────

export interface TTestResult {
	tStatistic: number;
	degreesOfFreedom: number;
	pValue: number;
}

/**
 * Welch's t-test from pre-computed group statistics.
 * Returns NaN p-value if either group has < 2 observations or zero variance.
 */
export function welchTTest(g1: GroupStats, g2: GroupStats): TTestResult {
	if (g1.n < 2 || g2.n < 2) {
		return { tStatistic: NaN, degreesOfFreedom: NaN, pValue: NaN };
	}

	const v1 = g1.variance / g1.n;
	const v2 = g2.variance / g2.n;
	const vSum = v1 + v2;

	if (vSum === 0) {
		if (g1.mean === g2.mean) {
			return { tStatistic: 0, degreesOfFreedom: g1.n + g2.n - 2, pValue: 1 };
		}
		// Perfect separation: means differ but zero variance in both groups
		const sign = g1.mean > g2.mean ? 1 : -1;
		return { tStatistic: sign * Infinity, degreesOfFreedom: g1.n + g2.n - 2, pValue: 0 };
	}

	const t = (g1.mean - g2.mean) / Math.sqrt(vSum);
	const df = (vSum * vSum) / (v1 * v1 / (g1.n - 1) + v2 * v2 / (g2.n - 1));
	const pValue = tDistPValue(t, df);

	return { tStatistic: t, degreesOfFreedom: df, pValue };
}

// ── Effect Size / Fold Change ────────────────────────────────────────────────

const PSEUDOCOUNT = 1e-9;

/**
 * Log2 fold change matching Scanpy's formula for log1p-normalized data:
 *   log2((expm1(meanTarget) + 1e-9) / (expm1(meanRef) + 1e-9))
 */
export function log2FoldChange(meanTarget: number, meanReference: number): number {
	const numFC = Math.expm1(meanTarget) + PSEUDOCOUNT;
	const denFC = Math.expm1(meanReference) + PSEUDOCOUNT;
	return Math.log2(numFC / denFC);
}

/**
 * Adaptive effect size that handles both log1p-normalized and z-scored data.
 *
 * - Log1p data (all non-negative): Scanpy's expm1-based log2 fold change
 * - Z-scored data (has negatives): simple mean difference (target - reference)
 */
export function computeEffectSize(meanTarget: number, meanReference: number, dataIsLog1p: boolean): number {
	if (dataIsLog1p) {
		return log2FoldChange(meanTarget, meanReference);
	}
	return meanTarget - meanReference;
}

/**
 * Detect whether expression data is log1p-normalized (all >= 0) or z-scored
 * (has negative values) by sampling values from a Float32Array.
 */
export function detectDataIsLog1p(values: Float32Array): boolean {
	for (let i = 0; i < values.length; i++) {
		if (values[i] < 0) return false;
	}
	return true;
}

// ── Benjamini-Hochberg Correction ───────────────────────────────────────────

/**
 * Benjamini-Hochberg multiple testing correction.
 * Returns adjusted p-values in the same order as input.
 */
export function benjaminiHochberg(pvalues: number[]): number[] {
	const n = pvalues.length;
	if (n === 0) return [];

	const indexed = pvalues.map((p, i) => ({ p, i }));
	indexed.sort((a, b) => a.p - b.p);

	const adjusted = new Array<number>(n);
	let cumMin = 1;

	for (let rank = n - 1; rank >= 0; rank--) {
		const corrected = (indexed[rank].p * n) / (rank + 1);
		cumMin = Math.min(cumMin, corrected);
		adjusted[indexed[rank].i] = Math.min(1, cumMin);
	}

	return adjusted;
}

// ── Per-Gene DGE Computation ────────────────────────────────────────────────

export interface GeneResult {
	gene: string;
	effectSize: number;
	pval: number;
	meanTarget: number;
	meanReference: number;
	tStatistic: number;
	degreesOfFreedom: number;
}

/**
 * Compute DGE statistics for a single gene from its expression values.
 *
 * @param geneName        Gene identifier
 * @param values          Float32Array of expression values
 * @param groupAssignments Uint8Array of group indices per cell
 * @param filterArray     Uint8Array where 0 = visible, >0 = filtered out
 * @param targetGroup     Index of the target group
 * @param referenceGroup  Index of reference group, or -1 for "rest"
 * @param dataIsLog1p     If true, use Scanpy log2FC; if false, use mean difference
 */
export function computeGeneStats(
	geneName: string,
	values: Float32Array,
	groupAssignments: Uint8Array,
	filterArray: Uint8Array,
	targetGroup: number,
	referenceGroup: number,
	dataIsLog1p = true,
): GeneResult {
	const accTarget = welfordCreate();
	const accRef = welfordCreate();
	const isRest = referenceGroup < 0;

	for (let i = 0; i < values.length; i++) {
		if (filterArray[i] !== 0) continue;
		const val = values[i];
		if (Number.isNaN(val)) continue;

		const group = groupAssignments[i];
		if (group === targetGroup) {
			welfordAccumulate(accTarget, val);
		} else if (isRest || group === referenceGroup) {
			welfordAccumulate(accRef, val);
		}
	}

	const statsTarget = welfordFinalize(accTarget);
	const statsRef = welfordFinalize(accRef);
	const ttest = welchTTest(statsTarget, statsRef);

	return {
		gene: geneName,
		effectSize: computeEffectSize(statsTarget.mean, statsRef.mean, dataIsLog1p),
		pval: ttest.pValue,
		meanTarget: statsTarget.mean,
		meanReference: statsRef.mean,
		tStatistic: ttest.tStatistic,
		degreesOfFreedom: ttest.degreesOfFreedom,
	};
}
