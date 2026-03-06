import { describe, expect, test, beforeAll } from "vitest";
import * as fs from "node:fs";
import * as path from "node:path";
import {
	welfordCreate,
	welfordAccumulate,
	welfordFinalize,
	welchTTest,
	log2FoldChange,
	benjaminiHochberg,
	computeGeneStats,
} from "@/datastore/dgeStats";

/**
 * Validates mdv-dge browser-side DGE against Scanpy reference results
 * using real PBMC3k data (100 gene subset from 2,638 cells, 9 leiden clusters).
 *
 * This test does NOT need a browser or Web Worker -- it tests the pure math
 * against Scanpy's rank_genes_groups output.
 */

interface ValidationData {
	n_cells: number;
	n_genes_total: number;
	n_genes_selected: number;
	genes: string[];
	leiden: string[];
	expression: Record<string, number[]>;
}

interface ReferenceData {
	method: string;
	groupby: string;
	reference: string;
	n_cells: number;
	n_genes: number;
	group_sizes: Record<string, number>;
	groups: Record<
		string,
		{
			genes: string[];
			scores: number[];
			log2fc: number[];
			pvals: number[];
			pvals_adj: number[];
		}
	>;
}

const FIXTURE_DIR = path.resolve(__dirname, "fixtures");

let validationData: ValidationData;
let referenceData: ReferenceData;

beforeAll(() => {
	validationData = JSON.parse(
		fs.readFileSync(path.join(FIXTURE_DIR, "dge_validation_data.json"), "utf-8"),
	);
	referenceData = JSON.parse(
		fs.readFileSync(path.join(FIXTURE_DIR, "dge_reference_ttest.json"), "utf-8"),
	);
});

/**
 * Run DGE for a specific target cluster vs rest, using the validation expression data.
 */
function runDGEForCluster(targetCluster: string) {
	const { leiden, genes, expression } = validationData;
	const nCells = leiden.length;

	const filterArray = new Uint8Array(nCells); // all zeros = no filtering
	const clusterValues = [...new Set(leiden)].sort();
	const groupAssignments = new Uint8Array(
		leiden.map((l) => clusterValues.indexOf(l)),
	);
	const targetGroupIdx = clusterValues.indexOf(targetCluster);

	const results = genes.map((gene) => {
		const values = new Float32Array(expression[gene]);
		return computeGeneStats(
			gene,
			values,
			groupAssignments,
			filterArray,
			targetGroupIdx,
			-1, // "rest"
		);
	});

	const pvals = results.map((r) => r.pval);
	const adjustedPvals = benjaminiHochberg(pvals);

	return results.map((r, i) => ({
		...r,
		pvalAdj: adjustedPvals[i],
	}));
}

describe("PBMC3k validation: mdv-dge vs Scanpy t-test", () => {
	test("fixture data loads correctly", () => {
		expect(validationData.n_cells).toBe(2638);
		expect(validationData.n_genes_selected).toBe(100);
		expect(validationData.genes).toHaveLength(100);
		expect(validationData.leiden).toHaveLength(2638);

		expect(referenceData.method).toBe("t-test");
		expect(referenceData.n_cells).toBe(2638);
	});

	const clustersToTest = ["0", "1", "2", "3"];

	for (const cluster of clustersToTest) {
		describe(`cluster ${cluster} vs rest`, () => {
			let mdvResults: ReturnType<typeof runDGEForCluster>;
			let refGroup: ReferenceData["groups"][string];

			beforeAll(() => {
				mdvResults = runDGEForCluster(cluster);
				refGroup = referenceData.groups[cluster];
			});

			test("log2FC direction agrees for genes with |log2FC| > 0.5", () => {
				let testedCount = 0;
				let agreementCount = 0;

				for (const mdvGene of mdvResults) {
					const refIdx = refGroup.genes.indexOf(mdvGene.gene);
					if (refIdx < 0) continue;
					const refLog2FC = refGroup.log2fc[refIdx];

					if (Math.abs(refLog2FC) < 0.5) continue;
					testedCount++;
					if (Math.sign(mdvGene.effectSize) === Math.sign(refLog2FC)) {
						agreementCount++;
					}
				}

				const agreementRate = testedCount > 0 ? agreementCount / testedCount : 1;
				expect(agreementRate).toBeGreaterThanOrEqual(0.95);
			});

			test("log2FC values are within tolerance for shared genes", () => {
				let testedCount = 0;
				let withinTolerance = 0;

				for (const mdvGene of mdvResults) {
					const refIdx = refGroup.genes.indexOf(mdvGene.gene);
					if (refIdx < 0) continue;
					const refLog2FC = refGroup.log2fc[refIdx];
					testedCount++;

					const absDiff = Math.abs(mdvGene.effectSize - refLog2FC);
					// Allow absolute tolerance of 0.15 (Scanpy computes on the full
					// normalized matrix while we compute on the same underlying data,
					// but Float32 precision and mean-of-log vs log-of-mean differences
					// can introduce small discrepancies)
					if (absDiff < 0.15 || absDiff / (Math.abs(refLog2FC) + 0.01) < 0.15) {
						withinTolerance++;
					}
				}

				const rate = testedCount > 0 ? withinTolerance / testedCount : 1;
				expect(rate).toBeGreaterThanOrEqual(0.9);
			});

			test("p-value ordering broadly agrees (top 20 overlap >= 60%)", () => {
				// Sort mdv results by p-value and get top 20 gene names
				const mdvTop20 = [...mdvResults]
					.filter((r) => !Number.isNaN(r.pval))
					.sort((a, b) => a.pval - b.pval)
					.slice(0, 20)
					.map((r) => r.gene);

				// Scanpy reference genes are already sorted by p-value
				// But only the 100 validation genes are in our analysis
				const validationGeneSet = new Set(validationData.genes);
				const refTop20 = refGroup.genes
					.filter((g) => validationGeneSet.has(g))
					.slice(0, 20);

				const overlap = mdvTop20.filter((g) => refTop20.includes(g)).length;
				const overlapRate = overlap / 20;

				expect(overlapRate).toBeGreaterThanOrEqual(0.6);
			});

			test("significant genes (p < 0.05) overlap with Scanpy", () => {
				const mdvSignificant = new Set(
					mdvResults.filter((r) => r.pval < 0.05).map((r) => r.gene),
				);

				const validationGeneSet = new Set(validationData.genes);
				const refSignificant = new Set(
					refGroup.genes.filter(
						(g, i) => refGroup.pvals[i] < 0.05 && validationGeneSet.has(g),
					),
				);

				if (refSignificant.size === 0) return;

				// How many of Scanpy's significant genes did we also call significant?
				let recalled = 0;
				for (const gene of refSignificant) {
					if (mdvSignificant.has(gene)) recalled++;
				}
				const recall = recalled / refSignificant.size;
				expect(recall).toBeGreaterThanOrEqual(0.7);
			});
		});
	}
});

describe("PBMC3k: per-gene sufficient statistics match", () => {
	test("Welford mean matches direct computation for a real gene", () => {
		const gene = validationData.genes[0];
		const values = validationData.expression[gene];
		const cluster0Indices: number[] = [];

		for (let i = 0; i < validationData.leiden.length; i++) {
			if (validationData.leiden[i] === "0") {
				cluster0Indices.push(i);
			}
		}

		// Direct mean
		let sum = 0;
		for (const idx of cluster0Indices) sum += values[idx];
		const directMean = sum / cluster0Indices.length;

		// Welford mean
		const acc = welfordCreate();
		for (const idx of cluster0Indices) welfordAccumulate(acc, values[idx]);
		const stats = welfordFinalize(acc);

		expect(stats.mean).toBeCloseTo(directMean, 6);
		expect(stats.n).toBe(cluster0Indices.length);
	});

	test("log2FC computed from group means matches Scanpy for a known gene", () => {
		const gene = "SNHG12"; // Known DE gene from fixtures
		if (!validationData.expression[gene]) return;

		const values = validationData.expression[gene];
		const leiden = validationData.leiden;

		// Compute group means
		const accTarget = welfordCreate();
		const accRef = welfordCreate();
		for (let i = 0; i < values.length; i++) {
			if (leiden[i] === "0") {
				welfordAccumulate(accTarget, values[i]);
			} else {
				welfordAccumulate(accRef, values[i]);
			}
		}
		const statsTarget = welfordFinalize(accTarget);
		const statsRef = welfordFinalize(accRef);

		const mdvLog2FC = log2FoldChange(statsTarget.mean, statsRef.mean);

		// Find Scanpy's log2FC for this gene in cluster 0
		const refGroup = referenceData.groups["0"];
		const refIdx = refGroup.genes.indexOf(gene);
		expect(refIdx).toBeGreaterThanOrEqual(0);
		const scanpyLog2FC = refGroup.log2fc[refIdx];

		expect(mdvLog2FC).toBeCloseTo(scanpyLog2FC, 1);
	});

	test("t-test p-value has correct order of magnitude for a known DE gene", () => {
		const gene = "SNHG12";
		if (!validationData.expression[gene]) return;

		const values = validationData.expression[gene];
		const leiden = validationData.leiden;

		const accTarget = welfordCreate();
		const accRef = welfordCreate();
		for (let i = 0; i < values.length; i++) {
			if (leiden[i] === "0") {
				welfordAccumulate(accTarget, values[i]);
			} else {
				welfordAccumulate(accRef, values[i]);
			}
		}
		const ttest = welchTTest(welfordFinalize(accTarget), welfordFinalize(accRef));

		const refGroup = referenceData.groups["0"];
		const refIdx = refGroup.genes.indexOf(gene);
		const scanpyPval = refGroup.pvals[refIdx];

		// Both should be significant
		expect(ttest.pValue).toBeLessThan(0.01);
		expect(scanpyPval).toBeLessThan(0.01);

		// P-values should be in a similar range (within a few orders of magnitude)
		if (scanpyPval > 1e-300 && ttest.pValue > 1e-300) {
			const log10Ratio = Math.abs(
				Math.log10(ttest.pValue) - Math.log10(scanpyPval),
			);
			expect(log10Ratio).toBeLessThan(5);
		}
	});
});
