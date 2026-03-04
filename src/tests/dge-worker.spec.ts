import { describe, expect, test } from "vitest";
import { processBatch, type DGEBatchInput } from "@/datastore/dgeWorker";
import { benjaminiHochberg } from "@/datastore/dgeStats";

/**
 * Tests for the DGE worker's batch processing logic.
 * Since vitest runs in jsdom (no real Web Workers), we test
 * processBatch() directly -- it contains all the computation.
 */

function createSharedBuffer<T extends Float32ArrayConstructor | Uint8ArrayConstructor>(
	ArrayType: T,
	values: number[],
): SharedArrayBuffer {
	const buf = new SharedArrayBuffer(values.length * ArrayType.BYTES_PER_ELEMENT);
	const arr = new (ArrayType as any)(buf);
	arr.set(values);
	return buf;
}

function makeTestData(opts: {
	nCells: number;
	nGroups: number;
	genes: Record<string, number[]>;
	filter?: number[];
}) {
	const { nCells, nGroups, genes, filter } = opts;

	const groupAssignments: number[] = [];
	const groupSize = Math.floor(nCells / nGroups);
	for (let i = 0; i < nCells; i++) {
		groupAssignments.push(Math.min(Math.floor(i / groupSize), nGroups - 1));
	}

	const filterArray = filter ?? new Array(nCells).fill(0);

	return {
		filterBuffer: createSharedBuffer(Uint8Array, filterArray),
		groupBuffer: createSharedBuffer(Uint8Array, groupAssignments),
		geneBuffers: Object.values(genes).map((v) =>
			createSharedBuffer(Float32Array, v),
		),
		geneNames: Object.keys(genes),
		groupAssignments,
	};
}

describe("processBatch", () => {
	test("returns correct results for clearly different groups", () => {
		const nCells = 10;
		// Group 0 low with some variance, group 1 high with some variance
		const gene1 = [0.1, 0.3, 0.2, 0.5, 0.4, 4.8, 5.2, 5.0, 4.9, 5.1];
		const gene2 = [1.0, 1.1, 0.9, 1.0, 1.0, 1.0, 0.9, 1.1, 1.0, 1.0]; // minimal difference

		const { filterBuffer, groupBuffer, geneBuffers, geneNames } = makeTestData({
			nCells,
			nGroups: 2,
			genes: { GeneA: gene1, GeneB: gene2 },
		});

		const result = processBatch({
			type: "run_batch",
			filterBuffer,
			groupBuffer,
			targetGroup: 0,
			referenceGroup: -1,
			geneBuffers,
			geneNames,
			batchIndex: 0,
			dataIsLog1p: true,
		});

		expect(result.type).toBe("batch_result");
		expect(result.batchIndex).toBe(0);
		expect(result.results).toHaveLength(2);

		const geneA = result.results[0];
		expect(geneA.gene).toBe("GeneA");
		expect(geneA.meanTarget).toBeCloseTo(0.3, 1);
		expect(geneA.meanReference).toBeCloseTo(5.0, 1);
		expect(geneA.effectSize).toBeLessThan(0);
		expect(geneA.pval).toBeLessThan(0.001);

		const geneB = result.results[1];
		expect(geneB.gene).toBe("GeneB");
		expect(geneB.pval).toBeGreaterThan(0.05);
	});

	test("correctly handles target vs specific reference group (not rest)", () => {
		const nCells = 15;
		// 3 groups of 5: group 0 low, group 1 medium, group 2 high
		const expression = [1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2];

		const filterBuffer = createSharedBuffer(Uint8Array, new Array(nCells).fill(0));
		const groupBuffer = createSharedBuffer(Uint8Array, groups);
		const geneBuffers = [createSharedBuffer(Float32Array, expression)];

		const result = processBatch({
			type: "run_batch",
			filterBuffer,
			groupBuffer,
			targetGroup: 0,
			referenceGroup: 2,
			geneBuffers,
			geneNames: ["TestGene"],
			batchIndex: 0,
			dataIsLog1p: true,
		});

		const gene = result.results[0];
		expect(gene.meanTarget).toBeCloseTo(1, 4);
		expect(gene.meanReference).toBeCloseTo(10, 4);
		// group 1 should be excluded from the reference
	});

	test("respects filter buffer", () => {
		const nCells = 10;
		// First cell has outlier value 100, but we filter it out
		const expression = [100, 2, 3, 4, 5, 10, 11, 12, 13, 14];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const filter = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]; // filter out index 0

		const filterBuffer = createSharedBuffer(Uint8Array, filter);
		const groupBuffer = createSharedBuffer(Uint8Array, groups);
		const geneBuffers = [createSharedBuffer(Float32Array, expression)];

		const result = processBatch({
			type: "run_batch",
			filterBuffer,
			groupBuffer,
			targetGroup: 0,
			referenceGroup: -1,
			geneBuffers,
			geneNames: ["Gene"],
			batchIndex: 0,
			dataIsLog1p: true,
		});

		// With outlier filtered: group 0 = [2,3,4,5], mean = 3.5
		expect(result.results[0].meanTarget).toBeCloseTo(3.5, 4);
	});

	test("handles multiple batches correctly", () => {
		const nCells = 20;
		const groups: number[] = [];
		for (let i = 0; i < nCells; i++) groups.push(i < 10 ? 0 : 1);

		const filterBuffer = createSharedBuffer(Uint8Array, new Array(nCells).fill(0));
		const groupBuffer = createSharedBuffer(Uint8Array, groups);

		// Batch 0: strongly differentially expressed gene
		const gene1Expr: number[] = [];
		for (let i = 0; i < 10; i++) gene1Expr.push(0.5 + Math.sin(i) * 0.2);
		for (let i = 0; i < 10; i++) gene1Expr.push(5.0 + Math.sin(i) * 0.2);

		const result0 = processBatch({
			type: "run_batch",
			filterBuffer,
			groupBuffer,
			targetGroup: 0,
			referenceGroup: -1,
			geneBuffers: [createSharedBuffer(Float32Array, gene1Expr)],
			geneNames: ["Gene1"],
			batchIndex: 0,
			dataIsLog1p: true,
		});

		// Batch 1: non-differentially expressed gene
		const gene2Expr = new Array(nCells).fill(0).map((_, i) => 3.0 + Math.sin(i) * 0.1);
		const result1 = processBatch({
			type: "run_batch",
			filterBuffer,
			groupBuffer,
			targetGroup: 0,
			referenceGroup: -1,
			geneBuffers: [createSharedBuffer(Float32Array, gene2Expr)],
			geneNames: ["Gene2"],
			batchIndex: 1,
			dataIsLog1p: true,
		});

		expect(result0.batchIndex).toBe(0);
		expect(result1.batchIndex).toBe(1);

		// Combine and apply BH
		const allResults = [...result0.results, ...result1.results];
		const pvals = allResults.map((r) => r.pval);
		const padj = benjaminiHochberg(pvals);

		expect(padj).toHaveLength(2);
		expect(padj[0]).toBeLessThan(padj[1]);
	});

	test("handles gene with all-zero expression", () => {
		const nCells = 10;
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const filterBuffer = createSharedBuffer(Uint8Array, new Array(nCells).fill(0));
		const groupBuffer = createSharedBuffer(Uint8Array, groups);

		const result = processBatch({
			type: "run_batch",
			filterBuffer,
			groupBuffer,
			targetGroup: 0,
			referenceGroup: -1,
			geneBuffers: [createSharedBuffer(Float32Array, new Array(nCells).fill(0))],
			geneNames: ["ZeroGene"],
			batchIndex: 0,
			dataIsLog1p: true,
		});

		const gene = result.results[0];
		expect(gene.meanTarget).toBe(0);
		expect(gene.meanReference).toBe(0);
		expect(gene.effectSize).toBeCloseTo(0, 6);
		// p-value should be 1 (no difference) or NaN (zero variance)
		expect(gene.pval === 1 || Number.isNaN(gene.pval)).toBe(true);
	});

	test("handles many genes in a single batch", () => {
		const nCells = 20;
		const groups: number[] = [];
		for (let i = 0; i < nCells; i++) groups.push(i < 10 ? 0 : 1);

		const filterBuffer = createSharedBuffer(Uint8Array, new Array(nCells).fill(0));
		const groupBuffer = createSharedBuffer(Uint8Array, groups);

		const nGenes = 100;
		const geneBuffers: SharedArrayBuffer[] = [];
		const geneNames: string[] = [];
		for (let g = 0; g < nGenes; g++) {
			const expr: number[] = [];
			for (let i = 0; i < nCells; i++) {
				const noise = Math.sin(i * 13 + g * 7) * 0.3;
				const groupEffect = i < 10 ? 0 : g * 0.1;
				expr.push(2.0 + noise + groupEffect);
			}
			geneBuffers.push(createSharedBuffer(Float32Array, expr));
			geneNames.push(`Gene${g}`);
		}

		const result = processBatch({
			type: "run_batch",
			filterBuffer,
			groupBuffer,
			targetGroup: 0,
			referenceGroup: -1,
			geneBuffers,
			geneNames,
			batchIndex: 0,
			dataIsLog1p: true,
		});

		expect(result.results).toHaveLength(nGenes);

		// Genes with higher indices (larger effects) should rank near the top
		const sortedByPval = [...result.results]
			.filter((r) => !Number.isNaN(r.pval))
			.sort((a, b) => a.pval - b.pval);
		const topGeneIndex = Number.parseInt(sortedByPval[0].gene.replace("Gene", ""));
		expect(topGeneIndex).toBeGreaterThan(80);
	});

	test("handles NaN values in expression", () => {
		const nCells = 10;
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const expression = [1, NaN, 3, 4, 5, 10, NaN, 12, 13, 14];

		const filterBuffer = createSharedBuffer(Uint8Array, new Array(nCells).fill(0));
		const groupBuffer = createSharedBuffer(Uint8Array, groups);
		const geneBuffers = [createSharedBuffer(Float32Array, expression)];

		const result = processBatch({
			type: "run_batch",
			filterBuffer,
			groupBuffer,
			targetGroup: 0,
			referenceGroup: -1,
			geneBuffers,
			geneNames: ["NaNGene"],
			batchIndex: 0,
			dataIsLog1p: true,
		});

		const gene = result.results[0];
		// Group 0: [1, 3, 4, 5] -> mean 3.25, Group 1: [10, 12, 13, 14] -> mean 12.25
		expect(gene.meanTarget).toBeCloseTo(3.25, 3);
		expect(gene.meanReference).toBeCloseTo(12.25, 3);
		expect(gene.pval).toBeLessThan(0.01);
	});
});
