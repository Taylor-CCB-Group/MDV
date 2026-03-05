/**
 * mdv-dge: Orchestrator that bridges MDV's DataStore column-loading
 * infrastructure with the DGE Web Worker.
 *
 * Loads gene expression columns in batches, posts them to dgeWorker,
 * collects results, applies BH correction, and emits the final table.
 */

import { benjaminiHochberg } from "./dgeStats";
import type { GeneResult } from "./dgeStats";
import type { DGEBatchInput, DGEBatchOutput } from "./dgeWorker";
import { processBatch } from "./dgeWorker";

export interface DGEConfig {
	groupColumn: string;
	targetGroup: string;
	referenceGroup: string | "rest";
	batchSize?: number;
}

export interface DGEFullResult extends GeneResult {
	pvalAdj: number;
	negLog10Pval: number;
	negLog10PvalAdj: number;
}

export type EffectSizeLabel = "log2fc" | "mean_diff";

export interface DGERunResult {
	config: DGEConfig;
	results: DGEFullResult[];
	elapsed: number;
	effectSizeLabel: EffectSizeLabel;
	skippedGenes: number;
}

type ColumnLoader = (columns: string[]) => Promise<void>;
type GetColumnBuffer = (field: string) => SharedArrayBuffer | null;

/**
 * Standalone DGE runner that works with any DataStore-like data source.
 *
 * Rather than extending Dimension (which is for filtering), this is a
 * utility class that orchestrates the full DGE pipeline:
 *   1. Resolve gene column field names from the rows_as_columns link
 *   2. Load gene columns in batches via the provided loader
 *   3. Run DGE computation (in-thread or via Worker)
 *   4. Apply BH correction
 *   5. Return the full results table
 */
export class DGERunner {
	private worker: Worker | null = null;

	/**
	 * Run DGE analysis.
	 *
	 * @param config           DGE configuration (group column, target, reference)
	 * @param filterBuffer     DataStore.filterBuffer (SharedArrayBuffer)
	 * @param size             Number of cells (DataStore.size)
	 * @param groupBuffer      SharedArrayBuffer for the group column (Uint8Array-typed)
	 * @param groupValues      Category labels for the group column
	 * @param geneFields       Array of gene column field names (gs|NAME (scores)|INDEX format)
	 * @param geneNames        Array of human-readable gene names (parallel to geneFields)
	 * @param loadColumns      Async function that loads columns into the DataStore
	 * @param getColumnBuffer  Function that returns the SharedArrayBuffer for a loaded column
	 * @param onProgress       Optional progress callback (batchesDone, totalBatches)
	 * @param useWorker        Whether to use a Web Worker (true) or compute in-thread (false)
	 * @param dataIsLog1p      Whether data is log1p-normalized (true) or z-scored (false)
	 */
	async run(
		config: DGEConfig,
		filterBuffer: SharedArrayBuffer,
		size: number,
		groupBuffer: SharedArrayBuffer,
		groupValues: string[],
		geneFields: string[],
		geneNames: string[],
		loadColumns: ColumnLoader,
		getColumnBuffer: GetColumnBuffer,
		onProgress?: (done: number, total: number) => void,
		useWorker = false,
		dataIsLog1p = true,
	): Promise<DGERunResult> {
		const t0 = performance.now();
		const batchSize = config.batchSize ?? 2000;

		const targetIdx = groupValues.indexOf(config.targetGroup);
		if (targetIdx < 0) {
			throw new Error(`Target group "${config.targetGroup}" not found in group values: ${groupValues.join(", ")}`);
		}
		let referenceIdx = -1;
		if (config.referenceGroup !== "rest") {
			referenceIdx = groupValues.indexOf(config.referenceGroup);
			if (referenceIdx < 0) {
				throw new Error(`Reference group "${config.referenceGroup}" not found in group values: ${groupValues.join(", ")}`);
			}
		}

			const allBatchResults: GeneResult[] = [];
		const totalBatches = Math.ceil(geneFields.length / batchSize);
		let skippedGenes = 0;

		for (let b = 0; b < totalBatches; b++) {
			const start = b * batchSize;
			const end = Math.min(start + batchSize, geneFields.length);
			const batchFields = geneFields.slice(start, end);
			const batchNames = geneNames.slice(start, end);

			await loadColumns(batchFields);

			const loadedNames: string[] = [];
			const geneBuffers: SharedArrayBuffer[] = [];
			for (let i = 0; i < batchFields.length; i++) {
				const buf = getColumnBuffer(batchFields[i]);
				if (!buf) {
					skippedGenes++;
					continue;
				}
				loadedNames.push(batchNames[i]);
				geneBuffers.push(buf);
			}

			if (geneBuffers.length === 0) {
				onProgress?.(b + 1, totalBatches);
				continue;
			}

			const input: DGEBatchInput = {
				type: "run_batch",
				filterBuffer,
				groupBuffer,
				targetGroup: targetIdx,
				referenceGroup: referenceIdx,
				geneBuffers,
				geneNames: loadedNames,
				batchIndex: b,
				dataIsLog1p,
			};

			let batchResult: DGEBatchOutput;
			if (useWorker) {
				batchResult = await this.runInWorker(input);
			} else {
				batchResult = processBatch(input);
			}

			allBatchResults.push(...batchResult.results);

			if (b === 0) {
				const top5 = [...batchResult.results]
					.filter(r => !Number.isNaN(r.pval))
					.sort((a, b) => a.pval - b.pval)
					.slice(0, 5);
				console.log(`[DGE diag] Batch 0 results (${batchResult.results.length} genes, ${loadedNames.length} loaded):`);
				for (const r of top5) {
					console.log(`[DGE diag]   ${r.gene}: effectSize=${r.effectSize.toFixed(4)}, pval=${r.pval.toExponential(3)}, meanTarget=${r.meanTarget.toFixed(4)}, meanRef=${r.meanReference.toFixed(4)}`);
				}
				const allPvals = batchResult.results.map(r => r.pval).filter(p => !Number.isNaN(p));
				const sigCount = allPvals.filter(p => p < 0.05).length;
				console.log(`[DGE diag]   Batch 0 sig (p<0.05): ${sigCount}/${allPvals.length}, NaN p-values: ${batchResult.results.length - allPvals.length}`);
			}

			onProgress?.(b + 1, totalBatches);
		}

		if (skippedGenes > 0) {
			console.warn(`DGE: skipped ${skippedGenes} genes due to failed column loads`);
		}

		const pvals = allBatchResults.map((r) => r.pval);
		const adjustedPvals = benjaminiHochberg(pvals);

		const results: DGEFullResult[] = allBatchResults.map((r, i) => ({
			...r,
			pvalAdj: adjustedPvals[i],
			negLog10Pval: r.pval > 0 ? -Math.log10(r.pval) : 300,
			negLog10PvalAdj: adjustedPvals[i] > 0 ? -Math.log10(adjustedPvals[i]) : 300,
		}));

		results.sort((a, b) => a.pval - b.pval);

		return {
			config,
			results,
			elapsed: performance.now() - t0,
			effectSizeLabel: dataIsLog1p ? "log2fc" : "mean_diff",
			skippedGenes,
		};
	}

	private runInWorker(input: DGEBatchInput): Promise<DGEBatchOutput> {
		if (!this.worker) {
			this.worker = new Worker(
				new URL("./dgeWorker.ts", import.meta.url),
			);
		}
		const w = this.worker;
		return new Promise((resolve, reject) => {
			w.onmessage = (e: MessageEvent<DGEBatchOutput>) => {
				resolve(e.data);
			};
			w.onerror = (e) => {
				reject(new Error(`DGE Worker error: ${e.message}`));
			};
			w.postMessage(input);
		});
	}

	destroy(): void {
		if (this.worker) {
			this.worker.terminate();
			this.worker = null;
		}
	}
}
