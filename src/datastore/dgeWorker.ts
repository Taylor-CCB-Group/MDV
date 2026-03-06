/**
 * mdv-dge Web Worker: processes gene expression batches and computes
 * per-gene DGE statistics using sufficient statistics (Welford + Welch t-test).
 *
 * Protocol:
 *   Main thread sends "run_batch" messages with SharedArrayBuffers.
 *   Worker responds with "batch_result" containing per-gene stats.
 *
 * All heavy computation is delegated to dgeStats.ts pure functions.
 */

import { computeGeneStats } from "./dgeStats";
import type { GeneResult } from "./dgeStats";

export interface DGEBatchInput {
	type: "run_batch";
	filterBuffer: SharedArrayBuffer;
	groupBuffer: SharedArrayBuffer;
	targetGroup: number;
	referenceGroup: number; // -1 means "rest"
	geneBuffers: SharedArrayBuffer[];
	geneNames: string[];
	batchIndex: number;
	dataType: "log1p" | "linear" | "zscored";
}

export interface DGEBatchOutput {
	type: "batch_result";
	batchIndex: number;
	results: GeneResult[];
}

/**
 * Core batch computation, exported for direct use in tests.
 */
export function processBatch(msg: DGEBatchInput): DGEBatchOutput {
	const filterArray = new Uint8Array(msg.filterBuffer);
	const groupAssignments = new Uint8Array(msg.groupBuffer);

	const results: GeneResult[] = [];
	for (let i = 0; i < msg.geneNames.length; i++) {
		const values = new Float32Array(msg.geneBuffers[i]);
		results.push(
			computeGeneStats(
				msg.geneNames[i],
				values,
				groupAssignments,
				filterArray,
				msg.targetGroup,
				msg.referenceGroup,
				msg.dataType,
			),
		);
	}

	return {
		type: "batch_result",
		batchIndex: msg.batchIndex,
		results,
	};
}

// Worker message handler -- only active when running as a Worker
// biome-ignore lint/suspicious/noGlobalAssign: Web Worker entry point
onmessage = (e: MessageEvent<DGEBatchInput>) => {
	if (e.data.type === "run_batch") {
		const result = processBatch(e.data);
		postMessage(result);
	}
};
