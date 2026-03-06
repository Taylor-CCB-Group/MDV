/**
 * qcmWorker: Web Worker for Quadrat Correlation Matrix (QCM) analysis.
 *
 * Offloads compute-heavy co-occurrence counting and permutations from the main thread.
 */

import { computeCoOccurrence, shuffle } from "./QCMStats";

export interface QCMBatchInput {
	type: "run_qcm";
	assignments: Int32Array;
	labels: Uint8Array;
	nLabels: number;
	nTiles: number;
	nPermutations: number;
}

export interface QCMBatchOutput {
	type: "qcm_result";
	observed: Uint32Array;
	sums: Float64Array;
	sumsSq: Float64Array;
	countsHigher: Uint32Array;
}

self.onmessage = (e: MessageEvent<QCMBatchInput>) => {
	const { assignments, labels, nLabels, nTiles, nPermutations } = e.data;

	const observed = computeCoOccurrence(assignments, labels, nLabels, nTiles);
	const sums = new Float64Array(nLabels * nLabels);
	const sumsSq = new Float64Array(nLabels * nLabels);
	const countsHigher = new Uint32Array(nLabels * nLabels);

	const shuffledLabels = new Uint8Array(labels);

	for (let p = 0; p < nPermutations; p++) {
		shuffle(shuffledLabels);
		const permObserved = computeCoOccurrence(assignments, shuffledLabels, nLabels, nTiles);

		for (let i = 0; i < nLabels * nLabels; i++) {
			const val = permObserved[i];
			sums[i] += val;
			sumsSq[i] += val * val;
			if (val >= observed[i]) {
				countsHigher[i]++;
			}
		}
		
		// Optional: report progress back to main thread
		if (p % 10 === 0) {
			self.postMessage({ type: "progress", done: p, total: nPermutations });
		}
	}

	self.postMessage({
		type: "qcm_result",
		observed,
		sums,
		sumsSq,
		countsHigher,
	});
};

// Also export the process function for in-thread execution (testing)
export function processQCM(input: QCMBatchInput): QCMBatchOutput {
	const { assignments, labels, nLabels, nTiles, nPermutations } = input;

	const observed = computeCoOccurrence(assignments, labels, nLabels, nTiles);
	const sums = new Float64Array(nLabels * nLabels);
	const sumsSq = new Float64Array(nLabels * nLabels);
	const countsHigher = new Uint32Array(nLabels * nLabels);

	const shuffledLabels = new Uint8Array(labels);

	for (let p = 0; p < nPermutations; p++) {
		shuffle(shuffledLabels);
		const permObserved = computeCoOccurrence(assignments, shuffledLabels, nLabels, nTiles);

		for (let i = 0; i < nLabels * nLabels; i++) {
			const val = permObserved[i];
			sums[i] += val;
			sumsSq[i] += val * val;
			if (val >= observed[i]) {
				countsHigher[i]++;
			}
		}
	}

	return {
		type: "qcm_result",
		observed,
		sums,
		sumsSq,
		countsHigher,
	};
}
