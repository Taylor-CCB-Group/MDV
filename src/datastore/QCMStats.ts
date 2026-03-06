/**
 * QCMStats: Core statistical functions for Quadrat Correlation Matrix (QCM) analysis.
 *
 * Modeled after MuSpAn's QCM, this provides grid generation, tile assignment,
 * and SES (Standardized Effect Size) calculation using a permutation-based null model.
 */

export interface BBox {
	minX: number;
	minY: number;
	maxX: number;
	maxY: number;
}

export interface QCMConfig {
	sideLength: number;
	method: "quadrats" | "hexgrid";
	nPermutations?: number;
	alpha?: number;
}

export interface QCMResult {
	ses: number[][]; // Standardized Effect Size matrix
	pvalues: number[][]; // P-value matrix
	labels: string[]; // Ordered label names
	nCells: number;
	nTiles: number;
	elapsed: number;
}

/**
 * Generate a grid of quadrats (rectangles) covering the bounding box.
 */
export function generateQuadrats(bbox: BBox, sideLength: number) {
	const nx = Math.ceil((bbox.maxX - bbox.minX) / sideLength);
	const ny = Math.ceil((bbox.maxY - bbox.minY) / sideLength);
	return { nx, ny, nTiles: nx * ny };
}

/**
 * Assign each cell to a quadrat index. Returns -1 if cell is filtered or outside bbox.
 */
export function assignToQuadrats(
	x: Float32Array,
	y: Float32Array,
	filter: Uint8Array,
	bbox: BBox,
	sideLength: number,
	nx: number,
	ny: number,
): Int32Array {
	const n = x.length;
	const assignments = new Int32Array(n).fill(-1);

	for (let i = 0; i < n; i++) {
		if (filter[i] !== 0) continue;

		const xi = Math.floor((x[i] - bbox.minX) / sideLength);
		const yi = Math.floor((y[i] - bbox.minY) / sideLength);

		if (xi >= 0 && xi < nx && yi >= 0 && yi < ny) {
			assignments[i] = yi * nx + xi;
		}
	}
	return assignments;
}

/**
 * Compute the co-occurrence matrix: O[a][b] = number of tiles containing both label a and b.
 */
export function computeCoOccurrence(
	assignments: Int32Array,
	labels: Uint8Array,
	nLabels: number,
	nTiles: number,
): Uint32Array {
	// For each tile, we track which labels are present using a bitset or boolean array.
	// Since nLabels is typically small (e.g. < 100), a Uint8Array per tile is fine.
	const presence = new Uint8Array(nTiles * nLabels);

	for (let i = 0; i < assignments.length; i++) {
		const tileIdx = assignments[i];
		if (tileIdx === -1) continue;
		const labelIdx = labels[i];
		presence[tileIdx * nLabels + labelIdx] = 1;
	}

	const coOccur = new Uint32Array(nLabels * nLabels);
	for (let t = 0; t < nTiles; t++) {
		const offset = t * nLabels;
		for (let a = 0; a < nLabels; a++) {
			if (presence[offset + a] === 0) continue;
			for (let b = a; b < nLabels; b++) {
				if (presence[offset + b] === 1) {
					coOccur[a * nLabels + b]++;
				}
			}
		}
	}

	// Fill symmetric part
	for (let a = 0; a < nLabels; a++) {
		for (let b = 0; b < a; b++) {
			coOccur[a * nLabels + b] = coOccur[b * nLabels + a];
		}
	}

	return coOccur;
}

/**
 * Fisher-Yates shuffle for Uint8Array.
 */
export function shuffle(array: Uint8Array): Uint8Array {
	for (let i = array.length - 1; i > 0; i--) {
		const j = Math.floor(Math.random() * (i + 1));
		[array[i], array[j]] = [array[j], array[i]];
	}
	return array;
}

/**
 * Compute SES and p-values using permutations.
 */
export function computeSES(
	assignments: Int32Array,
	labels: Uint8Array,
	nLabels: number,
	nTiles: number,
	nPermutations = 100,
): { ses: number[][]; pvalues: number[][] } {
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

	const ses: number[][] = Array.from({ length: nLabels }, () => new Array(nLabels).fill(0));
	const pvalues: number[][] = Array.from({ length: nLabels }, () => new Array(nLabels).fill(0));

	for (let a = 0; a < nLabels; a++) {
		for (let b = 0; b < nLabels; b++) {
			const idx = a * nLabels + b;
			const mean = sums[idx] / nPermutations;
			const variance = (sumsSq[idx] / nPermutations) - (mean * mean);
			const sd = Math.sqrt(Math.max(0, variance));

			if (sd > 0) {
				ses[a][b] = (observed[idx] - mean) / sd;
			} else {
				ses[a][b] = 0;
			}

			// Empirical p-value
			pvalues[a][b] = (countsHigher[idx] + 1) / (nPermutations + 1);
		}
	}

	return { ses, pvalues };
}
