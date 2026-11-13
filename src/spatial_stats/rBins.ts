import type { RSelection } from "./CrossPcfMatrix";

/** One r the chart shows a matrix for: a bin index, or the mean over all bins. */
export type RBin = number | "mean";

/** Band / matrix accent colours, distinct from the pair (curve) colours. */
export const R_COLOURS = ["#ffc800", "#5fd3ff", "#ff7ad9", "#a8ff5f", "#ffffff"];
export const rColour = (index: number) => R_COLOURS[index % R_COLOURS.length];

/** The configured r's, migrating a chart saved with the single `r_bin` of earlier versions. */
export function rBinsFrom(c: { r_bins?: readonly RBin[]; r_bin?: RBin }): RBin[] {
    if (c.r_bins && c.r_bins.length > 0) return [...c.r_bins];
    return c.r_bin !== undefined ? [c.r_bin] : ["mean"];
}

export function toSelection(b: RBin, nBins: number): RSelection {
    return b === "mean" ? { kind: "mean" } : { kind: "bin", k: Math.min(Math.max(0, b), nBins - 1) };
}

/** Index of the numeric r closest to bin `k`, or -1 when there is none. */
export function nearestBin(bins: readonly RBin[], k: number): number {
    let best = -1;
    let dist = Number.POSITIVE_INFINITY;
    bins.forEach((b, i) => {
        if (b === "mean") return;
        const d = Math.abs(b - k);
        if (d < dist) {
            dist = d;
            best = i;
        }
    });
    return best;
}

/** Move r `index` to bin `k` — unless another r is already there, so indices into `bins` (bands,
 *  matrices) stay one-to-one with distinct r's. */
export function moveBin(bins: readonly RBin[], index: number, k: number): RBin[] {
    if (bins.some((b, i) => b === k && i !== index)) return [...bins];
    return bins.map((b, i) => (i === index ? k : b));
}

/** `bins` plus `b`, unless it is already there. */
export function addBin(bins: readonly RBin[], b: RBin): RBin[] {
    return bins.includes(b) ? [...bins] : [...bins, b];
}

/** `bins` without entry `index`; the last one stays, so there is always a matrix. */
export function removeBin(bins: readonly RBin[], index: number): RBin[] {
    return bins.length <= 1 ? [...bins] : bins.filter((_, i) => i !== index);
}

/** Clamp numeric r's to a (possibly reduced) bin count, dropping duplicates that creates. */
export function clampBins(bins: readonly RBin[], nBins: number): RBin[] {
    const out: RBin[] = [];
    for (const b of bins) {
        const c = b === "mean" ? b : Math.min(b, nBins - 1);
        if (!out.includes(c)) out.push(c);
    }
    return out;
}

/** Columns × rows for `m` square matrices in a `width × height` area, maximising each one's side. */
export function matrixGrid(m: number, width: number, height: number): { cols: number; rows: number; side: number } {
    let best = { cols: 1, rows: Math.max(m, 1), side: 0 };
    for (let cols = 1; cols <= Math.max(m, 1); cols++) {
        const rows = Math.ceil(m / cols);
        const side = Math.min(width / cols, height / rows);
        if (side > best.side) best = { cols, rows, side };
    }
    return best;
}
