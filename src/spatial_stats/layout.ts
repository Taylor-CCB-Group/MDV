/**
 * Where the cross-PCF chart's parts go. The matrices share one set of row labels (left of the first
 * column) and one of column labels (above the first row), and are sized so every cell is the same, as
 * large as fits; the curves take a share of the width (or, in a tall area, of the height) that leaves
 * them readable without starving the matrices.
 */
export type PanelLayout = {
    /** Curves (+ legend) above the matrices rather than beside them. */
    stacked: boolean;
    curves: { width: number; height: number };
    legendWidth: number;
    /** Gutter for the type labels, shared by the grid's first row / column. */
    labelSize: number;
    grid: { cols: number; rows: number; cell: number; gap: number; strip: number };
};

/** Gutter for category labels: rotated column labels need about as much as row labels. */
export function labelGutter(maxLabelChars: number, fontPx = 11): number {
    return Math.round(Math.min(130, Math.max(22, maxLabelChars * fontPx * 0.62 + 10)));
}

const GAP = 10;
/** Below this a cell can't be read or hit; the chart scrolls instead. */
export const MIN_CELL = 6;
const STRIP = 22;

/** Best columns × rows and cell size for `m` matrices of `n × n` in a `width × height` area. */
export function matrixCells(m: number, n: number, width: number, height: number, labelSize: number) {
    let best = { cols: 1, rows: Math.max(1, m), cell: 0, gap: GAP, strip: STRIP };
    if (m === 0 || n === 0) return best;
    for (let cols = 1; cols <= m; cols++) {
        const rows = Math.ceil(m / cols);
        const cw = (width - labelSize - GAP * (cols - 1)) / (cols * n);
        const ch = (height - labelSize - STRIP * rows - GAP * (rows - 1)) / (rows * n);
        const cell = Math.floor(Math.min(cw, ch) * 4) / 4;
        if (cell > best.cell) best = { cols, rows, cell, gap: GAP, strip: STRIP };
    }
    return best;
}

const clamp = (v: number, lo: number, hi: number) => Math.min(hi, Math.max(lo, v));

/**
 * Lay out curves and `m` matrices of `n` types in `width × height`: both arrangements are tried and the
 * one giving the larger matrix cells wins (ties go to side by side), with the curves held within
 * readable bounds either way.
 */
export function panelLayout(width: number, height: number, m: number, n: number, labelSize: number): PanelLayout {
    // side by side: curves get more of the width when there is only one matrix to share it with
    const sideW = Math.round(m <= 1 ? clamp(width * 0.5, 280, 760) : clamp(width * 0.38, 280, 640));
    const side = matrixCells(m, n, width - sideW - GAP, height, labelSize);
    const stackH = Math.round(clamp(height * 0.38, 180, 420));
    const stack = matrixCells(m, n, width, height - stackH - GAP, labelSize);
    const stacked = width < 560 || stack.cell > side.cell;
    const curves = stacked ? { width, height: stackH } : { width: sideW, height };
    const grid = stacked ? stack : side;
    return {
        stacked,
        curves,
        legendWidth: Math.round(clamp(curves.width * 0.26, 110, 200)),
        labelSize,
        grid: { ...grid, cell: Math.max(grid.cell, MIN_CELL) },
    };
}
