import { describe, expect, test } from 'vitest';

/**
 * Tests for the guard logic used in HeatMap.js and MultiBoxPlot.js
 * to safely access col.quantiles[q][0] and col.quantiles[q][1]
 * 
 * The guard pattern is:
 *   const quantile = q && col.quantiles && col.quantiles !== "NA" && col.quantiles[q];
 * 
 * This handles:
 *   - q is null/undefined/falsy (no percentile trim requested)
 *   - col.quantiles is undefined (column has no quantiles data)
 *   - col.quantiles is "NA" (quantiles marked as not available)
 *   - col.quantiles[q] is undefined (specific percentile key missing)
 */

describe('Quantiles guard pattern', () => {
    // Helper function that mimics the guard pattern used in HeatMap.js and MultiBoxPlot.js
    function getQuantileWithGuard(q: string | null | undefined, col: { quantiles?: any }) {
        const quantile = q && col.quantiles && col.quantiles !== "NA" && col.quantiles[q];
        return quantile;
    }

    describe('when q (percentile) is falsy', () => {
        test('returns falsy when q is null', () => {
            const col = { quantiles: { '0.05': [0.1, 0.9] } };
            expect(getQuantileWithGuard(null, col)).toBeFalsy();
        });

        test('returns falsy when q is undefined', () => {
            const col = { quantiles: { '0.05': [0.1, 0.9] } };
            expect(getQuantileWithGuard(undefined, col)).toBeFalsy();
        });

        test('returns falsy when q is empty string', () => {
            const col = { quantiles: { '0.05': [0.1, 0.9] } };
            expect(getQuantileWithGuard('', col)).toBeFalsy();
        });
    });

    describe('when col.quantiles is undefined', () => {
        test('returns falsy when quantiles is undefined', () => {
            const col = {};
            expect(getQuantileWithGuard('0.05', col)).toBeFalsy();
        });

        test('returns falsy even with valid q', () => {
            const col = { quantiles: undefined };
            expect(getQuantileWithGuard('0.01', col)).toBeFalsy();
        });
    });

    describe('when col.quantiles is "NA"', () => {
        test('returns falsy when quantiles is marked as NA', () => {
            const col = { quantiles: 'NA' };
            expect(getQuantileWithGuard('0.05', col)).toBeFalsy();
        });

        test('returns falsy for all percentile values', () => {
            const col = { quantiles: 'NA' };
            expect(getQuantileWithGuard('0.001', col)).toBeFalsy();
            expect(getQuantileWithGuard('0.01', col)).toBeFalsy();
            expect(getQuantileWithGuard('0.05', col)).toBeFalsy();
        });
    });

    describe('when specific percentile key is missing', () => {
        test('returns falsy when specific percentile key does not exist', () => {
            const col = { quantiles: { '0.05': [0.1, 0.9] } };
            expect(getQuantileWithGuard('0.01', col)).toBeFalsy();
        });

        test('returns falsy for missing 0.001 key', () => {
            const col = { quantiles: { '0.05': [0.1, 0.9], '0.01': [0.05, 0.95] } };
            expect(getQuantileWithGuard('0.001', col)).toBeFalsy();
        });
    });

    describe('when quantiles are properly defined', () => {
        test('returns quantile array when all conditions are met', () => {
            const col = { quantiles: { '0.05': [0.1, 0.9] } };
            const result = getQuantileWithGuard('0.05', col);
            expect(result).toEqual([0.1, 0.9]);
        });

        test('allows access to quantile values', () => {
            const col = { quantiles: { '0.05': [0.1, 0.9], '0.01': [0.05, 0.95] } };
            const result = getQuantileWithGuard('0.01', col);
            expect(result?.[0]).toBe(0.05);
            expect(result?.[1]).toBe(0.95);
        });

        test('works with 0.001 percentile', () => {
            const col = { quantiles: { '0.001': [0.001, 0.999] } };
            const result = getQuantileWithGuard('0.001', col);
            expect(result).toEqual([0.001, 0.999]);
        });
    });
});

describe('Fallback to min/max behavior', () => {
    // Mimics the complete pattern used in HeatMap.js and MultiBoxPlot.js
    function getScaleValues(
        q: string | null | undefined,
        col: { quantiles?: any },
        min: number,
        max: number
    ): [number, number] {
        const quantile = q && col.quantiles && col.quantiles !== "NA" && col.quantiles[q];
        return [
            quantile ? quantile[0] : min,
            quantile ? quantile[1] : max,
        ];
    }

    test('falls back to min/max when q is null', () => {
        const col = { quantiles: { '0.05': [0.1, 0.9] } };
        expect(getScaleValues(null, col, 0, 100)).toEqual([0, 100]);
    });

    test('falls back to min/max when quantiles is undefined', () => {
        const col = {};
        expect(getScaleValues('0.05', col, 0, 100)).toEqual([0, 100]);
    });

    test('falls back to min/max when quantiles is NA', () => {
        const col = { quantiles: 'NA' };
        expect(getScaleValues('0.05', col, 0, 100)).toEqual([0, 100]);
    });

    test('falls back to min/max when specific percentile key is missing', () => {
        const col = { quantiles: { '0.05': [0.1, 0.9] } };
        expect(getScaleValues('0.01', col, 0, 100)).toEqual([0, 100]);
    });

    test('uses quantile values when available', () => {
        const col = { quantiles: { '0.05': [10, 90] } };
        expect(getScaleValues('0.05', col, 0, 100)).toEqual([10, 90]);
    });

    test('correctly handles edge case with zero values', () => {
        const col = { quantiles: { '0.05': [0, 0] } };
        expect(getScaleValues('0.05', col, 0, 100)).toEqual([0, 0]);
    });

    test('correctly handles negative min/max values', () => {
        const col = {};
        expect(getScaleValues('0.05', col, -50, 50)).toEqual([-50, 50]);
    });
});
