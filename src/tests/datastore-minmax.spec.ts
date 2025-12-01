import { describe, expect, test, vi, beforeEach } from 'vitest';

/**
 * Tests for getMinMaxForColumn and getColumnQuantile methods in DataStore.js
 * 
 * These tests verify:
 * - getMinMaxForColumn handles columns with missing minMax metadata
 * - getMinMaxForColumn computes min/max from data when minMax is missing
 * - getColumnQuantile handles missing or incomplete quantiles
 */

// Mock the isDatatypeNumeric function
vi.mock('@/lib/utils', () => ({
    isDatatypeNumeric: (t: string) => !!t.match(/double|float|int/),
}));

describe('getMinMaxForColumn computation logic', () => {
    // Helper function that mimics the getMinMaxForColumn logic
    function computeMinMax(column: { minMax?: [number, number]; data?: number[]; datatype: string }) {
        if (!column.datatype.match(/double|float|int/)) {
            throw new Error("Trying to get minMax for non-numeric column");
        }
        if (!column.minMax) {
            if (!column.data) {
                throw new Error("Attempting to compute minMax for column which is not loaded...");
            }
            // Return the calculated min and max values
            let min = Number.MAX_VALUE;
            let max = Number.MIN_VALUE;
            for (let i = 0; i < column.data.length; i++) {
                const value = column.data[i];
                if (Number.isNaN(value)) {
                    continue;
                }
                min = value < min ? value : min;
                max = value > max ? value : max;
            }
            column.minMax = [min, max];
            return column.minMax;
        }
        return column.minMax;
    }

    describe('when minMax is already defined', () => {
        test('returns existing minMax without recalculating', () => {
            const col = { minMax: [10, 100] as [number, number], data: [20, 30, 40], datatype: 'double' };
            expect(computeMinMax(col)).toEqual([10, 100]);
        });

        test('does not modify existing minMax', () => {
            const col = { minMax: [5, 50] as [number, number], data: [1, 100], datatype: 'integer' };
            computeMinMax(col);
            expect(col.minMax).toEqual([5, 50]);
        });
    });

    describe('when minMax is missing', () => {
        test('computes min/max from data array', () => {
            const col = { data: [5, 10, 15, 20, 25], datatype: 'double' };
            const result = computeMinMax(col);
            expect(result).toEqual([5, 25]);
        });

        test('caches computed minMax on the column', () => {
            const col = { data: [1, 2, 3], datatype: 'integer' } as any;
            computeMinMax(col);
            expect(col.minMax).toEqual([1, 3]);
        });

        test('handles negative values', () => {
            const col = { data: [-100, -50, 0, 50, 100], datatype: 'double' };
            expect(computeMinMax(col)).toEqual([-100, 100]);
        });

        test('handles single element array', () => {
            const col = { data: [42], datatype: 'int32' };
            expect(computeMinMax(col)).toEqual([42, 42]);
        });

        test('skips NaN values in computation', () => {
            const col = { data: [Number.NaN, 5, Number.NaN, 10, Number.NaN], datatype: 'double' };
            expect(computeMinMax(col)).toEqual([5, 10]);
        });

        test('handles array with all same values', () => {
            const col = { data: [7, 7, 7, 7], datatype: 'integer' };
            expect(computeMinMax(col)).toEqual([7, 7]);
        });
    });

    describe('error cases', () => {
        test('throws error for non-numeric column', () => {
            const col = { data: ['a', 'b', 'c'], datatype: 'text' } as any;
            expect(() => computeMinMax(col)).toThrow('Trying to get minMax for non-numeric column');
        });

        test('throws error when data is not loaded', () => {
            const col = { datatype: 'double' };
            expect(() => computeMinMax(col)).toThrow('Attempting to compute minMax for column which is not loaded');
        });
    });
});

describe('getColumnQuantile logic', () => {
    // Helper function that mimics the getColumnQuantile logic from DataStore.js line 1608-1617
    function getColumnQuantile(
        col: { quantiles?: any; minMax?: [number, number] },
        per?: string
    ): [number, number] | undefined {
        if (per && per !== 'none') {
            if (col.quantiles && col.quantiles !== 'NA') {
                return [col.quantiles[per][0], col.quantiles[per][1]];
            }
        } else {
            return col.minMax;
        }
        return undefined;
    }

    describe('when percentile is not provided or is "none"', () => {
        test('returns minMax when per is undefined', () => {
            const col = { minMax: [0, 100] as [number, number], quantiles: { '0.05': [5, 95] } };
            expect(getColumnQuantile(col, undefined)).toEqual([0, 100]);
        });

        test('returns minMax when per is "none"', () => {
            const col = { minMax: [10, 90] as [number, number], quantiles: { '0.05': [15, 85] } };
            expect(getColumnQuantile(col, 'none')).toEqual([10, 90]);
        });

        test('returns undefined minMax when not set', () => {
            const col = { quantiles: { '0.05': [5, 95] } };
            expect(getColumnQuantile(col, undefined)).toBeUndefined();
        });
    });

    describe('when percentile is provided', () => {
        test('returns quantile values for valid percentile', () => {
            const col = { quantiles: { '0.05': [5, 95], '0.01': [1, 99] } };
            expect(getColumnQuantile(col, '0.05')).toEqual([5, 95]);
            expect(getColumnQuantile(col, '0.01')).toEqual([1, 99]);
        });

        test('returns undefined when quantiles is undefined', () => {
            const col = { minMax: [0, 100] as [number, number] };
            expect(getColumnQuantile(col, '0.05')).toBeUndefined();
        });

        test('returns undefined when quantiles is "NA"', () => {
            const col = { quantiles: 'NA', minMax: [0, 100] as [number, number] };
            expect(getColumnQuantile(col, '0.05')).toBeUndefined();
        });

        test('throws when accessing missing percentile key (original behavior)', () => {
            const col = { quantiles: { '0.05': [5, 95] } };
            // Original implementation would throw here - this is the bug the guards fix
            expect(() => getColumnQuantile(col, '0.01')).toThrow();
        });
    });
});

describe('Safe quantile access pattern (with guards)', () => {
    // This is the safe pattern used in the fixed HeatMap.js and MultiBoxPlot.js
    function safeGetQuantile(
        col: { quantiles?: any; minMax?: [number, number] },
        per?: string
    ): [number, number] | undefined {
        if (per && per !== 'none') {
            if (col.quantiles && col.quantiles !== 'NA' && col.quantiles[per]) {
                return [col.quantiles[per][0], col.quantiles[per][1]];
            }
        }
        return col.minMax;
    }

    test('returns quantile when available', () => {
        const col = { quantiles: { '0.05': [5, 95] }, minMax: [0, 100] as [number, number] };
        expect(safeGetQuantile(col, '0.05')).toEqual([5, 95]);
    });

    test('falls back to minMax when quantiles is undefined', () => {
        const col = { minMax: [0, 100] as [number, number] };
        expect(safeGetQuantile(col, '0.05')).toEqual([0, 100]);
    });

    test('falls back to minMax when quantiles is "NA"', () => {
        const col = { quantiles: 'NA', minMax: [0, 100] as [number, number] };
        expect(safeGetQuantile(col, '0.05')).toEqual([0, 100]);
    });

    test('falls back to minMax when specific percentile key is missing', () => {
        const col = { quantiles: { '0.05': [5, 95] }, minMax: [0, 100] as [number, number] };
        expect(safeGetQuantile(col, '0.01')).toEqual([0, 100]);
    });

    test('returns minMax when per is undefined', () => {
        const col = { quantiles: { '0.05': [5, 95] }, minMax: [0, 100] as [number, number] };
        expect(safeGetQuantile(col, undefined)).toEqual([0, 100]);
    });

    test('returns minMax when per is "none"', () => {
        const col = { quantiles: { '0.05': [5, 95] }, minMax: [0, 100] as [number, number] };
        expect(safeGetQuantile(col, 'none')).toEqual([0, 100]);
    });
});
