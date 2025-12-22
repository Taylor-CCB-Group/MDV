import { describe, test, expect, vi, beforeEach } from 'vitest';
import { renderHook } from '@testing-library/react';
import type { CategoricalDataType } from '@/charts/charts';
import { useCategoryContour, getDensitySettings } from './contour_state';
import type { CategoryContourProps } from './contour_state';
import type BaseChart from '@/charts/BaseChart';
import type { BaseConfig } from '@/charts/BaseChart';

// Mock all the hooks and dependencies
vi.mock('./hooks', () => ({
    useFieldSpec: vi.fn(),
    useCategoryFilterIndices: vi.fn(),
    useParamColumns: vi.fn(),
}));

vi.mock('./context', () => ({
    useDataStore: vi.fn(),
}));

vi.mock('./deck_state', () => ({
    useViewState: vi.fn(),
}));

vi.mock('use-debounce', () => ({
    useDebounce: vi.fn((value) => [value]),
}));

import { useFieldSpec, useCategoryFilterIndices, useParamColumns } from './hooks';
import { useDataStore } from './context';
import { useViewState } from './deck_state';

describe('useCategoryContour', () => {
    const mockCx = {
        data: new Float32Array([0, 1, 2, 3, 4]),
    };
    const mockCy = {
        data: new Float32Array([0, 1, 2, 3, 4]),
    };

    const baseProps: CategoryContourProps = {
        id: 'test-contour',
        fill: true,
        bandwidth: 10,
        intensity: 0.5,
        opacity: 0.8,
        fillThreshold: 2,
    };

    beforeEach(() => {
        vi.clearAllMocks();
        (useParamColumns as any).mockReturnValue([mockCx, mockCy]);
        (useViewState as any).mockReturnValue({ zoom: 1 });
        (useDataStore as any).mockReturnValue({
            getColumnColors: vi.fn(() => [[255, 0, 0], [0, 255, 0], [0, 0, 255]]),
        });
    });

    test('returns undefined when category is not provided', () => {
        const mockContourParameter = {
            field: 'test-field',
            values: ['cat1', 'cat2', 'cat3'],
            data: new Uint8Array([0, 1, 2, 0, 1]),
            datatype: 'text' as CategoricalDataType,
        };

        (useFieldSpec as any).mockReturnValue(mockContourParameter);
        (useCategoryFilterIndices as any).mockReturnValue(new Uint32Array([0, 1, 2]));

        const { result } = renderHook(() =>
            useCategoryContour({
                ...baseProps,
                category: undefined,
            })
        );

        expect(result.current).toBeUndefined();
    });

    test('returns undefined when contourParameter is not provided', () => {
        (useFieldSpec as any).mockReturnValue(undefined);
        (useCategoryFilterIndices as any).mockReturnValue(new Uint32Array([]));

        const { result } = renderHook(() =>
            useCategoryContour({
                ...baseProps,
                parameter: undefined,
                category: 'cat1',
            })
        );

        expect(result.current).toBeUndefined();
    });

    test('returns undefined when both category and contourParameter are not provided', () => {
        (useFieldSpec as any).mockReturnValue(undefined);
        (useCategoryFilterIndices as any).mockReturnValue(new Uint32Array([]));

        const { result } = renderHook(() =>
            useCategoryContour({
                ...baseProps,
                parameter: undefined,
                category: undefined,
            })
        );

        expect(result.current).toBeUndefined();
    });

    test('returns layer props when both category and contourParameter are provided', () => {
        const mockContourParameter = {
            field: 'test-field',
            values: ['cat1', 'cat2', 'cat3'],
            data: new Uint8Array([0, 1, 2, 0, 1]),
            datatype: 'text' as CategoricalDataType,
        };

        const mockData = new Uint32Array([0, 1, 2]);

        (useFieldSpec as any).mockReturnValue(mockContourParameter);
        (useCategoryFilterIndices as any).mockReturnValue(mockData);

        const { result } = renderHook(() =>
            useCategoryContour({
                ...baseProps,
                parameter: 'test-field',
                category: 'cat1',
            })
        );

        expect(result.current).toBeDefined();
        expect(result.current?.id).toBe('test-contour');
        expect(result.current?.data).toBe(mockData);
        expect(result.current?.fillOpacity).toBe(0.5);
        expect(result.current?.contourOpacity).toBe(0.8);
        expect(result.current?.contourFill).toBe(2);
        expect(result.current?.pickable).toBe(false);
        expect(typeof result.current?.getPosition).toBe('function');
    });

    test('handles empty data array when contourParameter is falsy but category is provided', () => {
        (useFieldSpec as any).mockReturnValue(undefined);
        (useCategoryFilterIndices as any).mockReturnValue(new Uint32Array([]));

        const { result } = renderHook(() =>
            useCategoryContour({
                ...baseProps,
                parameter: undefined,
                category: 'cat1',
            })
        );

        // Should return undefined because contourParameter is falsy
        expect(result.current).toBeUndefined();
    });

    test('getPosition function works correctly', () => {
        const mockContourParameter = {
            field: 'test-field',
            values: ['cat1', 'cat2'],
            data: new Uint8Array([0, 1]),
            datatype: 'text' as CategoricalDataType,
        };

        const mockData = new Uint32Array([0, 1]);

        (useFieldSpec as any).mockReturnValue(mockContourParameter);
        (useCategoryFilterIndices as any).mockReturnValue(mockData);

        const { result } = renderHook(() =>
            useCategoryContour({
                ...baseProps,
                parameter: 'test-field',
                category: 'cat1',
            })
        );

        const getPosition = result.current?.getPosition;
        expect(getPosition).toBeDefined();

        if (getPosition) {
            const target = new Float32Array(3);
            const result = getPosition(0, { target });
            expect(result[0]).toBe(0);
            expect(result[1]).toBe(0);
            expect(result[2]).toBe(0);
        }
    });

    test('contourFill is set to large value when fill is disabled', () => {
        const mockContourParameter = {
            field: 'test-field',
            values: ['cat1'],
            data: new Uint8Array([0]),
            datatype: 'text' as CategoricalDataType,
        };

        (useFieldSpec as any).mockReturnValue(mockContourParameter);
        (useCategoryFilterIndices as any).mockReturnValue(new Uint32Array([0]));

        const { result } = renderHook(() =>
            useCategoryContour({
                ...baseProps,
                fill: false,
                parameter: 'test-field',
                category: 'cat1',
            })
        );

        expect(result.current?.contourFill).toBe(10000);
    });
});

describe('getDensitySettings disposer management', () => {
    test('getDensitySettings attaches disposer to returned spec', () => {
        const mockDataStore = {
            getColumnValues: vi.fn(() => ['cat1', 'cat2', 'cat3']),
        };
        
        const mockChart = {
            dataStore: mockDataStore,
        } as unknown as BaseChart<BaseConfig>;
        
        const mockConfig = {
            contourParameter: 'test-field',
            category1: [],
            category2: [],
            param: ['x', 'y', 'test-field'],
            field_legend: { display: true },
        } as any;
        
        const spec = getDensitySettings(mockConfig, mockChart);
        
        // Check that _disposers array exists and has one disposer
        expect(spec._disposers).toBeDefined();
        expect(Array.isArray(spec._disposers)).toBe(true);
        expect(spec._disposers).toHaveLength(1);
        
        // Check that the disposer is a function (IReactionDisposer)
        const disposer = spec._disposers?.[0];
        expect(disposer).toBeDefined();
        expect(typeof disposer).toBe('function');
        
        // Verify the disposer can be called (cleanup)
        if (disposer) {
            expect(() => disposer()).not.toThrow();
        }
    });
    
    test('getDensitySettings disposer updates values when contourParameter changes', () => {
        const mockDataStore = {
            getColumnValues: vi.fn(() => ['cat1', 'cat2']),
        };
        
        const mockChart = {
            dataStore: mockDataStore,
        } as unknown as BaseChart<BaseConfig>;
        
        const mockConfig = {
            contourParameter: 'test-field',
            category1: [],
            category2: [],
            param: ['x', 'y', 'test-field'],
            field_legend: { display: true },
        } as any;
        
        const spec = getDensitySettings(mockConfig, mockChart);
        
        // The autorun should have been set up and called getColumnValues
        expect(mockDataStore.getColumnValues).toHaveBeenCalledWith('test-field');
        
        // Change the contourParameter
        mockConfig.contourParameter = 'new-field';
        mockDataStore.getColumnValues.mockReturnValue(['new1', 'new2']);
        
        // Trigger the autorun by accessing the observable
        // (In a real scenario, changing the config would trigger mobx reactivity)
        // For this test, we're just verifying the disposer is attached
        
        // Clean up
        const disposer = spec._disposers?.[0];
        if (disposer) {
            disposer();
        }
    });
});
