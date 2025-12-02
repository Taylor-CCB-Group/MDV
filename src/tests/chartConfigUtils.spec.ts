import { describe, test, expect, vi, beforeEach } from 'vitest';
import {
  deserialiseParam,
  initialiseChartConfig,
  getConcreteFieldNames
} from '../charts/chartConfigUtils';
import { RowsAsColsQuery } from '../links/link_utils';
import type DataStore from '../datastore/DataStore';
import type BaseChartT from '../charts/BaseChart';
import type { BaseConfig } from '../charts/BaseChart';
type BaseChart = BaseChartT<BaseConfig>

// Mock the RowsAsColsQuery static method
vi.mock('../links/link_utils', async () => {
  const actual = await vi.importActual('../links/link_utils');
  return {
    ...actual,
    RowsAsColsQuery: {
      fromSerialized: vi.fn(),
    }
  };
});

describe('chartConfigUtils - deserialisation', () => {
  let mockDataStore: any;
  let mockChart: any;

  beforeEach(() => {
    vi.clearAllMocks();

    mockDataStore = {
      name: 'test-datastore',
      columnIndex: {
        'col1': { name: 'col1', datatype: 'double' },
        'col2': { name: 'col2', datatype: 'text' },
      },
    };

    mockChart = {
      dataStore: mockDataStore,
      deferredInit: vi.fn((fn) => fn()), // Execute immediately for tests
      setTitle: vi.fn(),
      colorByColumn: vi.fn(),
      setToolTipColumn: vi.fn(),
      mobxAutorun: vi.fn(),
    };
  });

  describe('deserialiseParam', () => {
    test('returns string unchanged', () => {
      const result = deserialiseParam(mockDataStore as DataStore, 'col1');
      expect(result).toBe('col1');
    });

    test('deserializes RowsAsColsQuery object', () => {
      const mockQuery = {
        linkedDsName: 'genes',
        maxItems: 10,
        type: 'RowsAsColsQuery' as const
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 10,
        fields: ['gene1', 'gene2'],
      };

      (RowsAsColsQuery.fromSerialized as any).mockReturnValue(mockQueryInstance);

      const result = deserialiseParam(mockDataStore as DataStore, mockQuery);

      expect(RowsAsColsQuery.fromSerialized).toHaveBeenCalledWith(mockDataStore, mockQuery);
      expect(result).toBe(mockQueryInstance);
    });

    test('throws on invalid param', () => {
      (RowsAsColsQuery.fromSerialized as any).mockReturnValue(null);

      expect(() => {
        deserialiseParam(mockDataStore as DataStore, { type: 'RowsAsColsQuery' } as any);
      }).toThrow('Failed to deserialise param');
    });
  });

  describe('getConcreteFieldNames', () => {
    test('handles string input', () => {
      expect(getConcreteFieldNames('col1')).toEqual(['col1']);
    });

    test('handles object with fields array', () => {
      const fieldSpec = { fields: ['col1', 'col2'] };
      expect(getConcreteFieldNames(fieldSpec as any)).toEqual(['col1', 'col2']);
    });

    test('handles array of strings', () => {
      expect(getConcreteFieldNames(['col1', 'col2'])).toEqual(['col1', 'col2']);
    });

    test('handles array with mixed strings and objects', () => {
      const fieldSpecs = ['col1', { fields: ['col2', 'col3'] }];
      expect(getConcreteFieldNames(fieldSpecs as any)).toEqual(['col1', 'col2', 'col3']);
    });
  });

  describe('initialiseChartConfig - recursive deserialization', () => {
    test('deserializes primitives unchanged', () => {
      const config = {
        id: 'test-id',
        title: 'Test Chart',
        type: 'box_plot',
        param: ['col1'],
        size: [400, 300],
      };

      const result = initialiseChartConfig(config as any, mockChart as BaseChart);

      expect(result.title).toBe('Test Chart');
      expect(result.type).toBe('box_plot');
      expect(result.param).toEqual(['col1']);
    });

    test('recursively deserializes nested RowsAsColsQuery in param array', () => {
      const mockQuerySerialized = {
        linkedDsName: 'genes',
        maxItems: 5,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 5,
        fields: ['gene1', 'gene2'],
      };

      (RowsAsColsQuery.fromSerialized as any).mockReturnValue(mockQueryInstance);

      const config = {
        id: 'test-id',
        title: 'Test Chart',
        type: 'box_plot',
        param: ['col1', mockQuerySerialized],
        size: [400, 300],
      };

      const result = initialiseChartConfig(config as any, mockChart as BaseChart);

      expect(result.param).toHaveLength(2);
      expect(result.param[0]).toBe('col1');
      expect(result.param[1]).toBe(mockQueryInstance);
      expect(RowsAsColsQuery.fromSerialized).toHaveBeenCalledWith(mockDataStore, mockQuerySerialized);
    });

    test('handles deeply nested RowsAsColsQuery objects', () => {
      const mockQuerySerialized = {
        linkedDsName: 'genes',
        maxItems: 3,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 3,
        fields: ['gene1'],
      };

      (RowsAsColsQuery.fromSerialized as any).mockReturnValue(mockQueryInstance);

      const config = {
        id: 'test-id',
        title: 'Test Chart',
        type: 'custom_chart',
        nested: {
          level1: {
            level2: {
              query: mockQuerySerialized,
            },
          },
        },
        param: [],
        size: [400, 300],
      };

      const result = initialiseChartConfig(config as any, mockChart as BaseChart);

      expect(result.nested.level1.level2.query).toBe(mockQueryInstance);
    });

    test('handles color_by deserialization and calls colorByColumn', () => {
      const mockQuerySerialized = {
        linkedDsName: 'genes',
        maxItems: 1,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 1,
        fields: ['gene1'],
      };

      (RowsAsColsQuery.fromSerialized as any).mockReturnValue(mockQueryInstance);

      const config = {
        id: 'test-id',
        title: 'Test Chart',
        type: 'scatter',
        param: [],
        color_by: mockQuerySerialized,
        size: [400, 300],
      };

      const result = initialiseChartConfig(config as any, mockChart as BaseChart);

      expect(result.color_by).toBeUndefined();
      expect(mockChart.deferredInit).toHaveBeenCalled();
      expect(mockChart.colorByColumn).toHaveBeenCalledWith(mockQueryInstance);
    });

    test('handles tooltip.column array deserialization', () => {
      const mockQuerySerialized = {
        linkedDsName: 'genes',
        maxItems: 2,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 2,
        fields: ['gene1', 'gene2'],
      };

      (RowsAsColsQuery.fromSerialized as any).mockReturnValue(mockQueryInstance);

      const config = {
        id: 'test-id',
        title: 'Test Chart',
        type: 'scatter',
        param: [],
        tooltip: {
          show: true,
          column: [mockQuerySerialized],
        },
        size: [400, 300],
      };

      const result = initialiseChartConfig(config as any, mockChart as BaseChart);

      expect(result.tooltip.column).toEqual(['gene1', 'gene2']);
    });

    test('handles tooltip.column single value and calls setToolTipColumn', () => {
      const mockQuerySerialized = {
        linkedDsName: 'genes',
        maxItems: 1,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 1,
        fields: ['gene1'],
      };

      (RowsAsColsQuery.fromSerialized as any).mockReturnValue(mockQueryInstance);

      const config = {
        id: 'test-id',
        title: 'Test Chart',
        type: 'scatter',
        param: [],
        tooltip: {
          show: true,
          column: mockQuerySerialized,
        },
        size: [400, 300],
      };

      const result = initialiseChartConfig(config as any, mockChart as BaseChart);

      expect(result.tooltip.column).toBe('gene1');
      expect(mockChart.setToolTipColumn).toHaveBeenCalledWith(mockQueryInstance);
    });
  });

  describe('config round-trip serialization', () => {
    test('config can be serialized and deserialized correctly', () => {
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 5,
        toJSON: () => ({
          linkedDsName: 'genes',
          maxItems: 5,
          type: 'RowsAsColsQuery' as const,
        }),
        fields: ['gene1', 'gene2'],
      };

      // Simulate serialization
      const serializedConfig = JSON.stringify({
        id: 'test-id',
        title: 'Test Chart',
        type: 'box_plot',
        param: ['col1', mockQueryInstance],
        size: [400, 300],
      });

      const parsedConfig = JSON.parse(serializedConfig);

      // The serialized form should have the query in plain object form
      expect(parsedConfig.param[1]).toEqual({
        linkedDsName: 'genes',
        maxItems: 5,
        type: 'RowsAsColsQuery',
      });

      // Now deserialize it back
      (RowsAsColsQuery.fromSerialized as any).mockReturnValue(mockQueryInstance);

      const result = initialiseChartConfig(parsedConfig, mockChart as BaseChart);

      expect(result.param[1]).toBe(mockQueryInstance);
    });
  });
});