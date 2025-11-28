import { describe, test, expect, vi, beforeEach } from 'vitest';
import { initialiseChartConfig } from '../charts/chartConfigUtils';
import { RowsAsColsQuery } from '../links/link_utils';

describe('chartConfigUtils - standalone methodsUsingColumns (no setParams redundancy)', () => {
  let mockDataStore: any;
  let mockChart: any;

  beforeEach(() => {
    vi.clearAllMocks();

    mockDataStore = {
      name: 'test-datastore',
      columnIndex: {
        'col1': { name: 'col1', datatype: 'text' },
        'col2': { name: 'col2', datatype: 'double' },
      },
    };

    mockChart = {
      dataStore: mockDataStore,
      deferredInit: vi.fn((fn) => fn()),
      setTitle: vi.fn(),
      mobxAutorun: vi.fn(),
    };
  });

  describe('ImageTableChart scenario - methods only called from UI', () => {
    test('deserializes image_label with RowsAsColsQuery', () => {
      const mockQuerySerialized = {
        linkedDsName: 'genes',
        maxItems: 1,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 1,
        fields: ['gene_label'],
      };

      vi.spyOn(RowsAsColsQuery, 'fromSerialized').mockReturnValue(mockQueryInstance as any);

      const config = {
        id: 'test-image-table',
        title: 'Image Table',
        type: 'image_table',
        param: ['id_col'],
        image_label: mockQuerySerialized,
        sortBy: 'col2',
        image_title: 'col1',
        size: [400, 300] as [number, number],
      };

      const result = initialiseChartConfig(config as any, mockChart as any);

      // The deserialized config should have the live query object
      expect(result.image_label).toBe(mockQueryInstance);
      expect(result.sortBy).toBe('col2');
      expect(result.image_title).toBe('col1');
    });

    test('deserializes all three ImageTableChart column configs', () => {
      const mockQuery1 = {
        linkedDsName: 'genes',
        maxItems: 1,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQuery2 = {
        linkedDsName: 'proteins',
        maxItems: 1,
        type: 'RowsAsColsQuery' as const,
      };

      const mockInstance1 = { linkedDsName: 'genes', maxItems: 1, fields: ['gene1'] };
      const mockInstance2 = { linkedDsName: 'proteins', maxItems: 1, fields: ['prot1'] };

      vi.spyOn(RowsAsColsQuery, 'fromSerialized')
        .mockReturnValueOnce(mockInstance1 as any)
        .mockReturnValueOnce(mockInstance2 as any);

      const config = {
        id: 'test-image-table',
        title: 'Image Table',
        type: 'image_table',
        param: ['id_col'],
        image_label: mockQuery1,
        sortBy: mockQuery2,
        image_title: 'col1', // Mixed: query + string
        size: [400, 300] as [number, number],
      };

      const result = initialiseChartConfig(config as any, mockChart as any);

      expect(result.image_label).toBe(mockInstance1);
      expect(result.sortBy).toBe(mockInstance2);
      expect(result.image_title).toBe('col1');
    });
  });

  describe('GenomeBrowser scenario - constructor calls method if config present', () => {
    test('deserializes feature_label for GenomeBrowser', () => {
      const mockQuerySerialized = {
        linkedDsName: 'genes',
        maxItems: 1,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 1,
        fields: ['gene_name'],
      };

      vi.spyOn(RowsAsColsQuery, 'fromSerialized').mockReturnValue(mockQueryInstance as any);

      const config = {
        id: 'test-genome-browser',
        title: 'Genome Browser',
        type: 'genome_browser',
        param: ['chr', 'start', 'end'],
        feature_label: mockQuerySerialized,
        size: [800, 400] as [number, number],
      };

      const result = initialiseChartConfig(config as any, mockChart as any);

      expect(result.feature_label).toBe(mockQueryInstance);
    });

    test('handles nested config properties for GenomeBrowser', () => {
      const mockQuerySerialized = {
        linkedDsName: 'genes',
        maxItems: 1,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 1,
        fields: ['present_col'],
      };

      vi.spyOn(RowsAsColsQuery, 'fromSerialized').mockReturnValue(mockQueryInstance as any);

      const config = {
        id: 'test-genome-browser',
        title: 'Genome Browser',
        type: 'genome_browser',
        param: ['chr', 'start', 'end'],
        feature_present_column: mockQuerySerialized,
        size: [800, 400] as [number, number],
      };

      const result = initialiseChartConfig(config as any, mockChart as any);

      expect(result.feature_present_column).toBe(mockQueryInstance);
    });
  });

  describe('CellRadialChart scenario - settings-only method', () => {
    test('deserializes link thickness column config', () => {
      const mockQuerySerialized = {
        linkedDsName: 'interactions',
        maxItems: 1,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'interactions',
        maxItems: 1,
        fields: ['weight_col'],
      };

      vi.spyOn(RowsAsColsQuery, 'fromSerialized').mockReturnValue(mockQueryInstance as any);

      const config = {
        id: 'test-cell-radial',
        title: 'Cell Radial Chart',
        type: 'cell_radial_chart',
        param: ['pivot', 'interaction1', 'interaction2', 'pivot', mockQuerySerialized],
        size: [600, 600] as [number, number],
      };

      const result = initialiseChartConfig(config as any, mockChart as any);

      // The param array should have the deserialized query
      expect(result.param[4]).toBe(mockQueryInstance);
    });
  });

  describe('Mixed array/object scenarios in configEntriesUsingColumns', () => {
    test('handles arrays within config entries', () => {
      const mockQuery1 = {
        linkedDsName: 'genes',
        maxItems: 2,
        type: 'RowsAsColsQuery' as const,
      };

      const mockInstance1 = {
        linkedDsName: 'genes',
        maxItems: 2,
        fields: ['gene1', 'gene2'],
      };

      vi.spyOn(RowsAsColsQuery, 'fromSerialized').mockReturnValue(mockInstance1 as any);

      const config = {
        id: 'test-chart',
        title: 'Test Chart',
        type: 'custom',
        param: [],
        custom_columns: [mockQuery1, 'col1', mockQuery1],
        size: [400, 300] as [number, number],
      };

      const result = initialiseChartConfig(config as any, mockChart as any);

      expect(result.custom_columns).toHaveLength(3);
      expect(result.custom_columns[0]).toBe(mockInstance1);
      expect(result.custom_columns[1]).toBe('col1');
      expect(result.custom_columns[2]).toBe(mockInstance1);
    });
  });

  describe('Deeply nested config properties', () => {
    test('deserializes queries in nested objects', () => {
      const mockQuerySerialized = {
        linkedDsName: 'genes',
        maxItems: 1,
        type: 'RowsAsColsQuery' as const,
      };
      const mockQueryInstance = {
        linkedDsName: 'genes',
        maxItems: 1,
        fields: ['nested_gene'],
      };

      vi.spyOn(RowsAsColsQuery, 'fromSerialized').mockReturnValue(mockQueryInstance as any);

      const config = {
        id: 'test-chart',
        title: 'Test Chart',
        type: 'custom',
        param: [],
        settings: {
          advanced: {
            column_mapping: {
              special_column: mockQuerySerialized,
            },
          },
        },
        size: [400, 300] as [number, number],
      };

      const result = initialiseChartConfig(config as any, mockChart as any);

      expect(result.settings.advanced.column_mapping.special_column).toBe(mockQueryInstance);
    });
  });
});