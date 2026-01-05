import { expect, test, vi, beforeEach } from 'vitest';
import TagModel from '../table/TagModel';
import DataStore from '../datastore/DataStore';
import { DataModel } from '../table/DataModel';
import { loadColumn } from '@/dataloaders/DataLoaderUtil';

vi.mock('../datastore/DataStore');
vi.mock('../table/DataModel');
vi.mock('@/dataloaders/DataLoaderUtil');

beforeEach(() => {
  vi.clearAllMocks();
});
test('TagModel can be instantiated and creates a __tags column if needed', async () => {
  const mockColumnIndex: Record<string, any> = {};
  const mockDataStoreInstance = {
    size: 10,
    name: 'test-store',
    columnIndex: mockColumnIndex,
    filterArray: new Uint8Array(10),
    filterSize: 10,
    addListener: vi.fn(),
    addColumn: vi.fn((spec, data) => {
      mockColumnIndex[spec.name] = {
        ...spec,
        data: new Uint16Array(data),
        // TagModel expects the column to be loaded.
        // In a real scenario, this would be handled by `loadColumn`.
        // For the mock, we can assume it's loaded.
        // Let's also mock the necessary properties to pass `isColumnLoaded`.
        values: spec.values,
        buffer: data,
        datatype: 'multitext',
      };
    }),
    dataChanged: vi.fn(),
  };
  (DataStore as any).mockImplementation(() => mockDataStoreInstance);

  const mockDataModelInstance = {
    updateModel: vi.fn(),
    addListener: vi.fn(),
    data: new Int32Array(10),
  };
  (DataModel as any).mockImplementation(() => mockDataModelInstance);

  const tagModel = await TagModel.create(mockDataStoreInstance as unknown as DataStore);

  expect(tagModel).toBeInstanceOf(TagModel);
  expect(DataModel).toHaveBeenCalledTimes(1);
  expect(mockDataModelInstance.updateModel).toHaveBeenCalledTimes(1);
  expect(mockDataModelInstance.addListener).toHaveBeenCalledWith(
    'tag',
    expect.any(Function),
  );
  expect(mockDataStoreInstance.addColumn).toHaveBeenCalledTimes(1);
  const addColumnCall = mockDataStoreInstance.addColumn.mock.calls[0][0];
  expect(addColumnCall.name).toBe('__tags');
  expect(addColumnCall.datatype).toBe('multitext');
});

test('getTags should handle undefined values in the column', async () => {
  // this test is added to simulate the situation that we had at runtime...
  // although I'm not sure that it really gets at the root of the problem.
  const mockColumnWithUndefined = {
    name: '__tags',
    datatype: 'multitext' as const,
    values: ['a', 'b', undefined] as any,
    delimiter: ';',
    data: new Uint16Array([0, 1, 2]),
    buffer: new SharedArrayBuffer(6),
  };

  const mockColumnIndex = {
    __tags: mockColumnWithUndefined,
  };

  const mockDataStoreInstance = {
    size: 3,
    name: 'test-store',
    columnIndex: mockColumnIndex,
    filterArray: new Uint8Array(3),
    filterSize: 3,
    addListener: vi.fn(),
    addColumn: vi.fn(),
    dataChanged: vi.fn(),
  };

  (DataStore as any).mockImplementation(() => mockDataStoreInstance);
  (loadColumn as any).mockResolvedValue(mockColumnWithUndefined);

  const mockDataModelInstance = {
    updateModel: vi.fn(),
    addListener: vi.fn(),
    data: new Int32Array(3),
  };
  (DataModel as any).mockImplementation(() => mockDataModelInstance);

  const tagModel = await TagModel.create(mockDataStoreInstance as any);

  // Allow for the promise in the constructor to resolve
  await new Promise(process.nextTick);

  expect(tagModel.isReady).toBe(true);

  // With the old implementation of splitTags, this would throw.
  // With the new implementation, it should work fine.
  let tags: Set<string> = new Set();
  expect(() => {
    tags = tagModel.getTags();
  }).not.toThrow();

  expect(tags).toEqual(new Set(['a', 'b']));
});