import * as hdf5 from 'h5wasm';

type GroupName = 'uns' | 'obs' | 'var' | 'X' | 'layers' | 'obsm' | 'varm' | 'obsp' | 'varp';
type H5DataType = number | string | boolean;
type DatasetValue = H5DataType[] | TypedArray | null;
type MatrixValue = number[][] | TypedArray | null;

type TypedArray = Int8Array | Uint8Array | Int16Array | Uint16Array | Int32Array | Uint32Array | Float32Array | Float64Array | BigInt64Array | BigUint64Array;

interface H5Metadata {
  uns: Record<string, DatasetValue>;
  obs: Record<string, DatasetValue>;
  var: Record<string, DatasetValue>;
  X: MatrixValue;
  layers: Record<string, MatrixValue>;
  obsm: Record<string, MatrixValue>;
  varm: Record<string, MatrixValue>;
  obsp: Record<string, MatrixValue>;
  varp: Record<string, MatrixValue>;
}

interface ProcessOptions {
  onProgress?: (progress: number) => void;
  chunkSize?: number;
  memoryThreshold?: number;
}

interface H5Entity {
  close: () => void;
}

interface H5Dataset extends H5Entity {
  shape: number[];
  value: TypedArray | H5DataType[];
}

interface ProgressTracker {
  current: number;
  weight: number;
  onProgress?: (progress: number) => void;
}

const CHUNK_SIZES = {
  SMALL: 1024,
  MEDIUM: 4096,
  LARGE: 16384
} as const;

const getOptimalChunkSize = (dataSize: number): number => {
  if (dataSize < 1e6) return CHUNK_SIZES.SMALL;
  if (dataSize < 1e8) return CHUNK_SIZES.MEDIUM;
  return CHUNK_SIZES.LARGE;
};

interface ExtendedPerformance extends Performance {
  memory?: {
    usedJSHeapSize: number;
    jsHeapSizeLimit: number;
  };
}

const checkMemoryUsage = (threshold: number): void => {
  const performance = window.performance as ExtendedPerformance;
  if (performance?.memory) {
    const usage = performance.memory.usedJSHeapSize / performance.memory.jsHeapSizeLimit;
    if (usage > threshold) {
      throw new Error(`Memory usage exceeded threshold: ${Math.round(usage * 100)}%`);
    }
  }
};

const withErrorContext = <T>(
  operation: () => Promise<T>,
  context: string
): Promise<T> => {
  return operation().catch(error => {
    throw new Error(`${context}: ${error.message}`);
  });
};

const isCloseable = (entity: unknown): entity is H5Entity => {
  return typeof entity === 'object' && 
         entity !== null && 
         'close' in entity && 
         typeof (entity as any).close === 'function';
};

const readMatrix = async (
  dataset: H5Dataset,
  chunkSize: number,
  tracker?: ProgressTracker
): Promise<MatrixValue> => {
  return withErrorContext(async () => {
    if (!dataset.shape || dataset.shape.length === 0) {
      return null;
    }

    const value = dataset.value;
    if (!ArrayBuffer.isView(value) && !Array.isArray(value)) {
      return null;
    }

    // Convert to 2D array if the data is a flat array
    if (dataset.shape.length === 1) {
      const result = Array.isArray(value) ? value : Array.from(new Float64Array(value.buffer));
      tracker?.onProgress?.(1);
      return [result.map(Number)];
    }

    // Handle 2D array
    const [rows, cols] = dataset.shape;
    const flatArray = Array.isArray(value) ? value : Array.from(new Float64Array(value.buffer));
    
    const matrix: number[][] = [];
    for (let i = 0; i < rows; i++) {
      matrix[i] = flatArray.slice(i * cols, (i + 1) * cols).map(Number);
      if (tracker) {
        tracker.onProgress?.((i + 1) / rows);
      }
    }

    return matrix;
  }, 'Failed to read matrix');
};

const readMatrixGroup = async (
  groupName: GroupName,
  h5File: hdf5.File,
  entities: H5Entity[],
  chunkSize: number,
  tracker: ProgressTracker
): Promise<Record<string, MatrixValue> | MatrixValue> => {
  return withErrorContext(async () => {
    const group = h5File.get(groupName);
    if (!group) return {};

    if (isCloseable(group)) {
      entities.push(group);
    }

    if (groupName === 'X') {
      if (!(group instanceof hdf5.Dataset)) return null;
      return await readMatrix(group as unknown as H5Dataset, chunkSize, tracker);
    }

    if (!(group instanceof hdf5.Group)) {
      throw new Error(`${groupName} is not a group`);
    }

    const result: Record<string, MatrixValue> = {};
    const keys = Array.from(group.keys());
    const keyProgress = 1 / keys.length;

    await Promise.all(keys.map(async (key, index) => {
      if (typeof key !== 'string') return;

      checkMemoryUsage(0.8);

      const dataset = group.get(key);
      if (!dataset || !(dataset instanceof hdf5.Dataset)) return;

      if (isCloseable(dataset)) {
        entities.push(dataset);
      }

      const datasetTracker: ProgressTracker = {
        current: tracker.current + (index * keyProgress),
        weight: keyProgress,
        onProgress: tracker.onProgress
      };

      result[key] = await readMatrix(
        dataset as unknown as H5Dataset,
        chunkSize,
        datasetTracker
      );
    }));

    return result;
  }, `Failed to read matrix group ${groupName}`);
};

const readDatasetChunked = async (
  dataset: H5Dataset,
  chunkSize: number,
  tracker?: ProgressTracker
): Promise<DatasetValue> => {
  return withErrorContext(async () => {
    if (!dataset.shape || dataset.shape.length === 0) {
      return coerceToDatasetValue(dataset.value);
    }

    const value = dataset.value;
    
    if (value instanceof BigInt64Array || value instanceof BigUint64Array) {
      const result = Array.from(value).map(Number);
      tracker?.onProgress?.(1);
      return result;
    }
    
    if (ArrayBuffer.isView(value)) {
      try {
        if (value.buffer.byteLength % 8 === 0) {
          const result = Array.from(new Float64Array(value.buffer));
          tracker?.onProgress?.(1);
          return result;
        }
        
        if (value instanceof Int8Array || value instanceof Uint8Array) {
          const result = Array.from(value).map(Number);
          tracker?.onProgress?.(1);
          return result;
        }
        
        if (value instanceof Int16Array || value instanceof Uint16Array) {
          const result = Array.from(value).map(Number);
          tracker?.onProgress?.(1);
          return result;
        }
        
        if (value instanceof Int32Array || value instanceof Uint32Array || value instanceof Float32Array) {
          const result = Array.from(value).map(Number);
          tracker?.onProgress?.(1);
          return result;
        }
        
        const result = Array.from(value).map(Number);
        tracker?.onProgress?.(1);
        return result;
      } catch (error) {
        console.warn('Failed to convert typed array:', error);
        const result = Array.from(value).map(Number);
        tracker?.onProgress?.(1);
        return result;
      }
    }
    
    return coerceToDatasetValue(value);
  }, 'Failed to read dataset');
};

function coerceToDatasetValue(value: unknown): DatasetValue {
  if (Array.isArray(value) && value.every(item => 
    typeof item === 'number' || 
    typeof item === 'string' || 
    typeof item === 'boolean'
  )) {
    return value as DatasetValue;
  }
  return null;
}

const cleanup = async (
  entities: H5Entity[], 
  virtualPath?: string, 
  FS?: typeof hdf5.FS
): Promise<void> => {
  return withErrorContext(async () => {
    const errors: Error[] = [];

    for (const entity of entities) {
      try {
        entity.close();
      } catch (error) {
        errors.push(new Error(`Failed to close entity: ${error}`));
      }
    }

    if (virtualPath && FS) {
      try {
        FS.unlink(virtualPath);
      } catch (error) {
        errors.push(new Error(`Failed to clean up virtual file: ${error}`));
      }
    }

    if (errors.length > 0) {
      throw new Error(`Cleanup failed: ${errors.join(', ')}`);
    }
  }, 'Cleanup failed');
};

const processH5File = async (
  file: File,
  options: ProcessOptions = {}
): Promise<H5Metadata> => {
  const {
    onProgress,
    memoryThreshold = 0.8
  } = options;

  const chunkSize = getOptimalChunkSize(file.size);
  const entities: H5Entity[] = [];
  const virtualPath = '/uploaded.h5';
  
  try {
    await hdf5.ready;
    const FS = await hdf5.FS;

    if (!FS) {
      throw new Error('Failed to initialize h5wasm filesystem');
    }

    const arrayBuffer = await file.arrayBuffer();
    const uint8Array = new Uint8Array(arrayBuffer);
    FS.writeFile(virtualPath, uint8Array);

    const h5File = new hdf5.File(virtualPath, 'r');
    entities.push(h5File);

    const metadata: H5Metadata = {
      uns: {},
      obs: {},
      var: {},
      X: null,
      layers: {},
      obsm: {},
      varm: {},
      obsp: {},
      varp: {}
    };

    const progressWeight = {
      uns: 0.1,
      obs: 0.1,
      var: 0.1,
      X: 0.2,
      layers: 0.1,
      obsm: 0.1,
      varm: 0.1,
      obsp: 0.1,
      varp: 0.1
    };

    let currentProgress = 0;
    const progressTrackers = Object.fromEntries(
      Object.entries(progressWeight).map(([key, weight]) => {
        const tracker = {
          current: currentProgress,
          weight,
          onProgress
        };
        currentProgress += weight;
        return [key, tracker];
      })
    );

    const readGroup = async (
      groupName: GroupName,
      target: Record<string, DatasetValue>,
      tracker: ProgressTracker
    ): Promise<void> => {
      return withErrorContext(async () => {
        const groupEntity = h5File.get(groupName);
        if (!groupEntity) return;

        if (isCloseable(groupEntity)) {
          entities.push(groupEntity);
        }

        if (!(groupEntity instanceof hdf5.Group)) {
          throw new Error(`${groupName} is not a group`);
        }

        const keys = Array.from(groupEntity.keys());
        const keyProgress = 1 / keys.length;
        
        await Promise.all(keys.map(async (key, index) => {
          if (typeof key !== 'string') return;

          checkMemoryUsage(memoryThreshold);

          const dataset = groupEntity.get(key);
          if (!dataset || !(dataset instanceof hdf5.Dataset)) return;

          if (isCloseable(dataset)) {
            entities.push(dataset);
          }
          
          const datasetTracker: ProgressTracker = {
            current: tracker.current + (index * keyProgress),
            weight: keyProgress,
            onProgress: tracker.onProgress
          };
          
          target[key] = await readDatasetChunked(
            dataset as unknown as H5Dataset,
            chunkSize,
            datasetTracker
          );
        }));
      }, `Failed to read group ${groupName}`);
    };

    await Promise.all([
      readGroup('uns', metadata.uns, progressTrackers.uns),
      readGroup('obs', metadata.obs, progressTrackers.obs),
      readGroup('var', metadata.var, progressTrackers.var)
    ]);

    const matrixData = await Promise.all([
      readMatrixGroup('X', h5File, entities, chunkSize, progressTrackers.X),
      readMatrixGroup('layers', h5File, entities, chunkSize, progressTrackers.layers),
      readMatrixGroup('obsm', h5File, entities, chunkSize, progressTrackers.obsm),
      readMatrixGroup('varm', h5File, entities, chunkSize, progressTrackers.varm),
      readMatrixGroup('obsp', h5File, entities, chunkSize, progressTrackers.obsp),
      readMatrixGroup('varp', h5File, entities, chunkSize, progressTrackers.varp)
    ]);
    
    metadata.X = matrixData[0] as MatrixValue;
    metadata.layers = matrixData[1] as Record<string, MatrixValue>;
    metadata.obsm = matrixData[2] as Record<string, MatrixValue>;
    metadata.varm = matrixData[3] as Record<string, MatrixValue>;
    metadata.obsp = matrixData[4] as Record<string, MatrixValue>;
    metadata.varp = matrixData[5] as Record<string, MatrixValue>;
    console.log(metadata);

    onProgress?.(1);
    return metadata;
  } catch (error) {
    throw new Error(`Failed to process H5 file: ${error}`);
  } finally {
    await cleanup(entities, virtualPath, await hdf5.FS);
  }
};


export default processH5File;