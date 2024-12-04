import * as hdf5 from 'h5wasm';

type GroupName = keyof H5Metadata;
type H5DataType = number | string | boolean;
type DatasetValue = H5DataType[] | null;
type TypedArray = Int8Array | Uint8Array | Int16Array | Uint16Array | Int32Array | Uint32Array | Float32Array | Float64Array | BigInt64Array | BigUint64Array;

interface H5Metadata {
  uns: Record<string, DatasetValue>;
  obs: Record<string, DatasetValue>;
  var: Record<string, DatasetValue>;
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
      const result = Array.from(new Float64Array(value.buffer));
      tracker?.onProgress?.(1);
      return result;
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
    };

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

    const progressTrackers = {
      uns: { current: 0, weight: 0.33, onProgress },
      obs: { current: 0.33, weight: 0.33, onProgress },
      var: { current: 0.66, weight: 0.34, onProgress }
    };

    await Promise.all([
      readGroup('uns', metadata.uns, progressTrackers.uns),
      readGroup('obs', metadata.obs, progressTrackers.obs),
      readGroup('var', metadata.var, progressTrackers.var)
    ]);

    onProgress?.(1);
    return metadata;
  } catch (error) {
    throw new Error(`Failed to process H5 file: ${error}`);
  } finally {
    await cleanup(entities, virtualPath, await hdf5.FS);
  }
};

export default processH5File;