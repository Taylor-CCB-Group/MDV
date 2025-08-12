import * as zarr from 'zarrita';

// Types for Zarr metadata structure
export interface ZarrMetadata {
  datasetStructure: {
    groups: string[];
    arrays: string[];
  };
  spatialData?: {
    coordinateSystems: Record<string, any>;
    elements: Record<string, any>;
    transformations: any[];
  };
  images?: {
    [key: string]: {
      shape: number[];
      dtype: string;
      chunks: number[];
      dimensions: string[];
      physicalSizeX?: number;
      physicalSizeY?: number;
      physicalSizeZ?: number;
      units?: string;
      channels?: string[];
    };
  };
  tables?: {
    [key: string]: {
      shape: number[];
      obsNames: string[];
      varNames: string[];
      observations: Record<string, any>;
      variables: Record<string, any>;
      embeddings?: string[];
    };
  };
  labels?: {
    [key: string]: {
      shape: number[];
      dtype: string;
      labelValues: number[];
      categories?: string[];
    };
  };
  groupAttributes: Record<string, any>;
  rawStructure: ZarrNode;
}

export interface ZarrNode {
  name: string;
  path: string;
  type: 'group' | 'array';
  attributes?: Record<string, any>;
  shape?: number[];
  dtype?: string;
  chunks?: number[];
  children?: ZarrNode[];
  arrayMetadata?: {
    compressor?: any;
    filters?: any[];
    fill_value?: any;
    order?: string;
    dimension_separator?: string;
  };
}

/**
 * Recursively explores a Zarr group structure and extracts metadata
 */
async function exploreZarrGroup(
  group: any,
  path = "",
  visited: Set<string> = new Set()
): Promise<ZarrNode> {
  const fullPath = path || '/';
  
  // Prevent infinite recursion
  if (visited.has(fullPath)) {
    return {
      name: path.split('/').pop() || 'root',
      path: fullPath,
      type: 'group',
      children: []
    };
  }
  visited.add(fullPath);

  const node: ZarrNode = {
    name: path.split('/').pop() || 'root',
    path: fullPath,
    type: 'group',
    attributes: {},
    children: []
  };

  try {
    // Get group attributes
    if (group.attrs) {
      node.attributes = { ...group.attrs };
    }

    // Get all keys in the group
    const keys = group.keys ? await group.keys() : [];
    
    for (const key of keys) {
      try {
        const childPath = path ? `${path}/${key}` : key;
        const child = group.resolve ? group.resolve(key) : await group.getItem(key);
        
        if (!child) continue;

        // Check if it's an array or group
        try {
          const opened = await zarr.open(child);
          
          // Type guard: check if 'opened' is an array by checking for 'shape' property
          if ('shape' in opened && Array.isArray((opened as any).shape)) {
            // It's an array
            const arrayNode: ZarrNode = {
              name: key,
              path: childPath,
              type: 'array',
              shape: (opened as any).shape,
              dtype: (opened as any).dtype,
              chunks: (opened as any).chunks,
              attributes: (opened as any).attrs || {},
              arrayMetadata: {
                compressor: (opened as any).compressor,
                filters: (opened as any).filters,
                fill_value: (opened as any).fill_value,
                order: (opened as any).order,
                dimension_separator: (opened as any).dimension_separator
              }
            };
            node.children?.push(arrayNode);
          } else {
            // It's a group, recurse
            const groupNode = await exploreZarrGroup(opened, childPath, visited);
            node.children?.push(groupNode);
          }
        } catch (arrayError) {
          // If we can't open as array, try as group
          try {
            const groupNode = await exploreZarrGroup(child, childPath, visited);
            node.children?.push(groupNode);
          } catch (groupError) {
            console.warn(`Could not explore ${childPath}:`, groupError);
          }
        }
      } catch (error) {
        console.warn(`Error processing key ${key}:`, error);
      }
    }
  } catch (error) {
    console.warn(`Error exploring group at ${path}:`, error);
  }

  return node;
}

/**
 * Flattens the tree structure to extract groups and arrays
 */
function flattenZarrStructure(node: ZarrNode): { groups: string[]; arrays: string[] } {
  const groups: string[] = [];
  const arrays: string[] = [];

  function traverse(n: ZarrNode) {
    if (n.type === 'group' && n.path !== '/') {
      groups.push(n.path);
    } else if (n.type === 'array') {
      arrays.push(n.path);
    }

    n.children?.forEach(traverse);
  }

  traverse(node);
  return { groups, arrays };
}

/**
 * Extract spatial data metadata following SpatialData conventions
 */
function extractSpatialMetadata(rootNode: ZarrNode): ZarrMetadata['spatialData'] {
  const spatialData: ZarrMetadata['spatialData'] = {
    coordinateSystems: {},
    elements: {},
    transformations: []
  };

  // Look for spatial data attributes in root
  if (rootNode.attributes) {
    if (rootNode.attributes.coordinateSystems) {
      spatialData.coordinateSystems = rootNode.attributes.coordinateSystems;
    }
    if (rootNode.attributes.elements) {
      spatialData.elements = rootNode.attributes.elements;
    }
    if (rootNode.attributes.transformations) {
      spatialData.transformations = rootNode.attributes.transformations;
    }
  }

  return Object.keys(spatialData.coordinateSystems).length > 0 ||
         Object.keys(spatialData.elements).length > 0 ||
         spatialData.transformations.length > 0
    ? spatialData
    : undefined;
}

/**
 * Extract image metadata
 */
function extractImageMetadata(node: ZarrNode): ZarrMetadata['images'] {
  const images: NonNullable<ZarrMetadata['images']> = {};

  function findImages(n: ZarrNode, basePath = '') {
    if (n.type === 'array' && n.shape && n.shape.length >= 2) {
      // Likely an image if it has 2+ dimensions
      const name = basePath ? `${basePath}/${n.name}` : n.name;
      
      images[name] = {
        shape: n.shape,
        dtype: n.dtype || 'unknown',
        chunks: n.chunks || [],
        dimensions: n.attributes?.axes || n.shape.map((_, i) => `dim_${i}`),
        physicalSizeX: n.attributes?.physicalSizeX,
        physicalSizeY: n.attributes?.physicalSizeY,
        physicalSizeZ: n.attributes?.physicalSizeZ,
        units: n.attributes?.units,
        channels: n.attributes?.channels
      };
    }

    n.children?.forEach(child => {
      const childPath = basePath ? `${basePath}/${child.name}` : child.name;
      findImages(child, childPath);
    });
  }

  findImages(node);
  return Object.keys(images).length > 0 ? images : undefined;
}

/**
 * Extract table/AnnData metadata
 */
function extractTableMetadata(node: ZarrNode): ZarrMetadata['tables'] {
  const tables: NonNullable<ZarrMetadata['tables']> = {};

  function findTables(n: ZarrNode, basePath = '') {
    // Look for AnnData-like structures
    const hasObs = n.children?.some(c => c.name === 'obs');
    const hasVar = n.children?.some(c => c.name === 'var');
    const hasX = n.children?.some(c => c.name === 'X');

    if (hasObs && hasVar && hasX) {
      const name = basePath ? `${basePath}/${n.name}` : n.name;
      const obsArray = n.children?.find(c => c.name === 'obs');
      const varArray = n.children?.find(c => c.name === 'var');
      const xArray = n.children?.find(c => c.name === 'X');

      tables[name] = {
        shape: xArray?.shape || [0, 0],
        obsNames: [], // Would need to actually read the data to get these
        varNames: [], // Would need to actually read the data to get these
        observations: obsArray?.attributes || {},
        variables: varArray?.attributes || {},
        embeddings: n.children
          ?.filter(c => c.name.startsWith('X_'))
          ?.map(c => c.name) || []
      };
    }

    n.children?.forEach(child => {
      const childPath = basePath ? `${basePath}/${child.name}` : child.name;
      findTables(child, childPath);
    });
  }

  findTables(node);
  return Object.keys(tables).length > 0 ? tables : undefined;
}

/**
 * Extract label metadata
 */
function extractLabelMetadata(node: ZarrNode): ZarrMetadata['labels'] {
  const labels: NonNullable<ZarrMetadata['labels']> = {};

  function findLabels(n: ZarrNode, basePath = '') {
    // Look for integer arrays that might be labels
    if (n.type === 'array' && n.dtype && 
        (n.dtype.includes('int') || n.dtype.includes('uint')) &&
        n.shape && n.shape.length >= 2) {
      
      const name = basePath ? `${basePath}/${n.name}` : n.name;
      
      labels[name] = {
        shape: n.shape,
        dtype: n.dtype,
        labelValues: [], // Would need to read data to get actual values
        categories: n.attributes?.categories
      };
    }

    n.children?.forEach(child => {
      const childPath = basePath ? `${basePath}/${child.name}` : child.name;
      findLabels(child, childPath);
    });
  }

  findLabels(node);
  return Object.keys(labels).length > 0 ? labels : undefined;
}

/**
 * Main function to extract Zarr metadata from a URL
 */
export async function extractZarrMetadata(url: string): Promise<ZarrMetadata> {
  try {
    // Ensure URL ends with proper separator
    const normalizedUrl = url.endsWith('/') ? url.slice(0, -1) : url;
    
    // Create store and try to open with consolidated metadata first
    const store = new zarr.FetchStore(normalizedUrl);
    let root;
    
    try {
      // Try with consolidated metadata first (more efficient)
      const consolidatedStore = await zarr.tryWithConsolidated(store);
      root = await zarr.open(consolidatedStore);
    } catch (consolidatedError) {
      console.log('Consolidated metadata not available, falling back to regular access');
      // Fall back to regular store
      root = await zarr.open(store);
    }

    // Recursively explore the structure
    console.log('Exploring Zarr structure...');
    const rawStructure = await exploreZarrGroup(root);
    
    // Extract flattened structure
    const { groups, arrays } = flattenZarrStructure(rawStructure);
    
    // Extract specialized metadata
    const spatialData = extractSpatialMetadata(rawStructure);
    const images = extractImageMetadata(rawStructure);
    const tables = extractTableMetadata(rawStructure);
    const labels = extractLabelMetadata(rawStructure);

    const metadata: ZarrMetadata = {
      datasetStructure: { groups, arrays },
      groupAttributes: rawStructure.attributes || {},
      rawStructure,
      spatialData,
      images,
      tables,
      labels
    };

    console.log('Zarr metadata extraction completed:', metadata);
    return metadata;

  } catch (error) {
    console.error('Error extracting Zarr metadata:', error);
    throw new Error(`Failed to extract metadata: ${error instanceof Error ? error.message : String(error)}`);
  }
}

/**
 * Debug function to explore Zarr dataset structure via HTTP
 */
export async function debugExploreWithHTTP(url: string): Promise<void> {
  try {
    // Try to fetch the directory listing
    const response = await fetch(url);
    const html = await response.text();
    
    console.log('HTTP directory listing:', html);
    
    // Try to fetch .zgroup and .zattrs if they exist
    try {
      const zgroupResponse = await fetch(`${url}/.zgroup`);
      if (zgroupResponse.ok) {
        const zgroupData = await zgroupResponse.json();
        console.log('.zgroup content:', zgroupData);
      }
    } catch (e) {
      console.log('No .zgroup file found or not accessible');
    }

    try {
      const zattrsResponse = await fetch(`${url}/.zattrs`);
      if (zattrsResponse.ok) {
        const zattrsData = await zattrsResponse.json();
        console.log('.zattrs content:', zattrsData);
      }
    } catch (e) {
      console.log('No .zattrs file found or not accessible');
    }

  } catch (error) {
    console.error('HTTP exploration failed:', error);
  }
}

/**
 * Debug function to explore Zarr dataset using zarrita
 */
export async function debugExploreZarrDataset(url: string): Promise<void> {
  try {
    console.log(`Debugging Zarr dataset at: ${url}`);
    
    const store = new zarr.FetchStore(url);
    console.log('Store created:', store);
    
    const root = await zarr.open(store);
    console.log('Root opened:', root);
    
    // Try to get keys if root is a group
    if ((root as any).keys && typeof (root as any).keys === 'function') {
      const keys = await (root as any).keys();
      console.log('Root keys:', keys);
      
      // Explore each key
      for (const key of keys.slice(0, 3)) { // Limit to first 3 for debugging
        try {
          const item = (root as any).resolve(key);
          const opened = await zarr.open(item);
          if ('shape' in opened) {
            console.log(`Item ${key}:`, {
              shape: (opened as any).shape,
              dtype: (opened as any).dtype,
              attrs: (opened as any).attrs
            });
          } else {
            console.log(`Item ${key} is a group, not an array.`);
          }
        } catch (itemError) {
          console.log(`Could not open item ${key}:`, itemError);
        }
      }
    }

    // Try to access attributes
    console.log('Root attributes:', root.attrs);
    
  } catch (error) {
    console.error('Zarr exploration failed:', error);
  }
}

export default extractZarrMetadata;