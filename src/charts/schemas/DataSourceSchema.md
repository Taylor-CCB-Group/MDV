# DataSource Schema

This document describes the schema for MDV datasource configurations. The schema is defined using Zod and provides comprehensive validation for all datasource configuration options.

## Overview

The datasource schema validates the configuration objects used to define data sources in MDV projects. Each datasource represents a collection of data with columns, metadata, and optional features like images, regions, and links to other datasources.

## Schema Files

- `DataSourceSchema.ts` - Main schema definitions
- `datasource-schema.json` - Generated JSON schema for single datasource
- `datasources-array-schema.json` - Generated JSON schema for arrays of datasources

## Core Structure

### DataSource

The main `DataSourceSchema` validates a complete datasource configuration:

```typescript
{
  name: string,                    // Human-readable name
  size: number,                    // Number of rows/items
  columns: DataSourceColumn[],     // Array of column definitions
  columnGroups?: ColumnGroup[],    // Optional column groupings
  links?: LinksConfig,             // Optional links to other datasources
  images?: ImagesConfig,           // Optional thumbnail images
  large_images?: ImagesConfig,     // Optional large images
  regions?: RegionsConfig,         // Optional spatial regions
  interactions?: InteractionsConfig, // Optional interaction data
  offsets?: OffsetsConfig,         // Optional coordinate offsets
  genome_browser?: GenomeBrowserConfig, // Optional genome browser
  tree_diagram?: Record<string, unknown>, // Optional tree diagram
  avivator?: boolean,              // Optional avivator flag
  row_data_loader?: boolean,       // Optional row data loader flag
  binary_data_loader?: boolean,    // Optional binary data loader flag
  deeptools?: Record<string, unknown> // Optional deeptools config
}
```

### DataSourceColumn

Each column in a datasource has the following structure:

```typescript
{
  name: string,                    // Column label for users
  field: string,                   // Unique column identifier
  datatype: DataType,              // Data type (text, integer, double, etc.)
  values?: string[],               // For text columns: possible values
  colors?: string[],               // Color mapping for values
  stringLength?: number,           // For unique/multitext: max string length
  colorLogScale?: boolean,         // For numeric: use log scale for colors
  editable?: boolean,              // Whether column can be edited
  is_url?: boolean,                // Whether column contains URLs
  minMax?: [number, number],       // For numeric: min/max values
  quantiles?: {                    // For numeric: quantile values
    "0.05": [number, number],
    "0.01": [number, number],
    "0.001": [number, number]
  },
  data?: unknown,                  // Actual data (SharedArrayBuffer or array)
  delimiter?: string,              // For multitext: delimiter character
  subgroup?: string,               // For linked columns: subgroup name
  sgindex?: string,                // For linked columns: subgroup index
  sgtype?: string                  // For linked columns: subgroup type
}
```

## Data Types

The schema supports the following data types (defined in `DataTypeSchema`):

- `text` - String with max 256 unique values
- `text16` - String with max 65536 unique values  
- `integer` - Integer (stored as float32)
- `double` - Floating point number (stored as float32)
- `unique` - Any text (variable length)
- `multitext` - Multiple values per item
- `int32` - Large integers (e.g., genomic coordinates)

## Advanced Features

### Column Groups

Logical groupings of columns for better organization:

```typescript
{
  name: string,        // Group name
  columns: string[]    // Array of column field names
}
```

### Images

Configuration for thumbnail and large images associated with data items:

```typescript
{
  [setName: string]: {
    base_url: string,      // Base URL for images
    key_column: string,    // Column containing image keys
    type: string          // Image file type (e.g., "png")
  }
}
```

### Regions

Spatial data configuration for region-based visualizations:

```typescript
{
  position_fields: string[],           // X,Y coordinate columns
  region_field: string,                // Region identifier column
  default_color?: string,              // Default color column
  scale_unit?: string,                 // Scale unit (e.g., "mm")
  scale?: number,                      // Scale factor
  base_url?: string,                   // Base URL for region images
  all_regions: {                       // Region definitions
    [regionName: string]: {
      roi: {                           // Region of interest bounds
        min_x: number, max_x: number,
        min_y: number, max_y: number
      },
      default_image?: string,          // Default background image
      images: {                        // Available images
        [imageName: string]: {
          file: string,                // Image file path
          position: [number, number],  // Image position
          height: number,              // Image height
          width: number,               // Image width
          name: string                 // Image name
        }
      },
      ome_tiff?: string               // OME-TIFF file path
    }
  }
}
```

### Links

Connections between datasources for data sharing and synchronization:

```typescript
{
  [targetDataSource: string]: {
    // Interactions link
    interactions?: {
      interaction_columns: string[],   // Columns defining interactions
      pivot_column: string,            // Shared grouping column
      is_single_region: boolean        // Single vs multi-region
    },
    
    // Rows as columns link
    rows_as_columns?: {
      name_column: string,             // Column identifying rows
      name: string,                    // Human-readable name
      subgroups: {                     // Data subgroups
        [subgroupKey: string]: {
          name: string,                // Subgroup name
          type?: string,               // Data type
          label: string                // User-facing label
        }
      }
    },
    
    // Column sharing link
    columns?: string[],                // Columns to share
    index?: string,                    // Index column for joining
    
    // Data access link
    access_data?: boolean,             // Allow data access
    
    // Color synchronization
    sync_column_colors?: {             // Color sync mappings
      link_to: string,                 // Source column
      col: string                      // Target column
    }[],
    
    // Valueset filtering
    valueset?: {
      source_column: string,           // Source column
      dest_column: string              // Target column
    }
  }
}
```

### Interactions

Configuration for interaction data between objects (e.g., cell interactions):

```typescript
{
  pivot_column: string,                // Grouping column
  interaction_columns: string[],       // Object identifier columns
  is_single_region: boolean,           // Single vs multi-region
  spatial_connectivity_map?: {         // Chart defaults
    link_length: string,
    link_thickness: string,
    link_color: string,
    node_size: string
  },
  interaction_matrix?: {
    groups: string[]
  },
  cell_radial_chart?: {
    link_thickness: string
  }
}
```

### Offsets

Coordinate transformation configuration for spatial alignment:

```typescript
{
  param: string[],                     // Columns to transform (e.g., ["x", "y"])
  groups: string,                      // Grouping column
  background_filter?: string,          // Background filter column
  values?: {                           // Transformations by filter/group
    [filterValue: string]: {
      [groupValue: string]: {
        rotation?: number,             // Rotation in degrees
        offset?: [number, number],     // X,Y translation
        rotation_center?: [number, number] // Rotation center
      }
    }
  }
}
```

### Genome Browser

Configuration for genome browser integration:

```typescript
{
  default_parameters?: Record<string, unknown>,  // Browser parameters
  default_track: {                               // Default feature track
    label: string,                               // Track label
    url: string                                  // Track URL
  },
  default_track_parameters?: Record<string, unknown>, // Track parameters
  location_fields: string[],                     // Chr, start, end columns
  default_tracks?: Record<string, unknown>[],    // Additional tracks
  atac_bam_track?: {                             // ATAC-seq BAM track
    url: string,                                 // BAM file URL
    cluster_read: string                         // Barcode column
  }
}
```

## Usage

### Validation

```typescript
import { validateDataSource, safeValidateDataSource } from './DataSourceSchema';

// Strict validation (throws on error)
const datasource = validateDataSource(config);

// Safe validation (returns null on error)
const datasource = safeValidateDataSource(config);
if (datasource) {
  // Use validated datasource
}
```

### TypeScript Types

```typescript
import type { DataSource, DataSourceColumn } from './DataSourceSchema';

function processDataSource(ds: DataSource) {
  // TypeScript knows the structure
  console.log(ds.name, ds.columns.length);
}
```

## Integration

The datasource schema integrates with:

- **Chart Config Schema**: Uses shared `DataTypeSchema`
- **Python Backend**: JSON schema can be used with Pydantic
- **Frontend Validation**: Runtime validation of datasource configs
- **API Validation**: Server-side validation of datasource requests

## Examples

See the `docs/extradocs/datasource.md` file for comprehensive examples of datasource configurations, including:

- Basic column definitions
- Text and numeric data types
- Image configurations
- Spatial region setups
- Inter-datasource links
- Interaction data
- Genome browser integration 