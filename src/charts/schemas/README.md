# MDV Chart Configuration Schemas

This directory contains Zod schemas for MDV chart configurations that can be used to derive JSON schema for Python and drive UI generation.

## Overview

The schemas define the structure of chart configurations in MDV, including:

- **BaseConfigSchema**: Common properties that all charts share (layout, basic styling)
- **ChartColorConfigSchema**: Color-related properties that some charts extend
- **Chart-specific schemas**: Extended schemas for each chart type (scatter plot, bar chart, etc.)
- **ChartConfigSchema**: Union of all chart configuration types
- **ViewConfigSchema**: Configuration for MDV views
- **ChartManagerConfigSchema**: Configuration for the Chart Manager

## Current State

This is in a experimental form to explore and guide how we shape the design in future.


### Schema Structure

- **BaseConfigSchema**: Contains layout, identification, and basic styling properties
- **ChartColorConfigSchema**: Separate schema for color-related properties (color_by, color_legend, etc.)
- **FieldSpecSchema**: Supports both direct column references and complex queries
- **Chart-specific schemas**: Each chart type extends BaseConfigSchema and optionally ChartColorConfigSchema

We also generate schema for `datasources`.

### Field Specifications

Field specifications can be:
- **String**: Direct column name (FieldName)
- **Object**: Serialized column query (e.g., RowsAsColsQuery)

## Usage

The `npm run build-schemas` script will output both JSON schema definitions & Python source-code.

### TypeScript/JavaScript

```typescript
import { ChartConfigSchema, validateChartConfig } from './ChartConfigSchema';

// Validate a chart configuration
const config = {
    id: "chart1",
    size: [800, 600],
    title: "My Scatter Plot",
    legend: "A scatter plot",
    type: "scatter_plot",
    param: ["x_column", "y_column"],
    opacity: 0.8,
    radius: 5.0
};

try {
    const validatedConfig = validateChartConfig(config);
    console.log("Valid configuration:", validatedConfig);
} catch (error) {
    console.error("Invalid configuration:", error);
}
```

### Generating JSON Schema for Python

Currently, the build process uses simple TypeScript code to generate Python strings. A future plan is to investigate richer code generation tools like `datamodel-code-generator` for more robust Python model generation.

**Current approach**: Simple TS string generation means that python doesn't need to be manually updated to reflect changes in model
**Future plan**: Use `datamodel-code-generator` for richer Python code generation

```bash
# Future: Install datamodel-code-generator using poetry
poetry add datamodel-code-generator

# Future: Generate Python models from JSON schema
datamodel-codegen --input chartConfig.json --output mdv_models.py
```

**Benefits of datamodel-code-generator**:
- Automatic Pydantic model generation
- Better type annotations and validation
- Support for complex schema features
- Integration with Python tooling (pyright, IDE support)
- More maintainable than manual string generation

## Future Plans

### 1. Schema-Driven UI Generation

**Goal**: Generate UI controls (AddChart dialog, settings dialog) directly from schemas.

**Benefits**:
- Reduce boilerplate code
- Ensure consistency between validation and UI
- Enable automatic UI updates when schemas change
- Support LLM generation of both configs and UIs

**Implementation**:
- Add metadata to schemas for UI hints (e.g., `.meta({ ui: { type: "slider", range: [0, 1] } })`)
- Create schema-to-UI generators for common control types
- Maintain custom settings for complex interactions

### 2. Explicit Field Properties

**Goal**: Move from generic `param` array to explicit, named properties.

**Current**: `param: ["category", "field1", "field2"]`
**Future**: `category: "category", fields: ["field1", "field2"]`

**Benefits**:
- Self-documenting configs
- Better type safety
- Easier LLM understanding and generation
- Eliminates complex array mapping logic

### 3. Field Type Constraints

**Goal**: Express column type constraints in schemas.

**Implementation**:
- Extend `FieldSpecSchema` with metadata for constraints
- Create specialized schemas (e.g., `CategoricalFieldSpecSchema`, `NumericFieldSpecSchema`)
- Use metadata to drive UI column filtering and validation

### 4. Graph-Based Data Transformations

**Goal**: Support complex data transformations through visual graph editing.

**Features**:
- React Flow UI for building transformation graphs
- Nodes for operations like normalization, scaling, custom expressions
- Schema-driven node definitions
- Support for both simple values and complex transformations

**Example Use Cases**:
- Scatter plot radius: fixed value, column values, or normalized column values
- Color mapping: direct column or transformed values
- Multi-step data processing pipelines

### 5. Migration Strategy

**Approach**: Parallel schema registry with migration functions.

**Components**:
- `BaseChart.schemas` (parallel to `BaseChart.types`)
- Migration functions for each chart type
- Version information in schemas
- Tests for migration correctness

**Benefits**:
- Backward compatibility with existing configs
- Gradual migration path
- Validation of migration functions

## Schema Structure

### Base Configuration

All chart configurations extend `BaseConfigSchema` which includes:

- `id`: Unique identifier for the chart
- `size`: Chart dimensions as [width, height]
- `title`: Chart title
- `legend`: Chart description
- `type`: Chart type identifier
- `param`: Array of field specifications (columns) - *to be replaced with explicit properties*
- `title_color`: Optional title bar color
- GridStack layout properties

### Color Configuration

Charts that support color mapping extend `ChartColorConfigSchema`:

- `color_by`: Field or column configuration for color mapping
- `color_legend`: Custom color legend configuration
- `log_color_scale`: Logarithmic scaling for color values
- `trim_color_scale`: Quantile trimming for outliers
- `color_overlay`: Opacity for color overlay effects
- `fallbackOnZero`: Fallback colors for zero values

### Chart Types

The schemas support the following chart types:
- Scatter Plot
- Bar Chart
- Histogram
- Heatmap
- Dot Plot
- Box Plot
- Violin Plot
- Pie Chart
- Table
- Text Box
- Wordcloud
- Sankey Diagram
- Multi-line Chart
- Multi-box Plot
- Abundance Box Plot
- Density Scatter
- Row Chart
- Stacked Row Chart
- Row Summary Box
- Selection Dialog

## Example JSON Schema Output

When converted to JSON schema, the `ChartConfigSchema` produces a schema like:

```json
{
  "$schema": "http://json-schema.org/draft-07/schema#",
  "title": "ChartConfig",
  "description": "Schema for chart configurations in MDV",
  "oneOf": [
    {
      "type": "object",
      "properties": {
        "id": { "type": "string" },
        "size": { 
          "type": "array", 
          "items": [{ "type": "number" }, { "type": "number" }],
          "minItems": 2,
          "maxItems": 2
        },
        "title": { "type": "string" },
        "legend": { "type": "string" },
        "type": { "const": "scatter_plot" },
        "param": {
          "type": "array",
          "items": {
            "oneOf": [
              { "type": "string" },
              {
                "type": "object",
                "properties": {
                  "type": { "const": "RowsAsColsQuery" },
                  "linkedDsName": { "type": "string" },
                  "maxItems": { "type": "integer", "minimum": 1 }
                },
                "required": ["type", "linkedDsName", "maxItems"]
              }
            ]
          }
        },
        "opacity": { "type": "number", "minimum": 0, "maximum": 1 },
        "radius": { "type": "number", "minimum": 0 }
      },
      "required": ["id", "size", "title", "legend", "type", "param"]
    }
    // ... other chart types
  ]
}
```

## Integration with Python

The generated JSON schema can be used with:

- **Pydantic**: For data validation and serialization
- **datamodel-code-generator**: For automatic model generation
- **jsonschema**: For runtime validation
- **pyright**: For static type checking

This provides a single source of truth for chart configuration validation across both TypeScript and Python codebases. 