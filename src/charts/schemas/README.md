# MDV Chart Configuration Schemas

This directory contains Zod schemas for MDV chart configurations that can be used to derive JSON schema for Python.

## Overview

The schemas define the structure of chart configurations in MDV, including:

- **BaseConfigSchema**: Common properties that all charts share
- **Chart-specific schemas**: Extended schemas for each chart type (scatter plot, bar chart, etc.)
- **ChartConfigSchema**: Union of all chart configuration types
- **AddChartDialogConfigSchema**: Configuration for the Add Chart Dialog
- **ViewConfigSchema**: Configuration for MDV views
- **ChartManagerConfigSchema**: Configuration for the Chart Manager

## Usage

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

To generate JSON schema from these Zod schemas for use in Python:

1. **Using Zod v4's built-in method**:
   ```typescript
   import { z } from 'zod/v4';
   import { ChartConfigSchema } from './ChartConfigSchema';
   
   const jsonSchema = z.toJSONSchema(ChartConfigSchema);
   console.log(JSON.stringify(jsonSchema, null, 2));
   ```

2. **Using zod-to-json-schema package**:
   ```bash
   npm install zod-to-json-schema
   ```
   
   ```typescript
   import { zodToJsonSchema } from 'zod-to-json-schema';
   import { ChartConfigSchema } from './ChartConfigSchema';
   
   const jsonSchema = zodToJsonSchema(ChartConfigSchema, {
       name: "ChartConfig",
       description: "Schema for chart configurations in MDV"
   });
   ```

3. **Generating Python Pydantic models**:
   ```bash
   # Install datamodel-code-generator
   pip install datamodel-code-generator
   
   # Generate Python models from JSON schema
   datamodel-codegen --input chartConfig.json --output mdv_models.py
   ```

## Schema Structure

### Base Configuration

All chart configurations extend `BaseConfigSchema` which includes:

- `id`: Unique identifier for the chart
- `size`: Chart dimensions as [width, height]
- `title`: Chart title
- `legend`: Chart description
- `type`: Chart type identifier
- `param`: Array of field specifications (columns)
- `title_color`: Optional title bar color
- Color configuration properties (color_by, color_legend, etc.)

### Field Specifications

Field specifications can be:
- **String**: Direct column name (FieldName)
- **Object**: Serialized column query (e.g., RowsAsColsQuery)

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
- **mypy**: For static type checking

This provides a single source of truth for chart configuration validation across both TypeScript and Python codebases. 