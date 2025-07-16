import * as z from "zod/v4";
import { ChartConfigSchema } from "../src/charts/schemas/ChartConfigSchema.ts";
import { DataSourceSchema, DataSourcesArraySchema } from "../src/charts/schemas/DataSourceSchema.ts";
import * as fs from "node:fs";
import * as path from "node:path";

// Function to generate and save schema
function generateAndSaveSchema(schema: z.ZodTypeAny, filename: string, outputDir = "./schemas") {
    try {
        // Create output directory if it doesn't exist
        if (!fs.existsSync(outputDir)) {
            fs.mkdirSync(outputDir, { recursive: true });
        }

        // Generate JSON schema
        // nb, built-in zod method has less nice internal names for $defs (__schemaN), but I think we prefer not needing to use zod-to-json-schema
        const jsonSchema = z.toJSONSchema(schema, { reused: "ref" });
        
        // Save to file
        const outputPath = path.join(outputDir, filename);
        fs.writeFileSync(outputPath, JSON.stringify(jsonSchema, null, 2));
        
        console.log(`✅ Generated schema: ${outputPath}`);
        return jsonSchema;
    } catch (error) {
        console.error(`❌ Error generating schema for ${filename}:`, error);
        throw error;
    }
}
// Function to extract enums from JSON schema
function extractEnums(jsonSchema: any): Record<string, string[]> {
    const enums: Record<string, string[]> = {};
    
    function traverse(obj: any, path: string[] = []) {
        if (obj && typeof obj === 'object') {
            if (obj.enum && Array.isArray(obj.enum)) {
                const enumName = path[path.length - 1] || 'UnknownEnum';
                enums[enumName] = obj.enum;
            }
            
            for (const [key, value] of Object.entries(obj)) {
                if (key === 'properties' && typeof value === 'object') {
                    for (const [propName, propSchema] of Object.entries(value as any)) {
                        traverse(propSchema, [...path, propName]);
                    }
                } else if (key === 'items' || key === 'anyOf' || key === 'allOf' || key === 'oneOf') {
                    traverse(value, path);
                }
            }
        }
    }
    
    traverse(jsonSchema);
    return enums;
}

// Function to generate Python enum
function generateEnum(enumName: string, values: string[]): string {
    let enumCode = `class ${enumName}(str, Enum):\n`;
    for (const value of values) {
        const constName = value.toUpperCase().replace(/[^A-Z0-9]/g, '_');
        enumCode += `    ${constName} = "${value}"\n`;
    }
    return enumCode;
}

// Function to generate Python class from JSON schema
function generateClassFromSchema(jsonSchema: any, className: string): string {
    let classCode = `class ${className}(BaseModel):\n`;
    
    if (jsonSchema.description) {
        classCode += `    """${jsonSchema.description}"""\n`;
    }
    
    if (jsonSchema.properties) {
        for (const [propName, propSchema] of Object.entries(jsonSchema.properties)) {
            const fieldType = getPythonType(propSchema as any, propName);
            const fieldDescription = (propSchema as any).description || '';
            const isRequired = jsonSchema.required?.includes(propName) || false;
            
            if (isRequired) {
                classCode += `    ${propName}: ${fieldType} = Field(description="${fieldDescription}")\n`;
            } else {
                classCode += `    ${propName}: Optional[${fieldType}] = Field(None, description="${fieldDescription}")\n`;
            }
        }
    }
    
    return classCode;
}

// Function to get Python type from JSON schema
function getPythonType(schema: any, fieldName: string): string {
    if (!schema) return 'Any';
    if (schema.type === 'string') {
        if (schema.enum) {
            return 'str';
        }
        return 'str';
    } else if (schema.type === 'number' || schema.type === 'integer') {
        return 'float';
    } else if (schema.type === 'boolean') {
        return 'bool';
    } else if (schema.type === 'array') {
        const itemType = getPythonType(schema.items, `${fieldName}_item`);
        return `List[${itemType}]`;
    } else if (schema.type === 'object') {
        return 'Dict[str, Any]';
    } else if (schema.anyOf || schema.oneOf) {
        const types = (schema.anyOf || schema.oneOf).map((s: any) => getPythonType(s, fieldName));
        return `Union[${types.join(', ')}]`;
    } else if (schema.$ref) {
        // Handle references - for now return Any
        return 'Any';
    }
    // Defensive: if schema.type is missing or unknown
    return 'Any';
}

// Function to generate Python Pydantic models file
async function generatePythonModels(): Promise<string> {
    // Import the schemas we need to process
    const { ChartConfigSchema } = await import("../src/charts/schemas/ChartConfigSchema.ts");
    const { DataSourceSchema, DataSourcesArraySchema } = await import("../src/charts/schemas/DataSourceSchema.ts");
    
    // Convert schemas to JSON Schema format using Zod's toJSONSchema method
    const chartConfigJsonSchema = z.toJSONSchema(ChartConfigSchema);
    const datasourceJsonSchema = z.toJSONSchema(DataSourceSchema);
    const datasourcesArrayJsonSchema = z.toJSONSchema(DataSourcesArraySchema);

    // Extract enums from schemas
    const chartEnums = extractEnums(chartConfigJsonSchema);
    const datasourceEnums = extractEnums(datasourceJsonSchema);
    const allEnums = { ...chartEnums, ...datasourceEnums };

    // Generate Python code
    let pythonCode = `"""
Auto-generated Pydantic models from Zod schemas.
This file is generated automatically - do not edit manually.
Generated from TypeScript Zod schemas to ensure consistency.
"""

from typing import List, Optional, Union, Dict, Any, Tuple
from pydantic import BaseModel, Field
from enum import Enum

`;

    // Generate enums
    for (const [enumName, values] of Object.entries(allEnums)) {
        pythonCode += generateEnum(enumName, values);
        pythonCode += "\n\n";
    }

    // Generate individual chart classes from the union schema
    const chartConfigAny = chartConfigJsonSchema as any;
    if (chartConfigAny.anyOf) {
        const chartClasses: string[] = [];
        const chartTypeNames: string[] = [];
        
        for (const schema of chartConfigAny.anyOf) {
            // Only generate a class if this schema has properties (is an object type)
            if ((schema as any).properties && typeof (schema as any).properties === 'object' && (schema as any).properties.type?.const) {
                const chartType = (schema as any).properties.type.const;
                // biome-ignore lint/style/useTemplate: wrapping complex expression would be harder to read.
                const className = chartType.split('_').map((word: string) => 
                    word.charAt(0).toUpperCase() + word.slice(1)
                ).join('') + 'Config';
                
                chartClasses.push(generateClassFromSchema(schema, className));
                chartTypeNames.push(className);
            }
            // else: skip schemas that are not object types (no properties)
        }
        
        // Add all chart classes
        pythonCode += chartClasses.join('\n\n');
        pythonCode += '\n\n';
        
        // Add the union type
        pythonCode += `# Union of all chart configuration types\nChartConfig = Union[\n    ${chartTypeNames.join(',\n    ')}\n]\n\n`;
    }

    // Generate datasource classes
    pythonCode += generateClassFromSchema(datasourceJsonSchema, "DataSource");
    pythonCode += "\n\n";

    // Generate datasources array class
    pythonCode += generateClassFromSchema(datasourcesArrayJsonSchema, "DataSourcesArray");
    pythonCode += "\n\n";

    // Add structured output classes for LLM integration
    pythonCode += `# Structured output for data analysis
class DataAnalysisResult(BaseModel):
    """Result of structured data analysis"""
    selected_columns: List[str] = Field(description="List of selected column names from the dataframes")
    chart_types: List[str] = Field(description="List of suitable chart types for the selected columns")
    reasoning: str = Field(description="Brief explanation of the column selection")
    data_summary: Dict[str, Any] = Field(description="Summary statistics of the selected data")

# Structured output for chart generation
class ChartGenerationResult(BaseModel):
    """Result of structured chart generation"""
    chart_config: ChartConfig = Field(description="Generated chart configuration")
    python_code: str = Field(description="Generated Python code to create the chart")
    view_name: str = Field(description="Name for the generated view")
    dependencies: List[str] = Field(description="List of required Python packages")

# Complete structured output for a query
class StructuredQueryResult(BaseModel):
    """Complete structured result for a user query"""
    question: str = Field(description="Original user question")
    data_analysis: DataAnalysisResult = Field(description="Analysis of the data")
    chart_generation: ChartGenerationResult = Field(description="Generated chart")
    execution_result: Optional[Dict[str, Any]] = Field(None, description="Result of code execution")
    error: Optional[str] = Field(None, description="Error message if any")
    success: bool = Field(description="Whether the query was successful")

# Utility functions
def get_chart_type_from_config(config: ChartConfig) -> str:
    """Extract chart type from configuration"""
    if hasattr(config, 'type'):
        return config.type
    return "unknown"

def validate_chart_config(config_dict: Dict[str, Any]) -> ChartConfig:
    """Validate a chart configuration dictionary"""
    # This would need to be implemented based on the specific chart type
    # For now, we'll use the base config
    # ChartConfig is a Union, so we try each possible config class (there may be a better way)
    last_error = None
    for config_type in ChartConfig.__args__:
        try:
            return config_type(**config_dict)
        except Exception as e:
            last_error = e
    raise ValueError(f"Could not validate chart config: {last_error}")

def create_chart_config(
    chart_type: str,
    title: str,
    legend: str,
    params: List[str],
    size: Tuple[float, float] = (800, 600),
    **kwargs
) -> ChartConfig:
    """Create a chart configuration based on type and parameters"""
    base_config = {
        "id": f"chart_{hash(title)}",
        "title": title,
        "legend": legend,
        "param": params,
        "size": size,
        **kwargs
    }
    
    # ChartConfig is a Union, so we try each possible config class (there may be a better way)
    last_error = None
    for config_type in ChartConfig.__args__:
        try:
            return config_type(**base_config)
        except Exception as e:
            last_error = e
    raise ValueError(f"Could not create chart config: {last_error}")
`;

    return pythonCode;
}

// Function to save Python models
async function savePythonModels(outputDir: string) {
    try {
        // Create output directory if it doesn't exist
        if (!fs.existsSync(outputDir)) {
            fs.mkdirSync(outputDir, { recursive: true });
        }

        const pythonCode = await generatePythonModels();
        const outputPath = path.join(outputDir, "structured_schemas.py");
        fs.writeFileSync(outputPath, pythonCode);
        
        console.log(`✅ Generated Python models: ${outputPath}`);
    } catch (error) {
        console.error(`❌ Error generating Python models: ${error}`);
        throw error;
    }
}

// Main execution
async function main() {
    console.log("🚀 Starting schema generation...");
    
    try {
        // Generate main chart config schema
        const chartConfigSchema = generateAndSaveSchema(
            ChartConfigSchema, 
            "chart-config-schema.json"
        );
        
        console.log("📊 Chart config schema generated successfully!");

        // Generate datasource schema
        const datasourceSchema = generateAndSaveSchema(
            DataSourceSchema,
            "datasource-schema.json"
        );

        console.log("📊 Datasource schema generated successfully!");
        console.log("📋 Schema includes validation for datasource configuration");

        // Generate datasources array schema
        const datasourcesArraySchema = generateAndSaveSchema(
            DataSourcesArraySchema,
            "datasources-array-schema.json"
        );

        console.log("📊 Datasources array schema generated successfully!");
        console.log("📋 Schema includes validation for arrays of datasources");
        
        // Generate Python Pydantic models
        await savePythonModels("./python/mdvtools/llm");
        
        console.log("🐍 Python Pydantic models generated successfully!");
        console.log("📋 Models include all chart types and structured output classes");
        
        console.log("\n✨ Schema generation complete!");
        console.log("🔒 Security note: Structured output approach enables secure server-side execution");
        console.log("   - No arbitrary code execution via PythonAstREPLTool");
        console.log("   - Type-safe schema validation");
        console.log("   - Controlled output generation");
        
    } catch (error) {
        console.error("💥 Schema generation failed:", error);
        process.exit(1);
    }
}

// Run if this file is executed directly
if (import.meta.url === `file://${process.argv[1]}`) {
    // we could have a watch mode...
    main();
}

export { generateAndSaveSchema, savePythonModels };