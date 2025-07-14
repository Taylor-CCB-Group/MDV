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
        const jsonSchema = z.toJSONSchema(schema);
        
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
        console.log("📋 Schema includes validation for all chart types:");
        
        // Log the chart types that are included in the schema
        if (chartConfigSchema.anyOf) {
            chartConfigSchema.anyOf.forEach((schema: any, index: number) => {
                if (schema.properties?.type?.const) {
                    console.log(`   - ${schema.properties.type.const}`);
                }
            });
        }

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
        
        console.log("\n✨ Schema generation complete!");
        
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

export { generateAndSaveSchema };