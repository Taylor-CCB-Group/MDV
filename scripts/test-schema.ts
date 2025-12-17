import * as z from "zod/v4";
import { 
    BaseConfigSchema, 
    ChartColorConfigSchema, 
    FieldSpecSchema,
    TooltipConfigSchema,
    ScatterPlotConfigSchema,
    ChartConfigSchema 
} from "../src/charts/schemas/ChartConfigSchema.ts";

// Test individual schemas to isolate the problem
function testSchema(name: string, schema: z.ZodTypeAny) {
    try {
        console.log(`Testing ${name}...`);
        const jsonSchema = z.toJSONSchema(schema);
        console.log(`✅ ${name} works`);
        return true;
    } catch (error) {
        console.error(`❌ ${name} failed:`, error);
        return false;
    }
}

console.log("🔍 Testing individual schemas...");

// Test schemas in dependency order
const results = [
    testSchema("FieldSpecSchema", FieldSpecSchema),
    testSchema("BaseConfigSchema", BaseConfigSchema),
    testSchema("ChartColorConfigSchema", ChartColorConfigSchema),
    testSchema("TooltipConfigSchema", TooltipConfigSchema),
    testSchema("ScatterPlotConfigSchema", ScatterPlotConfigSchema),
    testSchema("ChartConfigSchema", ChartConfigSchema),
];

console.log("\n📊 Results:");
results.forEach((result, index) => {
    const schemas = ["FieldSpecSchema", "BaseConfigSchema", "ChartColorConfigSchema", "TooltipConfigSchema", "ScatterPlotConfigSchema", "ChartConfigSchema"];
    console.log(`${result ? "✅" : "❌"} ${schemas[index]}`);
}); 