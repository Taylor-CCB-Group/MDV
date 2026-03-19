import type { z } from "zod/v4";

/**
 * Options for registering a chart config schema
 */
export interface ChartConfigSchemaOptions {
    /**
     * Version string for this schema (e.g., "1", "2.0")
     * Used for future versioning and migration support
     */
    version?: string;
    /**
     * Migration function to transform older config versions to this version
     * Not yet implemented, but reserved for future use
     */
    migrate?: (oldConfig: any) => any;
}

/**
 * Internal storage for registered schemas
 * Maps chart type -> version -> schema
 */
type SchemaRegistry = Map<string, Map<string, { schema: z.ZodTypeAny; options: ChartConfigSchemaOptions }>>;

const registry: SchemaRegistry = new Map();

/**
 * Register a chart configuration schema
 * 
 * @param type - The chart type identifier (e.g., "scatter_plot", "multi_line_chart")
 * @param schema - The Zod schema for this chart type
 * @param options - Optional version and migration function
 * 
 * @example
 * ```ts
 * registerChartConfigSchema("scatter_plot", ScatterPlotConfigSchema, { version: "1" });
 * ```
 */
export function registerChartConfigSchema(
    type: string,
    schema: z.ZodTypeAny,
    options: ChartConfigSchemaOptions = {}
): void {
    const version = options.version || "1";
    
    if (!registry.has(type)) {
        registry.set(type, new Map());
    }
    
    const typeRegistry = registry.get(type);
    // we don't have a type-guard for it, but we could assert here that typeRegistry is there.
    typeRegistry?.set(version, { schema, options });
}

/**
 * Get a chart configuration schema by type and optional version
 * 
 * @param type - The chart type identifier
 * @param version - Optional version string. If not provided, returns the latest registered version
 * @returns The Zod schema for the chart type, or undefined if not found
 */
export function getChartConfigSchema(
    type: string,
    version?: string
): z.ZodTypeAny | undefined {
    const typeRegistry = registry.get(type);
    if (!typeRegistry) {
        return undefined;
    }
    
    if (version) {
        // Return specific version if requested
        const entry = typeRegistry.get(version);
        return entry?.schema;
    }
    
    // Return the latest version (highest version string, or "1" if no version specified)
    // For now, we'll just return the first entry or the one with version "1"
    // In the future, we could implement proper version comparison
    let latestEntry: { schema: z.ZodTypeAny; options: ChartConfigSchemaOptions } | undefined;
    let latestVersion: string | undefined;
    
    for (const [v, entry] of typeRegistry.entries()) {
        if (!latestVersion || v > latestVersion) {
            latestVersion = v;
            latestEntry = entry;
        }
    }
    
    return latestEntry?.schema;
}

/**
 * Get all registered chart config schemas
 * Useful for building union types
 * 
 * @returns Array of all registered Zod schemas
 */
export function getAllChartConfigSchemas(): z.ZodTypeAny[] {
    const schemas: z.ZodTypeAny[] = [];
    
    for (const typeRegistry of registry.values()) {
        for (const entry of typeRegistry.values()) {
            schemas.push(entry.schema);
        }
    }
    
    return schemas;
}

/**
 * Check if a chart type is registered
 * 
 * @param type - The chart type identifier
 * @returns True if the chart type has been registered
 */
export function isChartTypeRegistered(type: string): boolean {
    return registry.has(type);
}

/**
 * Get all registered chart types
 * 
 * @returns Array of all registered chart type identifiers
 */
export function getRegisteredChartTypes(): string[] {
    return Array.from(registry.keys());
}

