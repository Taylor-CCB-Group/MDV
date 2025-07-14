import { z } from "zod/v4";

// Basic data types that can be stored in columns
export const DataTypeSchema = z.enum([
    "integer",
    "double", 
    "text",
    "text16",
    "unique",
    "multitext",
    "int32"
]);

// Parameter types for chart configuration
export const ParamTypeSchema = z.enum([
    "text",
    "number", 
    "multitext",
    "text16",
    "_multi_column:number",
    "_multi_column:all"
]);

// Serialized form of column queries
export const RowsAsColsQuerySerializedSchema = z.object({
    type: z.literal("RowsAsColsQuery"),
    linkedDsName: z.string(),
    // TODO - should have subgroup name
    maxItems: z.number().int().positive()
});

// Field specification can be a string (FieldName) or a serialized query
export const FieldSpecSchema = z.union([
    z.string(), // FieldName
    RowsAsColsQuerySerializedSchema
]);

// Field specifications (array of field specs)
export const FieldSpecsSchema = z.array(FieldSpecSchema);

export const GridStackConfigSchema = z.object({
    gssize: z.tuple([z.number().int().positive(), z.number().int().positive()]).optional(),
    gsposition: z.tuple([z.number().int().positive(), z.number().int().positive()]).optional(),
});

// Base configuration that all charts extend
export const BaseConfigSchema = GridStackConfigSchema.extend({
    id: z.string(),
    size: z.tuple([z.number(), z.number()]),
    title: z.string(),
    legend: z.string(),
    type: z.string(),
    param: FieldSpecsSchema,
    title_color: z.string().optional(),
    // Color configuration properties
    color_by: z.union([
        z.string(), // FieldName
        z.object({
            column: z.object({
                field: z.string(),
                name: z.string(),
                datatype: DataTypeSchema,
                // Other column properties can be added as needed
            }).loose()
        })
    ]).optional(),
    color_legend: z.any().optional(),
    log_color_scale: z.boolean().optional(),
    trim_color_scale: z.union([
        z.enum(["0.05", "0.01", "0.001"]),
        z.literal("none")
    ]).optional(),
    color_overlay: z.number().optional(),
    fallbackOnZero: z.boolean().optional()
});

// Tooltip configuration
export const TooltipConfigSchema = z.object({
    tooltip: z.object({
        show: z.boolean(),
        column: z.union([FieldSpecSchema, FieldSpecsSchema]).optional()
    })
});

// Chart-specific configuration schemas
export const ScatterPlotConfigSchema = BaseConfigSchema.extend({
    type: z.literal("scatter_plot"),
    opacity: z.number().min(0).max(1).optional(),
    radius: z.number().positive().optional(),
    // Additional scatter plot specific properties
}).loose();

export const BarChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("bar_chart"),
    // Additional bar chart specific properties
}).loose();

export const HistogramConfigSchema = BaseConfigSchema.extend({
    type: z.literal("histogram"),
    bins: z.number().int().positive().optional(),
    // Additional histogram specific properties
}).loose();

export const HeatmapConfigSchema = BaseConfigSchema.extend({
    type: z.literal("heatmap"),
    // Additional heatmap specific properties
}).loose();

export const DotPlotConfigSchema = BaseConfigSchema.extend({
    type: z.literal("dot_plot"),
    // Additional dot plot specific properties
}).loose();

export const BoxPlotConfigSchema = BaseConfigSchema.extend({
    type: z.literal("box_plot"),
    // Additional box plot specific properties
}).loose();

export const ViolinPlotConfigSchema = BaseConfigSchema.extend({
    type: z.literal("violin_plot"),
    // Additional violin plot specific properties
}).loose();

export const PieChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("pie_chart"),
    // Additional pie chart specific properties
}).loose();

export const TableConfigSchema = BaseConfigSchema.extend({
    type: z.literal("table"),
    // Additional table specific properties
}).loose();

export const TextBoxConfigSchema = BaseConfigSchema.extend({
    type: z.literal("text_box"),
    text: z.string(),
    // Additional text box specific properties
}).loose();

export const WordcloudConfigSchema = BaseConfigSchema.extend({
    type: z.literal("wordcloud"),
    // Additional wordcloud specific properties
}).loose();

export const SankeyConfigSchema = BaseConfigSchema.extend({
    type: z.literal("sankey"),
    // Additional sankey specific properties
}).loose();

export const MultiLineChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("multi_line_chart"),
    stacked: z.boolean().optional(),
    band_width: z.number().positive().optional(),
    intervals: z.number().int().positive().optional(),
    scaletrim: z.boolean().optional(),
    // Additional multi-line chart specific properties
}).loose();

export const MultiBoxPlotConfigSchema = BaseConfigSchema.extend({
    type: z.literal("multi_box_plot"),
    // Additional multi-box plot specific properties
}).loose();

export const AbundanceBoxPlotConfigSchema = BaseConfigSchema.extend({
    type: z.literal("abundance_box_plot"),
    // Additional abundance box plot specific properties
}).loose();

export const DensityScatterConfigSchema = BaseConfigSchema.extend({
    type: z.literal("density_scatter"),
    // Additional density scatter specific properties
}).loose();

export const RowChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("row_chart"),
    // Additional row chart specific properties
}).loose();

export const StackedRowChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("stacked_row_chart"),
    // Additional stacked row chart specific properties
}).loose();

export const RowSummaryBoxConfigSchema = BaseConfigSchema.extend({
    type: z.literal("row_summary_box"),
    // Additional row summary box specific properties
}).loose();

export const SelectionDialogConfigSchema = BaseConfigSchema.extend({
    type: z.literal("selection_dialog"),
    // Additional selection dialog specific properties
}).loose();

// Union of all chart configuration types
export const ChartConfigSchema = z.union([
    ScatterPlotConfigSchema,
    BarChartConfigSchema,
    HistogramConfigSchema,
    HeatmapConfigSchema,
    DotPlotConfigSchema,
    BoxPlotConfigSchema,
    ViolinPlotConfigSchema,
    PieChartConfigSchema,
    TableConfigSchema,
    TextBoxConfigSchema,
    WordcloudConfigSchema,
    SankeyConfigSchema,
    MultiLineChartConfigSchema,
    MultiBoxPlotConfigSchema,
    AbundanceBoxPlotConfigSchema,
    DensityScatterConfigSchema,
    RowChartConfigSchema,
    StackedRowChartConfigSchema,
    RowSummaryBoxConfigSchema,
    SelectionDialogConfigSchema,
    // Fallback for any other chart types
    BaseConfigSchema
]);

// Schema for the AddChartDialog configuration (intermediate form)
// export const AddChartDialogConfigSchema = z.object({
//     title: z.string(),
//     legend: z.string(),
//     param: z.optional(z.array(z.union([z.string(), z.array(z.string())]))),
//     type: z.string(),
//     extra: z.record(z.unknown()),
//     _updated: z.date()
// });

// Schema for view configuration
export const ViewConfigSchema = z.object({
    dataSources: z.record(z.string(), z.object({
        layout: z.enum(["absolute", "gridstack"]).optional(),
        panelWidth: z.number().min(0).max(100).optional()
    })),
    initialCharts: z.record(z.string(), z.array(ChartConfigSchema))
});

// Schema for chart manager configuration
export const ChartManagerConfigSchema = z.object({
    initialCharts: z.array(ChartConfigSchema).optional(),
    all_views: z.array(z.string()).optional(),
    current_view: z.string().optional(),
    permission: z.string().optional(),
    gridstack: z.boolean().optional(),
    chat_enabled: z.boolean().optional(),
    mdv_api_root: z.string().optional(),
    onlyView: ViewConfigSchema.optional()
});

// Export types
export type DataType = z.infer<typeof DataTypeSchema>;
export type ParamType = z.infer<typeof ParamTypeSchema>;
export type FieldSpec = z.infer<typeof FieldSpecSchema>;
export type FieldSpecs = z.infer<typeof FieldSpecsSchema>;
export type TooltipConfig = z.infer<typeof TooltipConfigSchema>;
export type BaseConfig = z.infer<typeof BaseConfigSchema>;
export type ChartConfig = z.infer<typeof ChartConfigSchema>;
// export type AddChartDialogConfig = z.infer<typeof AddChartDialogConfigSchema>;
export type ViewConfig = z.infer<typeof ViewConfigSchema>;
export type ChartManagerConfig = z.infer<typeof ChartManagerConfigSchema>;

// Utility function to validate a chart config
export function validateChartConfig(config: unknown): ChartConfig {
    return ChartConfigSchema.parse(config);
}

// Utility function to safely validate a chart config (returns null if invalid)
export function safeValidateChartConfig(config: unknown): ChartConfig | null {
    const result = ChartConfigSchema.safeParse(config);
    return result.success ? result.data : null;
}
