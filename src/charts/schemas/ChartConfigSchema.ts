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
]).describe("Data types for chart columns: integer/double for numbers, text/text16 for categorical data with value limits, unique for unlimited text, multitext for multiple values per item, int32 for large integers");

// Parameter types for chart configuration
export const ParamTypeSchema = z.enum([
    "text",
    "number", 
    "multitext",
    "text16",
    "_multi_column:number",
    "_multi_column:all"
]).describe("Parameter types for chart configuration: basic types plus special multi-column selectors for complex data relationships");

// Serialized form of column queries
export const RowsAsColsQuerySerializedSchema = z.object({
    type: z.literal("RowsAsColsQuery").describe("Identifies this as a rows-as-columns query"),
    linkedDsName: z.string().describe("Name of the linked datasource that contains the row data"),
    // TODO - should have subgroup name
    maxItems: z.number().int().positive().describe("Maximum number of items to retrieve from the linked datasource")
}).describe("Configuration for queries that treat rows from another datasource as columns in the current chart");

// Field specification can be a string (FieldName) or a serialized query
export const FieldSpecSchema = z.union([
    z.string().describe("Direct field name reference"), // FieldName
    RowsAsColsQuerySerializedSchema
]).describe("A field specification can be either a direct column name or a complex query to another datasource");

// Field specifications (array of field specs)
export const FieldSpecsSchema = z.array(FieldSpecSchema).describe("Array of field specifications defining which data columns to use in the chart");

export const GridStackConfigSchema = z.object({
    gssize: z.tuple([z.number().int().positive(), z.number().int().positive()]).optional().describe("GridStack size as [width, height] in grid units"),
    gsposition: z.tuple([z.number().int().positive(), z.number().int().positive()]).optional().describe("GridStack position as [x, y] coordinates"),
}).describe("GridStack layout configuration for dashboard positioning");

// Base configuration that all charts extend
export const BaseConfigSchema = GridStackConfigSchema.extend({
    id: z.string().describe("Unique identifier for the chart instance"),
    size: z.tuple([z.number(), z.number()]).describe("Chart dimensions as [width, height] in pixels"),
    title: z.string().describe("Display title for the chart"),
    legend: z.string().describe("Legend text describing the chart's purpose"),
    type: z.string().describe("Chart type identifier (e.g., 'scatter_plot', 'bar_chart')"),
    param: FieldSpecsSchema.describe("Array of field specifications defining the data columns used by this chart"),
    title_color: z.string().optional().describe("CSS color value for the chart title"),
}).describe("Base configuration shared by all chart types, including layout, styling, and data mapping");

// todo - this should be a union of all the color config options, not applied to all charts
export const ChartColorConfigSchema = z.object({
    color_by: FieldSpecSchema.optional().describe("Field or column configuration used to determine color mapping"),
    color_legend: z.record(z.string(), z.unknown()).optional().describe("Custom color legend configuration"),
    log_color_scale: z.boolean().optional().describe("Whether to use logarithmic scaling for color values"),
    trim_color_scale: z.enum(["0.05", "0.01", "0.001", "none"]).optional().describe("Quantile trimming for color scale to handle outliers"),
    color_overlay: z.number().optional().describe("Opacity value for color overlay effects"),
    fallbackOnZero: z.boolean().optional().describe("Whether to use fallback colors when values are zero")
}).describe("Configuration for chart color mapping and scaling");


// Tooltip configuration
export const TooltipConfigSchema = z.object({
    tooltip: z.object({
        show: z.boolean().describe("Whether to display tooltips on hover"),
        column: z.union([FieldSpecSchema, FieldSpecsSchema]).optional().describe("Field(s) to display in tooltip content")
    })
}).describe("Configuration for chart tooltip behavior and content");

// Chart-specific configuration schemas
export const ScatterPlotConfigSchema = BaseConfigSchema.extend({
    ...ChartColorConfigSchema.shape,
    type: z.literal("scatter_plot").describe("Scatter plot chart type"),
    opacity: z.number().min(0).max(1).optional().describe("Opacity of scatter plot points (0-1)"),
    radius: z.number().positive().optional().describe("Radius of scatter plot points (units should be better defined)"),
    // Additional scatter plot specific properties
}).describe("Configuration for scatter plot charts showing relationships between two variables");

export const BarChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("bar_chart").describe("Bar chart type"),
    // Additional bar chart specific properties
}).describe("Configuration for bar charts displaying categorical data");

export const HistogramConfigSchema = BaseConfigSchema.extend({
    type: z.literal("histogram").describe("Histogram chart type"),
    bins: z.number().int().positive().optional().describe("Number of bins for histogram data grouping"),
    // Additional histogram specific properties
}).describe("Configuration for histogram charts showing data distribution");

export const HeatmapConfigSchema = BaseConfigSchema.extend({
    type: z.literal("heatmap").describe("Heatmap chart type"),
    // Additional heatmap specific properties
}).describe("Configuration for heatmap charts displaying matrix data with color intensity");

export const DotPlotConfigSchema = BaseConfigSchema.extend({
    type: z.literal("dot_plot").describe("Dot plot chart type"),
    // Additional dot plot specific properties
}).describe("Configuration for dot plots showing individual data points");

export const BoxPlotConfigSchema = BaseConfigSchema.extend({
    ...ChartColorConfigSchema.shape,
    type: z.literal("box_plot").describe("Box plot chart type"),
    // Additional box plot specific properties
}).describe("Configuration for box plots showing data distribution statistics");

export const ViolinPlotConfigSchema = BaseConfigSchema.extend({
    ...ChartColorConfigSchema.shape,
    type: z.literal("violin_plot").describe("Violin plot chart type"),
    // Additional violin plot specific properties
}).describe("Configuration for violin plots showing data density distribution");

export const PieChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("pie_chart").describe("Pie chart type"),
    // Additional pie chart specific properties
}).describe("Configuration for pie charts showing proportional data");

export const TableConfigSchema = BaseConfigSchema.extend({
    type: z.literal("table").describe("Table chart type"),
    // Additional table specific properties
}).describe("Configuration for data tables displaying tabular information");

export const TextBoxConfigSchema = BaseConfigSchema.extend({
    type: z.literal("text_box").describe("Text box chart type"),
    text: z.string().describe("Text content to display in the text box"),
    // Additional text box specific properties
}).describe("Configuration for text boxes displaying static or dynamic text content");

export const WordcloudConfigSchema = BaseConfigSchema.extend({
    type: z.literal("wordcloud").describe("Wordcloud chart type"),
    // Additional wordcloud specific properties
}).describe("Configuration for wordcloud charts showing text frequency");

export const SankeyConfigSchema = BaseConfigSchema.extend({
    type: z.literal("sankey").describe("Sankey diagram type"),
    // Additional sankey specific properties
}).describe("Configuration for Sankey diagrams showing flow between categories");

export const MultiLineChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("multi_line_chart").describe("Multi-line chart type"),
    stacked: z.boolean().optional().describe("Whether to stack multiple lines on top of each other"),
    band_width: z.number().positive().optional().describe("Width of confidence bands around lines"),
    intervals: z.number().int().positive().optional().describe("Number of intervals for data aggregation"),
    scaletrim: z.boolean().optional().describe("Whether to trim the scale to focus on data range"),
    // Additional multi-line chart specific properties
}).describe("Configuration for multi-line charts showing multiple time series or categories");

export const MultiBoxPlotConfigSchema = BaseConfigSchema.extend({
    type: z.literal("multi_box_plot").describe("Multi-box plot chart type"),
    // Additional multi-box plot specific properties
}).describe("Configuration for multi-box plots comparing distributions across categories");

export const AbundanceBoxPlotConfigSchema = BaseConfigSchema.extend({
    type: z.literal("abundance_box_plot").describe("Abundance box plot chart type"),
    // Additional abundance box plot specific properties
}).describe("Configuration for abundance box plots showing distribution of abundance data");

export const DensityScatterConfigSchema = BaseConfigSchema.extend({
    ...ChartColorConfigSchema.shape,
    type: z.literal("density_scatter").describe("Density scatter plot chart type"),
    // Additional density scatter specific properties
}).describe("Configuration for density scatter plots showing point density with color intensity");

export const RowChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("row_chart").describe("Row chart type"),
    // Additional row chart specific properties
}).describe("Configuration for row charts displaying data in horizontal bars");

export const StackedRowChartConfigSchema = BaseConfigSchema.extend({
    type: z.literal("stacked_row_chart").describe("Stacked row chart type"),
    // Additional stacked row chart specific properties
}).describe("Configuration for stacked row charts with multiple data series per row");

export const RowSummaryBoxConfigSchema = BaseConfigSchema.extend({
    type: z.literal("row_summary_box").describe("Row summary box chart type"),
    // Additional row summary box specific properties
}).describe("Configuration for row summary boxes displaying aggregated statistics");

export const SelectionDialogConfigSchema = BaseConfigSchema.extend({
    type: z.literal("selection_dialog").describe("Selection dialog chart type"),
    // Additional selection dialog specific properties
}).describe("Configuration for selection dialogs allowing user interaction");

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
]).describe("Union of all possible chart configuration types, allowing for type-safe chart creation");

// Schema for view configuration
export const ViewConfigSchema = z.object({
    dataSources: z.record(z.string(), z.object({
        layout: z.enum(["absolute", "gridstack"]).optional().describe("Layout system for positioning charts: absolute for fixed positioning, gridstack for responsive grid"),
        panelWidth: z.number().min(0).max(100).optional().describe("Panel width as percentage of viewport (0-100)")
    })).describe("Configuration for each datasource in the view"),
    initialCharts: z.record(z.string(), z.array(ChartConfigSchema)).describe("Initial charts to display for each datasource")
}).describe("Configuration for a complete view containing multiple datasources and their charts");

// Schema for chart manager configuration
export const ChartManagerConfigSchema = z.object({
    initialCharts: z.array(ChartConfigSchema).optional().describe("Initial charts to create when the manager starts"),
    all_views: z.array(z.string()).optional().describe("List of all available view names"),
    current_view: z.string().optional().describe("Currently active view name"),
    permission: z.string().optional().describe("Permission level for chart operations"),
    gridstack: z.boolean().optional().describe("Whether to use GridStack for chart layout"),
    chat_enabled: z.boolean().optional().describe("Whether chat functionality is enabled"),
    mdv_api_root: z.string().optional().describe("Root URL for MDV API endpoints"),
    onlyView: ViewConfigSchema.optional().describe("Single view configuration when only one view is needed")
}).describe("Top-level configuration for the chart manager controlling the entire dashboard");

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
