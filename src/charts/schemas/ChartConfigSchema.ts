import { z } from "zod/v4";
import { registerChartConfigSchema, getAllChartConfigSchemas } from "./ChartConfigRegistry";

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
    subgroupName: z.string().optional().describe(
        "Canonical rows-as-columns subgroup key (see link `subgroups`). New `RowsAsColsQuery.toJSON` output always includes this; omit = first subgroup. Zod: keep optional until legacy configs are migrated, then make required.",
    ),
    maxItems: z.number().int().positive().describe("Maximum number of items to retrieve from the linked datasource"),
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
    gsposition: z.tuple([z.number().int().nonnegative(), z.number().int().nonnegative()]).optional().describe("GridStack position as [x, y] coordinates"),
}).describe("GridStack layout configuration for dashboard positioning");

// Allow `title` to be provided as either a string or a single-element string[]
const TitleSchema = z.preprocess((val) => {
    if (Array.isArray(val) && val.length === 1 && typeof val[0] === "string") {
        return val[0];
    }
    return val;
}, z.string()).describe("Display title for the chart");

// Base configuration that all charts extend
export const BaseConfigSchema = GridStackConfigSchema.extend({
    id: z.string().describe("Unique identifier for the chart instance"),
    size: z.tuple([z.number(), z.number()]).describe("Chart dimensions as [width, height] in pixels"),
    //nb: seems like tables in default view have e.g. `title: [\"cells\"]`
    title: TitleSchema,
    legend: z.string().optional().describe("Legend text describing the chart's purpose"),
    type: z.string().describe("Chart type identifier (e.g., 'scatter_plot', 'bar_chart')"),
    param: FieldSpecsSchema.describe("Array of field specifications defining the data columns used by this chart"),
    title_color: z.string().optional().describe("CSS color value for the chart title"),
    version: z.string().optional().describe("Schema version for this chart configuration (used for future migration support)"),
}).describe("Base configuration shared by all chart types, including layout, styling, and data mapping");

// todo - this should be a union of all the color config options, not applied to all charts
export const ChartColorConfigSchema = z.object({
    color_by: FieldSpecSchema.optional().describe("Field or column configuration used to determine color mapping"),
    color_legend: z
        .object({
            // Keep optional for backwards compatibility with legacy/partially-serialized
            // configs where `color_legend` exists but `display` may be missing/misspelled.
            display: z.boolean().optional().describe("Whether the color legend is visible"),
            pos: z.tuple([z.number(), z.number()]).optional().describe("Legend position in pixels [left, top]"),
        })
        .optional()
        .describe("Color legend display and position"),
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

const CategorySelectionSchema = z.union([
    z.string(),
    z.array(z.string()),
]).describe("One or more selected categorical values");

const CategoryFilterSchema = z.object({
    column: z.string().describe("Categorical column used for filtering"),
    category: CategorySelectionSchema.describe("Selected categorical values"),
}).describe("A categorical filter applied to a scatter plot");

const FieldLegendConfigSchema = z.object({
    display: z.boolean().describe("Whether to display the density field legend"),
}).describe("Configuration for density field legend visibility");

const ContourVisualConfigSchema = z.object({
    contour_fill: z.boolean().optional().describe("Whether to render filled contour areas"),
    contour_fillThreshold: z.number().optional().describe("Threshold controlling how much of the contour is filled"),
    contour_bandwidth: z.number().positive().optional().describe("Kernel density bandwidth used for contour generation"),
    contour_intensity: z.number().min(0).max(1).optional().describe("Opacity/intensity applied to contour fills"),
    contour_opacity: z.number().min(0).max(1).optional().describe("Opacity applied to contour lines"),
}).describe("Visual contour configuration shared by density-style charts");

const DensityFieldsConfigSchema = z.object({
    densityFields: FieldSpecsSchema.optional().describe("Numeric fields rendered as density overlays"),
}).describe("Density field selection shared by contour-based charts");

const LegacyContourCategorySelectionConfigSchema = z.object({
    contourParameter: FieldSpecSchema.optional().describe("Categorical field used to choose contour categories"),
    category1: CategorySelectionSchema.optional().describe("First contour category selection"),
    category2: CategorySelectionSchema.optional().describe("Second contour category selection"),
}).describe("Legacy contour category selection shared by deck scatter density charts");

const DensityFieldLegendConfigSchema = z.object({
    field_legend: FieldLegendConfigSchema.optional().describe("Legend settings for density field overlays"),
}).describe("Legend configuration for density field overlays");

const ContourScatterConfigSchema = z.object({
    ...ContourVisualConfigSchema.shape,
    ...LegacyContourCategorySelectionConfigSchema.shape,
    ...DensityFieldsConfigSchema.shape,
    ...DensityFieldLegendConfigSchema.shape,
}).describe("Contour-specific configuration shared by deck scatter charts and spatial viewers");

const ScatterAxisSchema = z.object({
    size: z.number().describe("Reserved size for this axis in pixels"),
    tickfont: z.number().describe("Axis tick font size"),
    rotate_labels: z.boolean().describe("Whether to rotate axis tick labels"),
}).describe("Configuration for a single scatter plot axis");

const ScatterAxis2DSchema = z.object({
    x: ScatterAxisSchema.describe("X axis settings"),
    y: ScatterAxisSchema.describe("Y axis settings"),
}).describe("Configuration for 2D scatter plot axes");

const OrthographicViewStateSchema = z.object({
    target: z.tuple([z.number(), z.number(), z.number()]).describe("Orthographic view target in world coordinates"),
    zoom: z.number().describe("Current zoom level"),
    minZoom: z.number().optional().describe("Minimum allowed zoom"),
    maxZoom: z.number().optional().describe("Maximum allowed zoom"),
}).passthrough().describe("Serializable orthographic deck.gl view state");

const OrbitViewStateSchema = OrthographicViewStateSchema.extend({
    rotationOrbit: z.number().optional().describe("Orbit rotation in degrees"),
    rotationX: z.number().optional().describe("X-axis rotation in degrees"),
}).passthrough().describe("Serializable orbit deck.gl view state");

const DeckScatterSharedConfigSchema = BaseConfigSchema.extend({
    ...ChartColorConfigSchema.shape,
    ...TooltipConfigSchema.shape,
    ...ContourScatterConfigSchema.shape,
    course_radius: z.number().positive().optional().describe("Coarse multiplier applied to the scatter point radius"),
    radius: z.number().positive().optional().describe("Scatter point radius"),
    opacity: z.number().min(0).max(1).optional().describe("Scatter point opacity"),
    category_filters: z.array(CategoryFilterSchema).optional().describe("Categorical filters applied before rendering"),
    on_filter: z.enum(["hide", "grey"]).optional().describe("How filtered-out points should be rendered"),
    zoom_on_filter: z.boolean().optional().describe("Whether filtering should auto-fit the current view"),
    point_shape: z.enum(["circle", "square", "gaussian"]).optional().describe("Scatter point rendering shape"),
    selectionFeatureCollection: z.unknown().optional().describe("Serialized selection overlay geometry"),
    region: z.string().optional().describe("Associated region identifier when rendering spatial scatter plots"),
    hideMissing: z.boolean().optional().describe("Whether to hide points with missing numeric color values"),
}).describe("Shared configuration for deck.gl scatter charts");

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

export const CategoryHeatmapConfigSchema = BaseConfigSchema.extend({
    type: z.literal("category_heatmap").describe("Category heatmap chart type"),
    x_display_categories: z.array(z.string()).optional().describe("Optional x-axis categories to display in the heatmap"),
    y_display_categories: z.array(z.string()).optional().describe("Optional y-axis categories to display in the heatmap"),
}).describe("Configuration for category heatmaps showing categorical co-occurrence counts");

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
    type: z.enum(["text_box_chart", "text_box"]).describe("Text box chart type"),
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

export const DeckScatterPlotConfigSchema = DeckScatterSharedConfigSchema.extend({
    type: z.literal("DeckScatter").describe("Deck.gl 2D scatter plot chart type"),
    dimension: z.literal("2d").optional().describe("Deck scatter plots default to 2D"),
    axis: ScatterAxis2DSchema.optional().describe("Axis configuration for the 2D scatter view"),
    viewState: OrthographicViewStateSchema.optional().describe("Current deck.gl orthographic view state"),
}).describe("Configuration for the deck.gl 2D scatter plot");

export const DeckContourScatterConfigSchema = DeckScatterSharedConfigSchema.extend({
    type: z.literal("DeckContourScatter").describe("Deck.gl contour-enabled scatter plot chart type"),
    dimension: z.literal("2d").optional().describe("Contour scatter plots render in 2D"),
    axis: ScatterAxis2DSchema.optional().describe("Axis configuration for the contour scatter view"),
    viewState: OrthographicViewStateSchema.optional().describe("Current deck.gl orthographic view state"),
}).describe("Configuration for the deck.gl scatter plot with contour overlays");

export const DeckDensityConfigSchema = DeckScatterSharedConfigSchema.extend({
    type: z.literal("DeckDensity").describe("Legacy alias for the deck.gl contour scatter chart type"),
    dimension: z.literal("2d").optional().describe("Contour scatter plots render in 2D"),
    axis: ScatterAxis2DSchema.optional().describe("Axis configuration for the contour scatter view"),
    viewState: OrthographicViewStateSchema.optional().describe("Current deck.gl orthographic view state"),
}).describe("Configuration for legacy deck.gl density scatter charts");

export const DeckScatter3DConfigSchema = DeckScatterSharedConfigSchema.extend({
    type: z.literal("DeckScatter3D").describe("Deck.gl 3D scatter plot chart type"),
    dimension: z.literal("3d").optional().describe("3D deck scatter plots render in orbit view"),
    viewState: OrbitViewStateSchema.optional().describe("Current deck.gl orbit view state"),
}).describe("Configuration for the deck.gl 3D scatter plot");

export const DeckSplatterConfigSchema = BaseConfigSchema.extend({
    type: z.literal("DeckSplatter").describe("Splatter plot chart type"),
    category: FieldSpecSchema.optional().describe("Categorical column used for the row layout"),
    ...DensityFieldsConfigSchema.shape,
    ...ContourVisualConfigSchema.shape,
}).describe("Configuration for splatter plots showing density fields across a category-by-field grid");

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

export const GeneInfoSchema = BaseConfigSchema.extend({
    type: z.literal("geneNetwork").describe("Gene information from genenetwork.nl API"),
    autoScroll: z.boolean().optional().describe("Automatically scroll to highlighted items")
}).describe("Configuration for chart allowing for selection of genes with information retrieved from genenetwork.nl API");

// Register all chart-specific schemas with the registry
// This establishes the pattern for future co-location of schemas with chart implementations
registerChartConfigSchema("scatter_plot", ScatterPlotConfigSchema, { version: "1" });
registerChartConfigSchema("bar_chart", BarChartConfigSchema, { version: "1" });
registerChartConfigSchema("histogram", HistogramConfigSchema, { version: "1" });
registerChartConfigSchema("heatmap", HeatmapConfigSchema, { version: "1" });
registerChartConfigSchema("category_heatmap", CategoryHeatmapConfigSchema, { version: "1" });
registerChartConfigSchema("dot_plot", DotPlotConfigSchema, { version: "1" });
registerChartConfigSchema("box_plot", BoxPlotConfigSchema, { version: "1" });
registerChartConfigSchema("violin_plot", ViolinPlotConfigSchema, { version: "1" });
registerChartConfigSchema("pie_chart", PieChartConfigSchema, { version: "1" });
registerChartConfigSchema("table", TableConfigSchema, { version: "1" });
registerChartConfigSchema("text_box_chart", TextBoxConfigSchema, { version: "1" });
registerChartConfigSchema("text_box", TextBoxConfigSchema, { version: "1" });
registerChartConfigSchema("wordcloud", WordcloudConfigSchema, { version: "1" });
registerChartConfigSchema("sankey", SankeyConfigSchema, { version: "1" });
registerChartConfigSchema("multi_line_chart", MultiLineChartConfigSchema, { version: "1" });
registerChartConfigSchema("multi_box_plot", MultiBoxPlotConfigSchema, { version: "1" });
registerChartConfigSchema("abundance_box_plot", AbundanceBoxPlotConfigSchema, { version: "1" });
registerChartConfigSchema("density_scatter", DensityScatterConfigSchema, { version: "1" });
registerChartConfigSchema("DeckScatter", DeckScatterPlotConfigSchema, { version: "1" });
registerChartConfigSchema("DeckContourScatter", DeckContourScatterConfigSchema, { version: "1" });
registerChartConfigSchema("DeckDensity", DeckDensityConfigSchema, { version: "1" });
registerChartConfigSchema("DeckScatter3D", DeckScatter3DConfigSchema, { version: "1" });
registerChartConfigSchema("DeckSplatter", DeckSplatterConfigSchema, { version: "1" });
registerChartConfigSchema("row_chart", RowChartConfigSchema, { version: "1" });
registerChartConfigSchema("stacked_row_chart", StackedRowChartConfigSchema, { version: "1" });
registerChartConfigSchema("row_summary_box", RowSummaryBoxConfigSchema, { version: "1" });
registerChartConfigSchema("selection_dialog", SelectionDialogConfigSchema, { version: "1" });
registerChartConfigSchema("geneNetwork", GeneInfoSchema, { version: "1" });

// Build union of all chart configuration types from the registry
// This ensures the union is derived from registered schemas, establishing the pattern
// for future co-location where schemas register themselves as side-effects
const registeredSchemas = getAllChartConfigSchemas();
const allSchemas = registeredSchemas.length > 0 
    ? [...registeredSchemas, BaseConfigSchema]
    : [BaseConfigSchema];
export const ChartConfigSchema = z.union(allSchemas as [z.ZodTypeAny, z.ZodTypeAny, ...z.ZodTypeAny[]])
    .describe("Union of all possible chart configuration types, allowing for type-safe chart creation");

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
    onlyView: ViewConfigSchema.optional().describe("Single view configuration when only one view is needed"),
    show_gallery_on_open: z.boolean().optional().describe("Whether to show the gallery view by default when the project is opened"),
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
