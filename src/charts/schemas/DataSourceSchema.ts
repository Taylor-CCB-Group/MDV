import { z } from "zod";
// @ts-expect-error - currently the way we use from node won't resolve without this
import { DataTypeSchema } from "./ChartConfigSchema.ts";

// Column schema for datasource columns
export const DataSourceColumnSchema = z.object({
    name: z.string().describe("Human-readable name for the column displayed to users"),
    field: z.string().describe("Unique identifier for the column used internally"),
    datatype: DataTypeSchema.describe("Data type defining how the column values are stored and processed"),
    values: z.array(z.string()).optional().describe("For text/multitext columns: array of possible values, with raw data containing indices into this array"),
    colors: z.array(z.string()).optional().describe("Array of hex color codes mapped to values or interpolated for numeric data"),
    stringLength: z.number().int().positive().optional().describe("For unique columns: length of longest value; for multitext: max values per item"),
    colorLogScale: z.boolean().optional().describe("Whether to use logarithmic scaling for color mapping of numeric columns"),
    editable: z.boolean().optional().describe("Whether users can edit values in this column"),
    is_url: z.boolean().optional().describe("Whether this column contains URL links (for unique and text columns)"),
    minMax: z.tuple([z.number(), z.number()]).optional().describe("For numeric columns: [minimum, maximum] values for optimization"),
    quantiles: z.object({
        "0.05": z.tuple([z.number(), z.number()]).describe("5th and 95th percentiles"),
        "0.01": z.tuple([z.number(), z.number()]).describe("1st and 99th percentiles"),
        "0.001": z.tuple([z.number(), z.number()]).describe("0.1st and 99.9th percentiles")
    }).optional().describe("Pre-calculated quantiles for efficient color scaling and outlier handling"),
    data: z.any().optional().describe("Column data as array or SharedArrayBuffer - usually loaded on demand"),
    delimiter: z.string().optional().describe("Delimiter character for multitext columns (default: comma)"),
    subgroup: z.string().optional().describe("Subgroup identifier for rows-as-columns queries"),
    sgindex: z.string().optional().describe("Unique identifier for rows-as-columns data"),
    sgtype: z.string().optional().describe("Data type for rows-as-columns (e.g., 'sparse')")
}).describe("Configuration for a single column in a datasource, defining its properties, data type, and behavior");

// Column group schema
export const ColumnGroupSchema = z.object({
    name: z.string().describe("Human-readable name for the column group"),
    columns: z.array(z.string()).describe("Array of column field names belonging to this group")
}).describe("Logical grouping of related columns for organizational purposes");

// Image set schema
export const ImageSetSchema = z.object({
    base_url: z.string().describe("Base URL where images are stored - image key and type will be appended"),
    key_column: z.string().describe("Column name containing the keys used to construct image URLs"),
    type: z.string().describe("Image file extension (e.g., 'png', 'jpg')")
}).describe("Configuration for a set of images associated with datasource rows");

// Images configuration schema
export const ImagesConfigSchema = z.record(z.string(), ImageSetSchema).describe("Multiple named image sets for the datasource");

// ROI schema for regions
export const ROISchema = z.object({
    min_x: z.number().describe("Minimum x coordinate of the region"),
    max_x: z.number().describe("Maximum x coordinate of the region"),
    min_y: z.number().describe("Minimum y coordinate of the region"),
    max_y: z.number().describe("Maximum y coordinate of the region")
}).describe("Bounding box defining the region of interest (ROI)");

// Image metadata schema
export const ImageMetadataSchema = z.object({
    file: z.string().describe("Filename of the image"),
    position: z.tuple([z.number(), z.number()]).describe("[x, y] offset position of the image within the region"),
    height: z.number().describe("Height of the image in pixels"),
    width: z.number().describe("Width of the image in pixels"),
    name: z.string().describe("Human-readable name for the image")
}).describe("Metadata for a background image within a spatial region");

// Region schema
export const RegionSchema = z.object({
    roi: ROISchema.describe("Bounding box defining the spatial extent of this region"),
    default_image: z.string().nullable().optional().describe("Name of the default background image to display (null = no image)"),
    images: z.record(z.string(), ImageMetadataSchema).describe("Dictionary of available background images for this region"),
    ome_tiff: z.string().optional().describe("Filename of the OME-TIFF image associated with this region")
}).describe("Configuration for a single spatial region with its boundaries and associated images");

// Regions configuration schema
export const RegionsConfigSchema = z.object({
    position_fields: z.array(z.string()).describe("Column names containing x,y coordinates for spatial data"),
    region_field: z.string().describe("Column name identifying which region each data point belongs to"),
    default_color: z.string().optional().describe("Default column to use for color mapping in spatial visualizations"),
    scale_unit: z.string().optional().describe("Unit of measurement for spatial scale (e.g., 'mm', 'μm')"),
    scale: z.number().optional().describe("Scale factor converting pixel coordinates to real-world units"),
    base_url: z.string().optional().describe("Base URL where region images are stored"),
    all_regions: z.record(z.string(), RegionSchema).describe("Dictionary of all regions, keyed by region identifier")
}).describe("Configuration for spatial data with multiple regions and associated background images");

// Interaction chart defaults schema
export const InteractionChartDefaultsSchema = z.object({
    spatial_connectivity_map: z.object({
        link_length: z.string().describe("Column name for link length in spatial connectivity maps"),
        link_thickness: z.string().describe("Column name for link thickness in spatial connectivity maps"),
        link_color: z.string().describe("Column name for link color in spatial connectivity maps"),
        node_size: z.string().describe("Column name for node size in spatial connectivity maps")
    }).optional().describe("Default column mappings for spatial connectivity map charts"),
    interaction_matrix: z.object({
        groups: z.array(z.string()).describe("Column names for grouping interactions in matrix charts")
    }).optional().describe("Default grouping configuration for interaction matrix charts"),
    cell_radial_chart: z.object({
        link_thickness: z.string().describe("Column name for link thickness in radial charts")
    }).optional().describe("Default column mapping for cell radial charts")
}).describe("Default chart configurations for interaction visualization types");

// Interactions configuration schema
export const InteractionsConfigSchema = z.object({
    pivot_column: z.string().describe("Column defining the scope of interactions (e.g., region, condition, sample)"),
    interaction_columns: z.array(z.string()).describe("Two columns specifying the interacting entities (e.g., cell types)"),
    is_single_region: z.boolean().describe("Whether interactions are within single regions (true) or across conditions (false)")
}).and(InteractionChartDefaultsSchema).describe("Configuration for interaction data between entities across spatial or experimental contexts");

// Offset values schema
export const OffsetValuesSchema = z.object({
    rotation: z.number().optional().describe("Rotation angle in degrees (positive = clockwise)"),
    offset: z.tuple([z.number(), z.number()]).optional().describe("[x, y] translation offset in coordinate units"),
    rotation_center: z.tuple([z.number(), z.number()]).optional().describe("[x, y] center point for rotation (relative to data bounds)")
}).describe("Spatial transformation parameters for aligning data groups");

// Offsets configuration schema
export const OffsetsConfigSchema = z.object({
    param: z.array(z.string()).describe("Column names that can be transformed (typically ['x', 'y'])"),
    groups: z.string().describe("Column name identifying groups that can have different transformations"),
    background_filter: z.string().optional().describe("Column name for filtering transformations by spatial context (e.g., 'ROI')"),
    values: z.record(z.string(), z.record(z.string(), OffsetValuesSchema)).optional().describe("Nested dictionary: {filter_value: {group_value: transformation}}")
}).describe("Configuration for spatial transformations to align different data groups or conditions");

// Genome browser track schema
export const GenomeBrowserTrackSchema = z.object({
    label: z.string().describe("Human-readable name for the track"),
    url: z.string().describe("URL to the track data file (gzipped BED format with tabix index)")
}).describe("Configuration for a single track in the genome browser");

// ATAC BAM track schema
export const ATACBAMTrackSchema = z.object({
    url: z.string().describe("URL to the BAM file containing ATAC-seq reads"),
    cluster_read: z.string().describe("Column name in linked datasource for clustering reads by barcode")
}).describe("Configuration for ATAC-seq BAM track with barcode-based clustering");

// Genome browser configuration schema
export const GenomeBrowserConfigSchema = z.object({
    default_parameters: z.record(z.string(), z.unknown()).optional().describe("Default parameters for the genome browser widget"),
    default_track: GenomeBrowserTrackSchema.describe("Primary track showing datasource features"),
    default_track_parameters: z.record(z.string(), z.unknown()).optional().describe("Parameters for the default feature track"),
    location_fields: z.array(z.string()).describe("Column names for [chromosome, start, end] genomic coordinates"),
    default_tracks: z.array(z.record(z.string(), z.unknown())).optional().describe("Additional tracks to display by default"),
    atac_bam_track: ATACBAMTrackSchema.optional().describe("ATAC-seq BAM track configuration if applicable")
}).describe("Configuration for genome browser integration with genomic coordinate data");

// Rows as columns subgroup schema
export const RowsAsColumnsSubgroupSchema = z.object({
    name: z.string().describe("Internal name for the subgroup"),
    type: z.string().optional().describe("Data type for the subgroup (e.g., 'sparse' for sparse matrices)"),
    label: z.string().describe("Human-readable label for the subgroup")
}).describe("Configuration for a subgroup within rows-as-columns data (e.g., gene expression scores)");

// Rows as columns link schema
export const RowsAsColumnsLinkSchema = z.object({
    name_column: z.string().describe("Column in linked datasource containing unique identifiers (e.g., gene names)"),
    name: z.string().describe("Human-readable name for the rows-as-columns data (e.g., 'Gene Scores')"),
    subgroups: z.record(z.string(), RowsAsColumnsSubgroupSchema).describe("Dictionary of available data subgroups")
}).describe("Configuration for treating rows from another datasource as columns in the current datasource");

// Sync column colors link schema
export const SyncColumnColorsLinkSchema = z.array(z.object({
    link_to: z.string().describe("Column name in the linked datasource to sync colors from"),
    col: z.string().describe("Column name in the current datasource to sync colors to")
})).describe("Array of column color synchronization mappings between datasources");

// Valueset link schema
export const ValuesetLinkSchema = z.object({
    source_column: z.string().describe("Column in source datasource containing values to filter by"),
    dest_column: z.string().describe("Column in destination datasource to apply filtering to")
}).describe("Configuration for filtering one datasource based on unique values from another");

// Link schema - union of all possible link types
export const LinkSchema = z.object({
    interactions: InteractionsConfigSchema.optional().describe("Configuration for interaction data between entities"),
    rows_as_columns: RowsAsColumnsLinkSchema.optional().describe("Configuration for treating rows as columns"),
    columns: z.array(z.string()).optional().describe("Array of column names to link from another datasource"),
    index: z.string().optional().describe("Column name used as index for linking datasources"),
    access_data: z.boolean().optional().describe("Whether this datasource can access data from the linked datasource"),
    sync_column_colors: SyncColumnColorsLinkSchema.optional().describe("Configuration for synchronizing column colors"),
    valueset: ValuesetLinkSchema.optional().describe("Configuration for valueset-based filtering")
}).describe("Configuration for different types of relationships between datasources");

// Links configuration schema
export const LinksConfigSchema = z.record(z.string(), LinkSchema).describe("Dictionary of datasource links, keyed by linked datasource name");

// Main datasource schema
export const DataSourceSchema = z.object({
    name: z.string().describe("Human-readable name for the datasource"),
    size: z.number().int().positive().describe("Number of rows in the datasource"),
    columns: z.array(DataSourceColumnSchema).describe("Array of column configurations defining the data structure"),
    columnGroups: z.array(ColumnGroupSchema).optional().describe("Logical groupings of related columns"),
    links: LinksConfigSchema.optional().describe("Relationships with other datasources"),
    images: ImagesConfigSchema.optional().describe("Thumbnail images associated with datasource rows"),
    large_images: ImagesConfigSchema.optional().describe("High-resolution images for detailed viewing"),
    regions: RegionsConfigSchema.optional().describe("Spatial regions configuration for spatial data"),
    interactions: InteractionsConfigSchema.optional().describe("Interaction data configuration"),
    offsets: OffsetsConfigSchema.optional().describe("Spatial transformation configuration"),
    genome_browser: GenomeBrowserConfigSchema.optional().describe("Genome browser integration configuration"),
    tree_diagram: z.record(z.string(), z.unknown()).optional().describe("Tree diagram visualization configuration"),
    avivator: z.boolean().optional().describe("Whether Avivator image viewer integration is enabled"),
    row_data_loader: z.boolean().optional().describe("Whether to use row-based data loading"),
    binary_data_loader: z.boolean().optional().describe("Whether to use binary data loading for performance"),
    deeptools: z.record(z.string(), z.unknown()).optional().describe("DeepTools integration configuration")
}).describe("Complete configuration for a datasource including data structure, relationships, and visualization options");

// Array of datasources schema
export const DataSourcesArraySchema = z.array(DataSourceSchema).describe("Array of multiple datasource configurations");

// Export types
export type DataSourceColumn = z.infer<typeof DataSourceColumnSchema>;
export type ColumnGroup = z.infer<typeof ColumnGroupSchema>;
export type ImageSet = z.infer<typeof ImageSetSchema>;
export type ImagesConfig = z.infer<typeof ImagesConfigSchema>;
export type ROI = z.infer<typeof ROISchema>;
export type ImageMetadata = z.infer<typeof ImageMetadataSchema>;
export type Region = z.infer<typeof RegionSchema>;
export type RegionsConfig = z.infer<typeof RegionsConfigSchema>;
export type InteractionChartDefaults = z.infer<typeof InteractionChartDefaultsSchema>;
export type InteractionsConfig = z.infer<typeof InteractionsConfigSchema>;
export type OffsetValues = z.infer<typeof OffsetValuesSchema>;
export type OffsetsConfig = z.infer<typeof OffsetsConfigSchema>;
export type GenomeBrowserTrack = z.infer<typeof GenomeBrowserTrackSchema>;
export type ATACBAMTrack = z.infer<typeof ATACBAMTrackSchema>;
export type GenomeBrowserConfig = z.infer<typeof GenomeBrowserConfigSchema>;
export type RowsAsColumnsSubgroup = z.infer<typeof RowsAsColumnsSubgroupSchema>;
export type RowsAsColumnsLink = z.infer<typeof RowsAsColumnsLinkSchema>;
export type SyncColumnColorsLink = z.infer<typeof SyncColumnColorsLinkSchema>;
export type ValuesetLink = z.infer<typeof ValuesetLinkSchema>;
export type Link = z.infer<typeof LinkSchema>;
export type LinksConfig = z.infer<typeof LinksConfigSchema>;
export type DataSource = z.infer<typeof DataSourceSchema>;
export type DataSourcesArray = z.infer<typeof DataSourcesArraySchema>;

// Utility function to validate a datasource
export function validateDataSource(config: unknown): DataSource {
    return DataSourceSchema.parse(config);
}

// Utility function to safely validate a datasource (returns null if invalid)
export function safeValidateDataSource(config: unknown): DataSource | null {
    const result = DataSourceSchema.safeParse(config);
    return result.success ? result.data : null;
}

// Utility function to validate an array of datasources
export function validateDataSourcesArray(config: unknown): DataSourcesArray {
    return DataSourcesArraySchema.parse(config);
}

// Utility function to safely validate an array of datasources (returns null if invalid)
export function safeValidateDataSourcesArray(config: unknown): DataSourcesArray | null {
    const result = DataSourcesArraySchema.safeParse(config);
    return result.success ? result.data : null;
} 