import { z } from "zod/v4";
// @ts-expect-error - currently the way we use from node won't resolve without this
import { DataTypeSchema } from "./ChartConfigSchema.ts";

// Column schema for datasource columns
export const DataSourceColumnSchema = z.object({
    name: z.string(),
    field: z.string(),
    datatype: DataTypeSchema,
    values: z.array(z.string()).optional(),
    colors: z.array(z.string()).optional(),
    stringLength: z.number().int().positive().optional(),
    colorLogScale: z.boolean().optional(),
    editable: z.boolean().optional(),
    is_url: z.boolean().optional(),
    minMax: z.tuple([z.number(), z.number()]).optional(),
    quantiles: z.object({
        "0.05": z.tuple([z.number(), z.number()]),
        "0.01": z.tuple([z.number(), z.number()]),
        "0.001": z.tuple([z.number(), z.number()])
    }).optional(),
    data: z.any().optional(), // SharedArrayBuffer or array
    delimiter: z.string().optional(), // for multitext columns
    subgroup: z.string().optional(),
    sgindex: z.string().optional(),
    sgtype: z.string().optional()
});

// Column group schema
export const ColumnGroupSchema = z.object({
    name: z.string(),
    columns: z.array(z.string())
});

// Image set schema
export const ImageSetSchema = z.object({
    base_url: z.string(),
    key_column: z.string(),
    type: z.string()
});

// Images configuration schema
export const ImagesConfigSchema = z.record(z.string(), ImageSetSchema);

// ROI schema for regions
export const ROISchema = z.object({
    min_x: z.number(),
    max_x: z.number(),
    min_y: z.number(),
    max_y: z.number()
});

// Image metadata schema
export const ImageMetadataSchema = z.object({
    file: z.string(),
    position: z.tuple([z.number(), z.number()]),
    height: z.number(),
    width: z.number(),
    name: z.string()
});

// Region schema
export const RegionSchema = z.object({
    roi: ROISchema,
    default_image: z.string().nullable().optional(),
    images: z.record(z.string(), ImageMetadataSchema),
    ome_tiff: z.string().optional()
});

// Regions configuration schema
export const RegionsConfigSchema = z.object({
    position_fields: z.array(z.string()),
    region_field: z.string(),
    default_color: z.string().optional(),
    scale_unit: z.string().optional(),
    scale: z.number().optional(),
    base_url: z.string().optional(),
    all_regions: z.record(z.string(), RegionSchema)
});

// Interaction chart defaults schema
export const InteractionChartDefaultsSchema = z.object({
    spatial_connectivity_map: z.object({
        link_length: z.string(),
        link_thickness: z.string(),
        link_color: z.string(),
        node_size: z.string()
    }).optional(),
    interaction_matrix: z.object({
        groups: z.array(z.string())
    }).optional(),
    cell_radial_chart: z.object({
        link_thickness: z.string()
    }).optional()
});

// Interactions configuration schema
export const InteractionsConfigSchema = z.object({
    pivot_column: z.string(),
    interaction_columns: z.array(z.string()),
    is_single_region: z.boolean()
}).and(InteractionChartDefaultsSchema);

// Offset values schema
export const OffsetValuesSchema = z.object({
    rotation: z.number().optional(),
    offset: z.tuple([z.number(), z.number()]).optional(),
    rotation_center: z.tuple([z.number(), z.number()]).optional()
});

// Offsets configuration schema
export const OffsetsConfigSchema = z.object({
    param: z.array(z.string()),
    groups: z.string(),
    background_filter: z.string().optional(),
    values: z.record(z.string(), z.record(z.string(), OffsetValuesSchema)).optional()
});

// Genome browser track schema
export const GenomeBrowserTrackSchema = z.object({
    label: z.string(),
    url: z.string()
});

// ATAC BAM track schema
export const ATACBAMTrackSchema = z.object({
    url: z.string(),
    cluster_read: z.string()
});

// Genome browser configuration schema
export const GenomeBrowserConfigSchema = z.object({
    default_parameters: z.record(z.string(), z.unknown()).optional(),
    default_track: GenomeBrowserTrackSchema,
    default_track_parameters: z.record(z.string(), z.unknown()).optional(),
    location_fields: z.array(z.string()),
    default_tracks: z.array(z.record(z.string(), z.unknown())).optional(),
    atac_bam_track: ATACBAMTrackSchema.optional()
});

// Rows as columns subgroup schema
export const RowsAsColumnsSubgroupSchema = z.object({
    name: z.string(),
    type: z.string().optional(),
    label: z.string()
});

// Rows as columns link schema
export const RowsAsColumnsLinkSchema = z.object({
    name_column: z.string(),
    name: z.string(),
    subgroups: z.record(z.string(), RowsAsColumnsSubgroupSchema)
});

// Sync column colors link schema
export const SyncColumnColorsLinkSchema = z.array(z.object({
    link_to: z.string(),
    col: z.string()
}));

// Valueset link schema
export const ValuesetLinkSchema = z.object({
    source_column: z.string(),
    dest_column: z.string()
});

// Link schema - union of all possible link types
export const LinkSchema = z.object({
    interactions: InteractionsConfigSchema.optional(),
    rows_as_columns: RowsAsColumnsLinkSchema.optional(),
    columns: z.array(z.string()).optional(),
    index: z.string().optional(),
    access_data: z.boolean().optional(),
    sync_column_colors: SyncColumnColorsLinkSchema.optional(),
    valueset: ValuesetLinkSchema.optional()
});

// Links configuration schema
export const LinksConfigSchema = z.record(z.string(), LinkSchema);

// Main datasource schema
export const DataSourceSchema = z.object({
    name: z.string(),
    size: z.number().int().positive(),
    columns: z.array(DataSourceColumnSchema),
    columnGroups: z.array(ColumnGroupSchema).optional(),
    links: LinksConfigSchema.optional(),
    images: ImagesConfigSchema.optional(),
    large_images: ImagesConfigSchema.optional(),
    regions: RegionsConfigSchema.optional(),
    interactions: InteractionsConfigSchema.optional(),
    offsets: OffsetsConfigSchema.optional(),
    genome_browser: GenomeBrowserConfigSchema.optional(),
    tree_diagram: z.record(z.string(), z.unknown()).optional(),
    avivator: z.boolean().optional(),
    row_data_loader: z.boolean().optional(),
    binary_data_loader: z.boolean().optional(),
    deeptools: z.record(z.string(), z.unknown()).optional()
});

// Array of datasources schema
export const DataSourcesArraySchema = z.array(DataSourceSchema);

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