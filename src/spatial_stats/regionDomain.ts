import { rasterLevelsFromTree, type SpatialData } from "@spatialdata/core";
import { type Domain, polygonDomain, rectDomain } from "intraspatial/spatial/domain";
import { getChildGroup } from "zarrextra";
import type DataStore from "@/datastore/DataStore";
import { type AnalysisFrame, affine2, applyAffine2 } from "./analysisFrame";
import type { DomainSource } from "./domainSource";
import { readRegionStore } from "./regionSpatialData";
import { type ConvertedRegion, regionImageKey } from "./spatialDataFrame";

/** Axis-aligned bounds in the data-column frame. */
export type Bounds = { minX: number; minY: number; maxX: number; maxY: number };

/** Level-0 pixel extent from an NGFF array shape: spatial axes come last, ordered (z,) y, x. */
export function imageExtentFromShape(shape: readonly number[]): Bounds | undefined {
    const w = shape.at(-1);
    const h = shape.at(-2);
    if (!w || !h) return undefined;
    return { minX: 0, minY: 0, maxX: w, maxY: h };
}

/**
 * The extent of the region's image in the data-column frame — its pixels, since that is the frame the
 * converter writes x/y in (see `regionImageKey`). Read from the store's consolidated metadata; no
 * array data is fetched.
 */
export function regionImageExtent(sdata: SpatialData, region: ConvertedRegion): Bounds | undefined {
    const key = regionImageKey(region);
    const image = key ? sdata.images?.[key] : undefined;
    const node = key ? getChildGroup(sdata.rootStore.tree, "images", key) : undefined;
    if (!image || !node) return undefined;
    const level0 = rasterLevelsFromTree(node, image.scaleLevels)[0];
    return level0 ? imageExtentFromShape(level0.shape) : undefined;
}

/** A data-frame rectangle through the frame: a `rectDomain` when the frame keeps it axis-aligned,
 *  otherwise the exact quadrilateral as a polygon. */
export function boundsDomain(b: Bounds, frame: AnalysisFrame): Domain {
    const m = affine2(frame);
    const corners = [
        applyAffine2(m, b.minX, b.minY),
        applyAffine2(m, b.maxX, b.minY),
        applyAffine2(m, b.maxX, b.maxY),
        applyAffine2(m, b.minX, b.maxY),
    ];
    if (m[1] !== 0 || m[2] !== 0) return polygonDomain(corners);
    const xs = corners.map((c) => c[0]);
    const ys = corners.map((c) => c[1]);
    return rectDomain([Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)]);
}

/** Where a region's whole-extent bounds come from. */
export type RegionBounds = { bounds: Bounds; from: string };

/** The region image's extent, else the range of the position columns (all regions' rows). */
export async function resolveRegionBounds(
    dataStore: DataStore,
    columns: readonly [string, string],
    region?: string,
): Promise<RegionBounds> {
    const { entry, sdata } = await readRegionStore(dataStore, region);
    const image = entry && sdata ? regionImageExtent(sdata, entry) : undefined;
    if (image && entry) return { bounds: image, from: `image ${regionImageKey(entry)} extent` };
    const [mx, my] = columns.map((c) => dataStore.columnIndex[c]?.minMax);
    if (!mx || !my) throw new Error(`No extent for ${columns.join("/")}: no region image and no column range`);
    return { bounds: { minX: mx[0], minY: my[0], maxX: mx[1], maxY: my[1] }, from: "x/y column range" };
}

/**
 * The whole region — its coordinate-system extent — as the observation window: the default when there
 * is no selection or gate. Cells are restricted to the region's rows.
 */
export function regionDomainSource(dataStore: DataStore, region?: string): DomainSource {
    const columns = positionColumns(dataStore);
    const label = `whole region${region ? ` "${region}"` : ""}`;
    return {
        label,
        columns,
        region,
        toDomain: async (frame) => boundsDomain((await resolveRegionBounds(dataStore, columns, region)).bounds, frame),
    };
}

export function positionColumns(dataStore: DataStore): readonly [string, string] {
    const [x, y] = dataStore.regions?.position_fields ?? [];
    if (!x || !y) throw new Error("Datastore has no regions.position_fields");
    return [x, y];
}
