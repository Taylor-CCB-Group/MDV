import { Matrix4 } from "@math.gl/core";
import { micrometresPer } from "intraspatial/spatial/ngffTransform";
import { type AnalysisFrame, frameFromMatrix } from "./analysisFrame";

/** The bits of an sd.js `BaseTransformation` this needs. */
type Transformation = {
    toMatrix(): Matrix4;
    readonly output?: { name: string; axes?: readonly { name: string; unit?: string }[] };
};
type TransformResult = { ok: true; value: Transformation } | { ok: false; error: unknown };
/** The bits of an sd.js `SpatialData` this needs (structural, so tests need no zarr store). */
export type SpatialDataTransforms = {
    images?: Record<string, { getTransformation(cs?: string): TransformResult }>;
};

/** An MDV region entry as the SpatialData converter writes it. */
export type ConvertedRegion = {
    spatial?: { coordinate_system?: string; file?: string };
    viv_image?: { file?: string };
};

/**
 * The image element whose pixel frame the converter wrote the table's x/y in.
 *
 * `python/mdvtools/spatial/conversion.py` `_transform_table_coordinates` maps `obsm["spatial"]`
 * through the annotated element → image transform (point_transform modes image/auto/xenium), so
 * the datastore's position columns are that image's intrinsic pixels, and `viv_image.file` is
 * `<sdata>/images/<key>`. (`--point-transform identity` leaves them in the annotated element's frame
 * instead, and the mode is not recorded in datasources.json — see notes.)
 */
export function regionImageKey(region: ConvertedRegion): string | undefined {
    const file = region.viv_image?.file;
    const sd = region.spatial?.file;
    if (!file || !sd) return undefined;
    const prefix = `${sd}/images/`;
    return file.startsWith(prefix) ? file.slice(prefix.length) : undefined;
}

/**
 * Data columns → µm in the region's coordinate system: the region image's own SpatialData
 * transformation to `spatial.coordinate_system` (any NGFF scale/translation/affine/sequence, via
 * sd.js `toMatrix()`), followed by the coordinate system's unit as µm.
 *
 * Returns an error string when the frame cannot be established, rather than guessing a scale: a
 * length in the wrong unit looks just as plausible as the right one.
 */
export function spatialDataRegionFrame(
    sdata: SpatialDataTransforms,
    region: ConvertedRegion,
): AnalysisFrame | { error: string } {
    const cs = region.spatial?.coordinate_system;
    const key = regionImageKey(region);
    if (!cs || !key) return { error: "region has no SpatialData coordinate system / image" };
    const image = sdata.images?.[key];
    if (!image) return { error: `image element "${key}" not in the SpatialData store` };
    const t = image.getTransformation(cs);
    if (!t.ok) return { error: `image "${key}" has no transformation to "${cs}"` };
    const unit = t.value.output?.axes?.find((a) => a.name === "x")?.unit;
    const um = micrometresPer(unit);
    const toCs = t.value.toMatrix();
    const provenance = `sd.js images/${key} → ${cs}`;
    if (um === undefined) {
        // SpatialData writes the literal "unit" when none was given; "pixel" is not a length either.
        const stated = unit && unit !== "unit" ? unit : "unit unstated";
        return frameFromMatrix(toCs, `${cs} units (${stated})`, provenance);
    }
    return frameFromMatrix(new Matrix4().scale([um, um, um]).multiplyRight(toCs), "µm", provenance);
}
