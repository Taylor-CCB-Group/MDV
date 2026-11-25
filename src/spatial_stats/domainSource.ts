import type { SimpleFeatureCollection } from "@deck.gl-community/editable-layers";
import { type Domain, polygonDomain } from "intraspatial/spatial/domain";
import type { Gate } from "@/react/gates/types";
import { type AnalysisFrame, affine2, applyAffine2 } from "./analysisFrame";

/**
 * Where an observation window comes from — a saved gate, the live lasso/polygon selection, later a
 * mask from a wand or a segmentation — separated from what is computed on it.
 *
 * A source carries its geometry in its own frame and produces the intraspatial `Domain` only when
 * handed the `AnalysisFrame` the cells are also going through, so the two cannot disagree. Sources
 * are plain values: a dataflow node can hold one and recompute when it changes.
 */
export interface DomainSource {
    /** For debug output / UI. */
    readonly label: string;
    /** The position columns the geometry is expressed in (and the cells are read from). */
    readonly columns: readonly [string, string];
    /** Restrict cells to rows whose `regions.region_field` equals this. */
    readonly region?: string;
    /** May load (a region's extent, later a mask) — hence the promise. */
    toDomain(frame: AnalysisFrame): Domain | Promise<Domain>;
}

type Ring = [number, number][];

/**
 * Every ring of every Polygon / MultiPolygon feature, in the geometry's own frame. Holes and
 * disjoint parts are kept; `polygonDomain` treats the set as an even–odd region, so ring order and
 * winding don't matter.
 */
export function featureRings(geometry: SimpleFeatureCollection): Ring[] {
    const rings: Ring[] = [];
    const add = (ring: readonly (readonly number[])[]) => {
        if (ring.length >= 3) rings.push(ring.map(([x, y]) => [x, y]));
    };
    for (const f of geometry.features) {
        const g = f.geometry;
        if (g?.type === "Polygon") for (const ring of g.coordinates) add(ring);
        else if (g?.type === "MultiPolygon") for (const poly of g.coordinates) for (const ring of poly) add(ring);
    }
    return rings;
}

/** Polygon geometry in the data-column frame → `polygonDomain` in the analysis frame. An affine maps
 *  polygons to polygons, so this is exact under rotation and shear too. */
export function polygonGeometryDomain(geometry: SimpleFeatureCollection, frame: AnalysisFrame): Domain {
    const m = affine2(frame);
    const rings = featureRings(geometry).map((ring) => ring.map(([x, y]) => applyAffine2(m, x, y)));
    if (rings.length === 0) throw new Error("Geometry has no polygon rings");
    return polygonDomain(rings);
}

export function geometryDomainSource(
    label: string,
    geometry: SimpleFeatureCollection,
    columns: readonly [string, string],
    region?: string,
): DomainSource {
    return { label, columns, region, toDomain: (frame) => polygonGeometryDomain(geometry, frame) };
}

/** A saved gate. */
export function gateDomainSource(gate: Gate): DomainSource {
    return geometryDomainSource(`gate "${gate.name}"`, gate.geometry, gate.columns, gate.region);
}

/**
 * The live lasso / polygon / rectangle selection of a deck scatter chart: `config.selectionFeatureCollection`,
 * which `unproject` stores in the data-column frame (modelMatrix already inverted), over `config.param`.
 */
export function selectionDomainSource(config: {
    selectionFeatureCollection?: SimpleFeatureCollection;
    param: readonly unknown[];
    region?: string;
    title?: string;
}): DomainSource | undefined {
    const geometry = config.selectionFeatureCollection;
    if (!geometry || featureRings(geometry).length === 0) return undefined;
    const [x, y] = config.param;
    if (typeof x !== "string" || typeof y !== "string") return undefined;
    return geometryDomainSource(`selection in "${config.title ?? "chart"}"`, geometry, [x, y], config.region);
}
