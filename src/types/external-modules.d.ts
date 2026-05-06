/**
 * Temporary module shims for tsgo + pnpm resolution on this branch.
 *
 * Why this exists:
 * - Turf v6 packages provide declarations, but TypeScript fails to resolve them
 *   cleanly via package exports in our current setup.
 * - pako currently has no resolved declaration file in this workspace.
 *
 * Follow-up cleanup:
 * - Prefer importing shared geometry types from `geojson` at call sites.
 * - Try `pnpm add -D @types/pako` and remove the pako ambient module if it works.
 * - Revisit when Turf is upgraded and exports/type resolution improves.
 */
declare module "pako";

// Minimal surface used by this repo; intentionally not a full Turf typing model.
declare module "@turf/helpers" {
    export type Position = number[];

    export type Geometry = {
        type: string;
        coordinates?: any;
        [key: string]: unknown;
    };

    export type Polygon = {
        type: "Polygon";
        coordinates: Position[][];
    };

    export type Feature<G = Geometry, P = Record<string, unknown>> = {
        type: "Feature";
        geometry: G;
        properties: P;
        id?: string | number;
    };

    export type FeatureCollection<G = Geometry, P = Record<string, unknown>> = {
        type: "FeatureCollection";
        features: Array<Feature<G, P>>;
    };

    export type AllGeoJSON = Feature | FeatureCollection | Geometry;

    export function point<P = Record<string, unknown>>(
        coordinates: Position,
        properties?: P,
        options?: { bbox?: number[]; id?: string | number },
    ): Feature<{ type: "Point"; coordinates: Position }, P>;
}

// Mirrors our current use: clone returns same geojson shape.
declare module "@turf/clone" {
    import type { AllGeoJSON } from "@turf/helpers";
    declare function clone<T extends AllGeoJSON>(geojson: T): T;
    export default clone;
}
