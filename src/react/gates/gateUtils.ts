import type { FeatureCollection, Polygon } from "@turf/helpers";
import type { Gate } from "./types";
import { v4 as uuid } from "uuid";

export const DEFAULT_GATE_COLOR: [number, number, number] = [76, 175, 80];
export const SELECTION_FILL_COLOR: [number, number, number, number] = [200, 200, 200, 50];
export const SELECTION_LINE_COLOR: [number, number, number, number] = [150, 150, 150, 255];

/**
 * Extract the coordinates of the gate selection from it's feature collection
 *
 * @param featureCollection - Feature Collection of the gate
 * @returns - Extracted coordinates of the gate
 */
export function extractCoords(featureCollection: FeatureCollection): [number, number][] {
    const feature = featureCollection.features[0];
    if (!feature || !feature.geometry || feature.geometry.type !== "Polygon") {
        return [];
    }

    const polygon = feature.geometry as Polygon;
    if (!polygon.coordinates || polygon.coordinates.length === 0) {
        return [];
    }

    const coords = polygon.coordinates[0];
    if (!Array.isArray(coords) || coords.length === 0) {
        return [];
    }

    return coords.map((coord: number[]) => [coord[0], coord[1]] as [number, number]);
}

/**
 * Check if the point is in the polygon or not
 * (Taken from the RangeDimension.filterPoly)
 */
export function isPointInPolygon(point: [number, number], polygon: [number, number][]): boolean {
    const [x, y] = point;
    let minX = Number.MAX_VALUE;
    let minY = Number.MAX_VALUE;
    let maxX = Number.NEGATIVE_INFINITY;
    let maxY = Number.NEGATIVE_INFINITY;

    for (const [px, py] of polygon) {
        minX = Math.min(minX, px);
        minY = Math.min(minY, py);
        maxX = Math.max(maxX, px);
        maxY = Math.max(maxY, py);
    }

    if (x < minX || x > maxX || y < minY || y > maxY || !Number.isFinite(x) || !Number.isFinite(y)) {
        return false;
    }

    let inside = false;

    for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
        const [xi, yi] = polygon[i];
        const [xj, yj] = polygon[j];

        const intersect = yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi;
        if (intersect) inside = !inside;
    }
    return inside;
}

/**
 * Area-weighted centroid of a polygon using the shoelace formula.
 * Uses the polygon's actual shape and area so the result is stable and lies
 * inside the polygon; ignores closing duplicate vertex and uneven vertex density.
 * @param geometry - Feature collection
 * @returns [cx, cy] or [0, 0] if empty or degenerate (zero area)
 */
export function computeCentroid(geometry: FeatureCollection): [number, number] {
    const polygonCoords = extractCoords(geometry);
    if (polygonCoords.length === 0) return [0, 0];

    // LLM generated shoelace formula code for computing centroid
    const n = polygonCoords.length;
    let signedArea = 0;
    let cx = 0;
    let cy = 0;

    for (let i = 0; i < n; i++) {
        const j = (i + 1) % n;
        const [xi, yi] = polygonCoords[i];
        const [xj, yj] = polygonCoords[j];
        const cross = xi * yj - xj * yi;
        signedArea += cross;
        cx += (xi + xj) * cross;
        cy += (yi + yj) * cross;
    }

    signedArea *= 0.5;
    if (Math.abs(signedArea) < 1e-10) return [0, 0];
    cx /= 6 * signedArea;
    cy /= 6 * signedArea;
    return [cx, cy] as [number, number];
}

export function generateGateId(): string {
    return `gate-${uuid()}-${Date.now()}`;
}

type RelevantGatesArgs = {
    gates: Gate[];
    xField?: string;
    yField?: string;
    region?: string;
};

export function getRelevantGates({
    gates,
    xField,
    yField,
    region,
}: RelevantGatesArgs): Gate[] {
    if (!xField || !yField) return [];

    return gates.filter((gate) => {
        const sameAxes = gate.columns[0] === xField && gate.columns[1] === yField;
        if (!sameAxes) return false;

        if (region) {
            // For viv plot, include both global gates and region-scoped gates for the active region.
            return gate.region === undefined || gate.region === region;
        }

        // For deck scatterplot, include only global gates.
        return gate.region === undefined;
    });
}
