import type { FeatureCollection, Polygon } from '@turf/helpers';
import type { Gate } from './types';
import { v4 as uuid } from 'uuid';

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
export function isPointInPolygon(
    point: [number, number],
    polygon: [number, number][]
): boolean {
    const [x, y] = point;
    let minX = Number.MAX_VALUE;
    let minY = Number.MAX_VALUE;
    let maxX = Number.MIN_VALUE;
    let maxY = Number.MIN_VALUE;

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

export function isPointInGate(x: number, y: number, gate: Gate): boolean {
    const polygonCoords = extractCoords(gate.geometry);
    if (polygonCoords.length === 0) return false;
    return isPointInPolygon([x, y], polygonCoords);
}

/**
 * Calculate the average of the coordinates to get the centroid of the gate
 */
export function computeCentroid(geometry: FeatureCollection): [number, number] {
    const polygonCoords = extractCoords(geometry);
    if (polygonCoords.length === 0) return [0, 0];

    let sumX = 0;
    let sumY = 0;
    for (const [x, y] of polygonCoords) {
        sumX += x;
        sumY += y;
    }

    return [sumX / polygonCoords.length, sumY / polygonCoords.length] as [number, number];
}

export function generateGateId(): string {
    return `gate-${uuid()}-${Date.now()}`;
}

