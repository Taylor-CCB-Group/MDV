/**
 * Synchronous, CPU-side picks for the selection (editable GeoJSON) layer.
 *
 * editable-layers decides what a press does (translate, move a vertex, start drawing) from
 * the picks at that instant, and asks for them synchronously. deck cannot pick synchronously
 * on WebGPU, and an async pick can only answer for where the pointer *was*: a quick move then
 * press would draw a new shape instead of moving the existing one. The selection layer holds
 * a handful of vertices, so its picks are cheap to compute exactly, in screen space.
 *
 * The picks have the fields the edit modes read: `object`, `index`, and `isGuide` for edit
 * handles (whose `object.properties` carry `guideType`, `featureIndex`, `positionIndexes`).
 * Edit handles come first, nearest first, as the layer draws them on top.
 */

type Position = readonly number[];
type Geometry = { type: string; coordinates?: unknown };
type Feature = { type: "Feature"; geometry: Geometry | null; properties?: Record<string, unknown> | null };

export type CpuPick = { index: number; object: Feature; isGuide: boolean; picked: true };

/** Screen (CSS pixel) position of a coordinate, or null if it cannot be projected. */
export type Project = (position: Position) => readonly number[] | null;

function isPosition(value: unknown): value is Position {
    return Array.isArray(value) && value.length >= 2 && typeof value[0] === "number" && typeof value[1] === "number";
}

/** The polygons of a Polygon or MultiPolygon, each as rings of positions. */
function polygonsOf(geometry: Geometry | null): Position[][][] {
    if (!geometry || !Array.isArray(geometry.coordinates)) return [];
    if (geometry.type === "Polygon") return [geometry.coordinates.filter(Array.isArray)];
    if (geometry.type === "MultiPolygon") return geometry.coordinates.filter(Array.isArray);
    return [];
}

function project2(project: Project, position: Position): [number, number] | null {
    const p = project(position);
    return p && Number.isFinite(p[0]) && Number.isFinite(p[1]) ? [p[0], p[1]] : null;
}

/** Distance from p to the segment a–b. */
function segmentDistance(p: readonly number[], a: readonly number[], b: readonly number[]): number {
    const dx = b[0] - a[0];
    const dy = b[1] - a[1];
    const lengthSq = dx * dx + dy * dy;
    const t = lengthSq > 0 ? Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / lengthSq)) : 0;
    return Math.hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy));
}

/** Inside (even-odd, so holes count) or within `radius` of an edge, in screen space. */
function hitsPolygon(rings: Position[][], point: readonly number[], project: Project, radius: number): boolean {
    let inside = false;
    for (const ring of rings) {
        const positions = ring.filter(isPosition);
        const pts = positions
            .map((position) => project2(project, position))
            .filter((p): p is [number, number] => p !== null);
        if (pts.length !== positions.length) continue;
        for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
            const [xi, yi] = pts[i];
            const [xj, yj] = pts[j];
            if (segmentDistance(point, pts[j], pts[i]) <= radius) return true;
            if (yi > point[1] !== yj > point[1] && point[0] < ((xj - xi) * (point[1] - yi)) / (yj - yi) + xi) {
                inside = !inside;
            }
        }
    }
    return inside;
}

export function cpuEditablePicks({
    features,
    guides,
    point,
    project,
    radius,
}: {
    features: readonly Feature[];
    guides: readonly Feature[];
    point: readonly number[];
    project: Project;
    radius: number;
}): CpuPick[] {
    const handles: Array<{ pick: CpuPick; distance: number }> = [];
    guides.forEach((guide, index) => {
        const coordinates = guide.geometry?.type === "Point" ? guide.geometry.coordinates : null;
        if (!isPosition(coordinates) || guide.properties?.guideType !== "editHandle") return;
        const p = project2(project, coordinates);
        if (!p) return;
        const distance = Math.hypot(p[0] - point[0], p[1] - point[1]);
        if (distance <= radius) handles.push({ pick: { index, object: guide, isGuide: true, picked: true }, distance });
    });
    handles.sort((a, b) => a.distance - b.distance);
    const picks = handles.map(({ pick }) => pick);
    features.forEach((feature, index) => {
        if (polygonsOf(feature.geometry).some((rings) => hitsPolygon(rings, point, project, radius))) {
            picks.push({ index, object: feature, isGuide: false, picked: true });
        }
    });
    return picks;
}
