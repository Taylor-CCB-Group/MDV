export type ThumbnailCenter = { x: number; y: number };

/**
 * Place filtered cell centroids into thumbnail squares.
 * Raster row 0 is the top of each image, so the smaller y sits at the top of the square.
 * `regionCode` null means every row belongs to the single thumbnail.
 * A filter value other than 0 is hidden, matching `DataStore.filterArray`.
 */
export function layoutCellPositions(input: {
    x: ArrayLike<number>;
    y: ArrayLike<number>;
    regionCode: ArrayLike<number> | null;
    slotCodes: number[];
    filter: ArrayLike<number> | null;
    centers: ReadonlyArray<ThumbnailCenter>;
    size: number;
}): Float32Array {
    const n = Math.min(input.x.length, input.y.length);
    const slots = input.centers.length;
    if (slots === 0 || n === 0) return new Float32Array(0);
    const useRegion = input.regionCode !== null;
    if (!useRegion && slots !== 1) return new Float32Array(0);

    const count = new Uint32Array(slots);
    const minX = new Float64Array(slots);
    const maxX = new Float64Array(slots);
    const minY = new Float64Array(slots);
    const maxY = new Float64Array(slots);
    minX.fill(Number.POSITIVE_INFINITY);
    maxX.fill(Number.NEGATIVE_INFINITY);
    minY.fill(Number.POSITIVE_INFINITY);
    maxY.fill(Number.NEGATIVE_INFINITY);

    const slotOf = new Int32Array(n);
    slotOf.fill(-1);
    const codeToSlot = new Map<number, number>();
    if (useRegion) {
        for (let slot = 0; slot < slots; slot++) {
            const code = input.slotCodes[slot];
            if (code >= 0) codeToSlot.set(code, slot);
        }
    }

    const regionCode = input.regionCode;
    const filter = input.filter;
    for (let i = 0; i < n; i++) {
        if (filter && filter[i] !== 0) continue;
        const x = input.x[i];
        const y = input.y[i];
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        let slot = 0;
        if (useRegion) {
            if (!regionCode) continue;
            const found = codeToSlot.get(regionCode[i]);
            if (found === undefined) continue;
            slot = found;
        }
        slotOf[i] = slot;
        count[slot] += 1;
        if (x < minX[slot]) minX[slot] = x;
        if (x > maxX[slot]) maxX[slot] = x;
        if (y < minY[slot]) minY[slot] = y;
        if (y > maxY[slot]) maxY[slot] = y;
    }

    let total = 0;
    for (let slot = 0; slot < slots; slot++) total += count[slot];
    const positions = new Float32Array(total * 2);
    const cursor = new Uint32Array(slots);
    let offset = 0;
    for (let slot = 0; slot < slots; slot++) {
        cursor[slot] = offset;
        offset += count[slot];
    }

    const size = input.size;
    for (let i = 0; i < n; i++) {
        const slot = slotOf[i];
        if (slot < 0) continue;
        const spanX = maxX[slot] - minX[slot] || 1;
        const spanY = maxY[slot] - minY[slot] || 1;
        const u = (input.x[i] - minX[slot]) / spanX;
        const v = (input.y[i] - minY[slot]) / spanY;
        const center = input.centers[slot];
        const at = cursor[slot] * 2;
        cursor[slot] += 1;
        positions[at] = center.x + (u - 0.5) * size;
        positions[at + 1] = center.y + (0.5 - v) * size;
    }
    return positions;
}
