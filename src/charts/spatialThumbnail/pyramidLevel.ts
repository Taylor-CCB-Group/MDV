/** One OME-Zarr multiscale level, finest first (same order as `loadOmeZarr`). */
export type PyramidLevelExtent = {
    width: number;
    height: number;
};

/**
 * Finest level whose longer side is at most `textureSize`.
 * If every level is larger, the coarsest (last) level.
 */
export function pickPyramidLevel(levels: PyramidLevelExtent[], textureSize: number): number {
    if (levels.length === 0) return -1;
    for (let i = 0; i < levels.length; i++) {
        const longer = Math.max(levels[i].width, levels[i].height);
        if (longer <= textureSize) return i;
    }
    return levels.length - 1;
}

export function sourceExtent(source: {
    shape?: number[];
    labels?: string[];
    width?: number;
    height?: number;
}): PyramidLevelExtent {
    const labels = source.labels;
    const shape = source.shape;
    if (labels && shape) {
        const x = labels.indexOf("x");
        const y = labels.indexOf("y");
        if (x >= 0 && y >= 0) {
            return { width: shape[x], height: shape[y] };
        }
    }
    return { width: source.width ?? 0, height: source.height ?? 0 };
}

/** Nearest-neighbor scale of one channel plane into a square thumbnail. */
export function scalePlane(
    data: ArrayLike<number>,
    srcWidth: number,
    srcHeight: number,
    size: number,
): Float32Array {
    const out = new Float32Array(size * size);
    if (srcWidth <= 0 || srcHeight <= 0 || data.length === 0) return out;
    for (let y = 0; y < size; y++) {
        const sy = Math.min(srcHeight - 1, Math.floor(((y + 0.5) * srcHeight) / size));
        for (let x = 0; x < size; x++) {
            const sx = Math.min(srcWidth - 1, Math.floor(((x + 0.5) * srcWidth) / size));
            out[y * size + x] = Number(data[sy * srcWidth + sx]) || 0;
        }
    }
    return out;
}

export function planeDomain(data: ArrayLike<number>): [number, number] {
    let min = Infinity;
    let max = -Infinity;
    for (let i = 0; i < data.length; i++) {
        const v = data[i];
        if (!Number.isFinite(v)) continue;
        if (v < min) min = v;
        if (v > max) max = v;
    }
    if (!Number.isFinite(min) || !Number.isFinite(max)) return [0, 1];
    if (min === max) return [min, min + 1];
    return [min, max];
}

function clamp01(v: number) {
    if (v < 0) return 0;
    if (v > 1) return 1;
    return v;
}

/** Same bias/gain as VivContrastExtension. Tone is clamped off 0 and 1 so log stays finite. */
export function applyBrightnessContrast(intensity: number, brightness: number, contrast: number): number {
    const b = Math.min(0.99, Math.max(0.01, brightness));
    const g = Math.min(0.99, Math.max(0.01, contrast));
    const bias = (biasValue: number, t: number) => Math.pow(t, Math.log(biasValue) / Math.log(0.5));
    const gain = (gainValue: number, t: number) => {
        if (t < 0.5) return bias(1 - gainValue, 2 * t) / 2;
        return 1 - bias(1 - gainValue, 2 - 2 * t) / 2;
    };
    return bias(b, gain(g, clamp01(intensity)));
}

export type ChannelComposite = {
    plane: Float32Array;
    color: [number, number, number];
    visible: boolean;
    contrastLimits: [number, number];
    brightness: number;
    contrast: number;
};

/** Additive color composite of channel planes into one RGBA thumbnail. */
export function compositePlanes(size: number, channels: ChannelComposite[]): Uint8ClampedArray {
    const rgba = new Uint8ClampedArray(size * size * 4);
    for (let i = 0; i < size * size; i++) {
        let r = 0;
        let g = 0;
        let b = 0;
        for (const channel of channels) {
            if (!channel.visible) continue;
            const [lo, hi] = channel.contrastLimits;
            const raw = channel.plane[i] ?? 0;
            const span = hi - lo;
            const limited = span === 0 ? 0 : clamp01((raw - lo) / span);
            const intensity = applyBrightnessContrast(limited, channel.brightness, channel.contrast);
            r += channel.color[0] * intensity;
            g += channel.color[1] * intensity;
            b += channel.color[2] * intensity;
        }
        const o = i * 4;
        rgba[o] = Math.min(255, r);
        rgba[o + 1] = Math.min(255, g);
        rgba[o + 2] = Math.min(255, b);
        rgba[o + 3] = 255;
    }
    return rgba;
}
