/**
 * Flags image/labels data whose on-disk layout will make the spatial chart slow, so the
 * chart can say "this data needs re-encoding" instead of just being slow.
 *
 * Two layouts hurt, independently of anything the client does:
 *
 * - **No usable pyramid.** Channel stats and the histogram read the coarsest level, and
 *   the zoomed-out view draws from it. When that level is full resolution, every channel
 *   load scans every pixel (viv's `getChannelStats` alone is seconds on the main thread
 *   for a ~60 MP plane) and the overview decodes the whole image.
 * - **Oversized chunks.** A tile decodes every chunk it touches. Full-width strips or
 *   very large chunks mean each tile decodes many times the pixels it shows, and only a
 *   handful fit in the shared chunk cache, so panning re-fetches.
 *
 * Only array metadata is read (`zarr.json` / `.zarray` per level), never pixels.
 */

/** Coarsest-level plane size above which stats and overview rendering get slow. */
export const COARSEST_LEVEL_PIXEL_BUDGET = 2048 * 2048;
/** Decoded size of one chunk above which tiles do far more work than they show. */
export const CHUNK_DECODED_BYTES_BUDGET = 16 * 1024 * 1024;

export type RasterLevelInfo = {
    path: string;
    shape: number[];
    /** Inner (decode-unit) chunk shape: the sharding codec's inner chunks when sharded. */
    chunks: number[];
    bytesPerElement: number;
};

export type RasterIssue =
    | { kind: "no-pyramid"; width: number; height: number }
    | { kind: "coarsest-level-large"; levels: number; width: number; height: number }
    | { kind: "oversized-chunks"; chunkWidth: number; chunkHeight: number; decodedBytes: number };

const product = (values: number[]) => values.reduce((a, b) => a * b, 1);

/** OME-NGFF puts spatial axes last, so y and x are the final two dimensions. */
function planeSize(shape: number[]) {
    return { height: shape.at(-2) ?? 1, width: shape.at(-1) ?? 1 };
}

export function assessRasterLevels(levels: RasterLevelInfo[]): RasterIssue[] {
    if (levels.length === 0) return [];
    const issues: RasterIssue[] = [];
    // Levels are listed finest first; the coarsest is what stats and the overview read.
    const coarsest = levels[levels.length - 1];
    const { width, height } = planeSize(coarsest.shape);
    if (width * height > COARSEST_LEVEL_PIXEL_BUDGET) {
        issues.push(
            levels.length === 1
                ? { kind: "no-pyramid", width, height }
                : { kind: "coarsest-level-large", levels: levels.length, width, height },
        );
    }
    const finest = levels[0];
    const decodedBytes = product(finest.chunks) * finest.bytesPerElement;
    if (decodedBytes > CHUNK_DECODED_BYTES_BUDGET) {
        const chunk = planeSize(finest.chunks);
        issues.push({
            kind: "oversized-chunks",
            chunkWidth: chunk.width,
            chunkHeight: chunk.height,
            decodedBytes,
        });
    }
    return issues;
}

const megapixels = (width: number, height: number) => `${((width * height) / 1e6).toFixed(1)} MP`;
const mebibytes = (bytes: number) => `${Math.round(bytes / (1024 * 1024))} MB`;

export function describeRasterIssue(issue: RasterIssue): string {
    switch (issue.kind) {
        case "no-pyramid":
            return `No resolution pyramid: a single ${issue.width}×${issue.height} level (${megapixels(issue.width, issue.height)}). Channel stats and the zoomed-out view read every pixel.`;
        case "coarsest-level-large":
            return `Pyramid stops too early: the coarsest of ${issue.levels} levels is ${issue.width}×${issue.height} (${megapixels(issue.width, issue.height)}). Channel stats and the zoomed-out view read all of it.`;
        case "oversized-chunks":
            return `Chunks are ${issue.chunkWidth}×${issue.chunkHeight} (${mebibytes(issue.decodedBytes)} decoded each), so every tile decodes far more than it shows.`;
    }
}

export const RASTER_REENCODE_ADVICE =
    "Re-encode with a multiscale pyramid (downsample by 2 until the coarsest level is about 2048 px or smaller) and 512–1024 px chunks.";

const V2_DTYPE_BYTES = /^[<>|][a-z](\d+)$/;
const V3_DTYPE_BYTES: Record<string, number> = {
    bool: 1,
    int8: 1,
    uint8: 1,
    int16: 2,
    uint16: 2,
    float16: 2,
    int32: 4,
    uint32: 4,
    float32: 4,
    int64: 8,
    uint64: 8,
    float64: 8,
};

function isNumberArray(value: unknown): value is number[] {
    return Array.isArray(value) && value.every((v) => typeof v === "number");
}

function isRecord(value: unknown): value is Record<string, unknown> {
    return typeof value === "object" && value !== null;
}

/** Shape, inner chunk shape and element size from a v3 `zarr.json` or v2 `.zarray`. */
export function parseArrayMetadata(path: string, meta: unknown): RasterLevelInfo | null {
    if (!isRecord(meta) || !isNumberArray(meta.shape)) return null;
    // v2
    if (isNumberArray(meta.chunks) && typeof meta.dtype === "string") {
        const bytes = V2_DTYPE_BYTES.exec(meta.dtype)?.[1];
        if (!bytes) return null;
        return { path, shape: meta.shape, chunks: meta.chunks, bytesPerElement: Number(bytes) };
    }
    // v3
    if (typeof meta.data_type !== "string") return null;
    const bytesPerElement = V3_DTYPE_BYTES[meta.data_type];
    const grid = isRecord(meta.chunk_grid) && isRecord(meta.chunk_grid.configuration) ? meta.chunk_grid.configuration : null;
    if (!bytesPerElement || !grid || !isNumberArray(grid.chunk_shape)) return null;
    let chunks = grid.chunk_shape;
    // A sharded array's grid chunk is the shard; the inner chunk is what gets decoded.
    const codecs = Array.isArray(meta.codecs) ? meta.codecs : [];
    for (const codec of codecs) {
        if (
            isRecord(codec) &&
            codec.name === "sharding_indexed" &&
            isRecord(codec.configuration) &&
            isNumberArray(codec.configuration.chunk_shape)
        ) {
            chunks = codec.configuration.chunk_shape;
        }
    }
    return { path, shape: meta.shape, chunks, bytesPerElement };
}

type RasterStore = { get(key: `/${string}`): Promise<Uint8Array | undefined> };
export type RasterElementLike = { scaleLevels: string[]; getStore(): RasterStore };

const decoder = new TextDecoder();

async function readJson(store: RasterStore, key: `/${string}`) {
    const bytes = await store.get(key);
    return bytes ? JSON.parse(decoder.decode(bytes)) : undefined;
}

export async function loadRasterLevels(element: RasterElementLike): Promise<RasterLevelInfo[]> {
    const store = element.getStore();
    const levels = await Promise.all(
        element.scaleLevels.map(async (path) => {
            const meta = (await readJson(store, `/${path}/zarr.json`)) ?? (await readJson(store, `/${path}/.zarray`));
            return parseArrayMetadata(path, meta);
        }),
    );
    return levels.filter((level): level is RasterLevelInfo => level !== null);
}
