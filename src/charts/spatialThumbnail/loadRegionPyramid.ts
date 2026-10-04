import { loadOmeZarr } from "@hms-dbmi/viv";
import { loadOmeZarrMultiscalesFromStore, type VivCompatiblePixelSource } from "zarrextra";
import { pickPyramidLevel, scalePlane, sourceExtent, type PyramidLevelExtent } from "./pyramidLevel";

export type LoadedPyramid = {
    channelNames: string[];
    /** OME channel color, or null when the store does not specify one. */
    channelColors: Array<[number, number, number] | null>;
    /** Scaled planes keyed by channel index `c`. */
    planes: Map<number, Float32Array>;
    level: PyramidLevelExtent;
};

type RasterSelection = { c: number; z?: number; t?: number };

type RasterSource = {
    shape?: number[];
    labels?: string[];
    width?: number;
    height?: number;
    getRaster: (args: { selection: RasterSelection }) => Promise<{
        data: ArrayLike<number>;
        width?: number;
        height?: number;
    }>;
};

/** Include z/t only when the multiscale actually has those axes. Zarr v3 SpatialData images are often c,y,x. */
function rasterSelection(source: { labels?: string[] }, c: number): RasterSelection {
    const labels = new Set((source.labels ?? ["c", "z", "t"]).map((label) => label.toLowerCase()));
    const selection: RasterSelection = { c };
    if (labels.has("z")) selection.z = 0;
    if (labels.has("t")) selection.t = 0;
    return selection;
}

/**
 * SpatialData image groups are Zarr v3 (`zarr.json`). Viv's `loadOmeZarr` only opens Zarr v2 (`.zattrs`).
 * `zarrextra` reads both and returns the same multiscale pixel sources.
 */
function fetchStore(url: string) {
    const root = url.replace(/\/$/, "");
    return {
        async get(key: string) {
            const path = key.startsWith("/") ? key : `/${key}`;
            const response = await fetch(`${root}${path}`);
            if (!response.ok) return undefined;
            return new Uint8Array(await response.arrayBuffer());
        },
    };
}

type OmeChannel = { label?: string; color?: unknown; Color?: unknown };

function parseChannelColor(value: unknown): [number, number, number] | null {
    if (Array.isArray(value) && value.length >= 3) {
        const rgb = value.slice(0, 3).map((part) => Number(part));
        if (rgb.every((part) => Number.isFinite(part))) return [rgb[0], rgb[1], rgb[2]];
    }
    if (typeof value !== "string") return null;
    const hex = value.trim().replace(/^#/, "");
    if (!/^[0-9a-fA-F]{6}$/.test(hex)) return null;
    return [Number.parseInt(hex.slice(0, 2), 16), Number.parseInt(hex.slice(2, 4), 16), Number.parseInt(hex.slice(4, 6), 16)];
}

async function channelMetaFromZarr(url: string): Promise<{ names: string[]; colors: Array<[number, number, number] | null> }> {
    try {
        const response = await fetch(`${url.replace(/\/$/, "")}/zarr.json`);
        if (!response.ok) return { names: [], colors: [] };
        const meta = (await response.json()) as {
            attributes?: { ome?: { omero?: { channels?: OmeChannel[] } }; omero?: { channels?: OmeChannel[] } };
        };
        const channels = meta.attributes?.ome?.omero?.channels ?? meta.attributes?.omero?.channels ?? [];
        return {
            names: channels.map((channel, index) => channel.label || `Channel ${index + 1}`),
            colors: channels.map((channel) => parseChannelColor(channel.color ?? channel.Color)),
        };
    } catch {
        return { names: [], colors: [] };
    }
}

async function planesFromSources(
    sources: Array<RasterSource | VivCompatiblePixelSource>,
    textureSize: number,
    channels: number[],
): Promise<{ planes: Map<number, Float32Array>; level: PyramidLevelExtent } | null> {
    const extents = sources.map((source) => sourceExtent(source));
    const levelIndex = pickPyramidLevel(extents, textureSize);
    if (levelIndex < 0) return null;
    const source = sources[levelIndex];
    const planes = new Map<number, Float32Array>();
    for (const c of channels) {
        try {
            const raster = await source.getRaster({ selection: rasterSelection(source, c) });
            const width = raster.width ?? extents[levelIndex].width;
            const height = raster.height ?? extents[levelIndex].height;
            planes.set(c, scalePlane(raster.data as ArrayLike<number>, width, height, textureSize));
        } catch (error) {
            console.warn("Thumbnail channel plane failed", c, error);
        }
    }
    return { planes, level: extents[levelIndex] };
}

/**
 * Open one OME-Zarr multiscale, copy the chosen pyramid level for `channels`, then drop the loader.
 * Returns null when the URL is not an OME-Zarr pyramid.
 */
export async function loadRegionPyramid(
    url: string,
    textureSize: number,
    channels: number[],
): Promise<LoadedPyramid | null> {
    const meta = await channelMetaFromZarr(url);
    try {
        const sources = await loadOmeZarrMultiscalesFromStore(fetchStore(url));
        const loaded = await planesFromSources(sources, textureSize, channels);
        if (loaded && loaded.planes.size > 0) {
            return { channelNames: meta.names, channelColors: meta.colors, ...loaded };
        }
    } catch (error) {
        console.warn("Zarr thumbnail read failed", url, error);
    }

    let opened: { data?: RasterSource[]; metadata?: { omero?: { channels?: OmeChannel[] } } };
    const abort = new AbortController();
    const timer = setTimeout(() => abort.abort(), 20000);
    try {
        opened = await loadOmeZarr(url, { type: "multiscales", fetchOptions: { signal: abort.signal } });
    } catch {
        return null;
    } finally {
        clearTimeout(timer);
    }
    const sources = opened.data;
    if (!sources?.length) return null;
    const loaded = await planesFromSources(sources, textureSize, channels);
    if (!loaded) return null;
    const vivChannels = opened.metadata?.omero?.channels ?? [];
    const vivNames = vivChannels.map((channel, index) => channel.label || `Channel ${index + 1}`);
    const vivColors = vivChannels.map((channel) => parseChannelColor(channel.color ?? channel.Color));
    return {
        channelNames: meta.names.length ? meta.names : vivNames,
        channelColors: meta.colors.length ? meta.colors : vivColors,
        ...loaded,
    };
}
