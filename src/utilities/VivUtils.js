import { getChannelStats } from "@hms-dbmi/viv";
import { Matrix4 } from "@math.gl/core";

/** copied (not quite verbatim) from avivator utils */
export async function getSingleSelectioStats3D(
    loader,
    selection = { c: 0, t: 0 },
) {
    const lowResSource = loader[loader.length - 1];
    const { shape, labels } = lowResSource;
    const sizeZ = shape[labels.indexOf("z")] >> (loader.length - 1);
    const raster0 = await lowResSource.getRaster({
        selection: { ...selection, z: 0 },
    });
    const rasterMid = await lowResSource.getRaster({
        selection: { ...selection, z: Math.floor(sizeZ / 2) },
    });
    const rasterTop = await lowResSource.getRaster({
        selection: { ...selection, z: Math.max(0, sizeZ - 1) },
    });
    const stats0 = getChannelStats(raster0.data);
    const statsMid = getChannelStats(rasterMid.data);
    const statsTop = getChannelStats(rasterTop.data);
    return {
        domain: [
            Math.min(stats0.domain[0], statsMid.domain[0], statsTop.domain[0]),
            Math.max(stats0.domain[1], statsMid.domain[1], statsTop.domain[1]),
        ],
        contrastLimits: [
            Math.min(
                stats0.contrastLimits[0],
                statsMid.contrastLimits[0],
                statsTop.contrastLimits[0],
            ),
            Math.max(
                stats0.contrastLimits[1],
                statsMid.contrastLimits[1],
                statsTop.contrastLimits[1],
            ),
        ],
    };
}

export async function getMultiSelectionStats(
    loader,
    selections = [{ c: 0, t: 0 }],
) {
    const stats = await Promise.all(
        selections.map((selection) =>
            getSingleSelectioStats3D(loader, selection),
        ),
    );
    const domains = stats.map((stat) => stat.contrastLimits);
    const contrastLimits = stats.map((stat) => stat.contrastLimits);
    return { domains, contrastLimits };
}

export function getDefaultSelectionStats(n) {
    const domains = new Array(n).fill([0, 1000]);
    const contrastLimits = domains;
    //nb, not clear that adding _id here is necessarily ideal, seems to be working for now
    //key for keeping track of multiple 'channels' with selections on the same 'c'.
    const selections = new Array(n).fill().map((_, i) => {
        return { c: i, t: 0, z: 0, _id: i };
    });
    const colors = getDefaultChannelColors(n);
    const channelsVisible = new Array(n).fill(true);
    return { domains, contrastLimits, selections, colors, channelsVisible };
}

export function getDefaultChannelColors(n) {
    if (n === 1) return [[255, 255, 255]];
    return new Array(n).fill([0, 0, 0]).map((_, i) => {
        //TODO: non-shit algorithm / use a nice palette.
        const a = i / n;
        return [Math.floor(a * 255), Math.floor((1 - a) * 255), 0];
    });
}

/**
 * Get physical size scaling Matrix4
 * @param {Object} loader PixelSource
 */
export function getPhysicalSizeScalingMatrix(loader) {
    const { x, y, z } = loader?.meta?.physicalSizes ?? {};
    if (x?.size && y?.size && z?.size) {
        const min = Math.min(z.size, x.size, y.size);
        const ratio = [x.size / min, y.size / min, z.size / min];
        return new Matrix4().scale(ratio);
    }
    return new Matrix4().identity();
}

export function getBoundingCube(loader) {
    const source = Array.isArray(loader) ? loader[0] : loader;
    const { shape, labels } = source;
    const physicalSizeScalingMatrix = getPhysicalSizeScalingMatrix(source);
    const xSlice = [
        0,
        physicalSizeScalingMatrix[0] * shape[labels.indexOf("x")],
    ];
    const ySlice = [
        0,
        physicalSizeScalingMatrix[5] * shape[labels.indexOf("y")],
    ];
    const zSlice = [
        0,
        physicalSizeScalingMatrix[10] * shape[labels.indexOf("z")],
    ];
    return [xSlice, ySlice, zSlice];
}
