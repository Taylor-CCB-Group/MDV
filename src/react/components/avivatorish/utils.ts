import { useState, useEffect } from "react";
import { fromBlob, fromUrl } from "geotiff";
import { Matrix4 } from "@math.gl/core";

import {
    loadOmeTiff,
    loadBioformatsZarr,
    loadOmeZarr,
    loadMultiTiff,
    getChannelStats,
    RENDERING_MODES,
    ColorPalette3DExtensions,
    AdditiveColormap3DExtensions,
} from "@vivjs-experimental/viv";

import { GLOBAL_SLIDER_DIMENSION_FIELDS } from "./constants";
import { isArray } from "@/lib/utils";
import type { OME_TIFF, PixelSource } from "./state";

// Import and register JPEG2000 decoder
import "@/utilities/jpeg2000/Jpeg2000Decoder";

const MAX_CHANNELS_FOR_SNACKBAR_WARNING = 40;

/**
 * Guesses whether string URL or File is for an OME-TIFF image.
 * @param {string | File} urlOrFile
 */
function isOmeTiff(urlOrFile: string | File) {
    if (Array.isArray(urlOrFile)) return false; // local Zarr is array of File Objects
    const name = typeof urlOrFile === "string" ? urlOrFile : urlOrFile.name;
    return (
        name.includes("ome.tiff") ||
        name.includes("ome.tif") ||
        name.includes(".companion.ome")
    );
}
export type UrlOrFiles = string | File | File[];
/**
 * Gets an array of filenames for a multi tiff input.
 * @param {string | File | File[]} urlOrFiles
 */
function getMultiTiffFilenames(urlOrFiles: UrlOrFiles) {
    if (Array.isArray(urlOrFiles)) {
        return urlOrFiles.map((f) => f.name);
    }
    if (urlOrFiles instanceof File) {
        return [urlOrFiles.name];
    }
    return urlOrFiles.split(",");
}

/**
 * Guesses whether string URL or File is one or multiple standard TIFF images.
 */
function isMultiTiff(urlOrFiles: UrlOrFiles) {
    const filenames = getMultiTiffFilenames(urlOrFiles);
    for (const filename of filenames) {
        const lowerCaseName = filename.toLowerCase();
        if (
            !(lowerCaseName.includes(".tiff") || lowerCaseName.includes(".tif"))
        )
            return false;
    }
    return true;
}

/**
 * Turns an input string of one or many urls, file, or file array into a uniform array.
 */
async function generateMultiTiffFileArray(urlOrFiles: UrlOrFiles) {
    if (Array.isArray(urlOrFiles)) {
        return urlOrFiles;
    }
    if (urlOrFiles instanceof File) {
        return [urlOrFiles];
    }
    return urlOrFiles.split(",");
}

/**
 * Gets the basic image count for a TIFF using geotiff's getImageCount.
 */
async function getTiffImageCount(src: string | File) {
    const from = typeof src === "string" ? fromUrl : fromBlob;
    //@ts-ignore - wtf why is the inference not working here?
    const tiff = await from(src);
    return tiff.getImageCount();
}
/** addresses `c` channel, `z` z-stack index, `t`: time index */
export type VivSelection = { c: number; z: number; t: number };
// all very clever but unnecessary abstraction and also wrong:
// type VivSelection = { [Key in typeof GLOBAL_SLIDER_DIMENSION_FIELDS[number]]?: number };
/**
 * Guesses whether string URL or File is one or multiple standard TIFF images.
 */
async function generateMultiTiffSources(urlOrFiles: UrlOrFiles) {
    const multiTiffFiles = await generateMultiTiffFileArray(urlOrFiles);
    const sources: [VivSelection[], any][] = [];
    let c = 0;
    for (const tiffFile of multiTiffFiles) {
        const selections: VivSelection[] = [];
        const numImages = await getTiffImageCount(tiffFile);
        for (let i = 0; i < numImages; i++) {
            selections.push({ c, z: 0, t: 0 });
            c += 1;
        }
        sources.push([selections, tiffFile]);
    }
    return sources;
}

class UnsupportedBrowserError extends Error {
    constructor(message: string) {
        super(message);
        this.name = "UnsupportedBrowserError";
    }
}
/** todo - clarify types / zarr... */
type OmeTiffImage = OME_TIFF;
async function getTotalImageCount(sources: OmeTiffImage[]) {
    const firstOmeTiffImage = sources[0];
    const firstPixelSource = firstOmeTiffImage.data[0];
    //@ts-ignore - this is a private method
    const representativeGeoTiffImage = await firstPixelSource._indexer({
        c: 0,
        z: 0,
        t: 0,
    });
    const hasSubIFDs = Boolean(
        representativeGeoTiffImage?.fileDirectory?.SubIFDs,
    );

    // Non-Bioformats6 pyramids use Image tags for pyramid levels and do not have offsets
    // built in to the format for them, hence the ternary.

    if (hasSubIFDs) {
        return sources.reduce((sum, { metadata }) => {
            const { SizeC, SizeT, SizeZ } = metadata.Pixels;
            const numImagesPerResolution = SizeC * SizeT * SizeZ;
            return numImagesPerResolution + sum;
        }, 1);
    }

    const levels = firstOmeTiffImage.data.length;
    const { SizeC, SizeT, SizeZ } = firstOmeTiffImage.metadata.Pixels;
    const numImagesPerResolution = SizeC * SizeT * SizeZ;
    return numImagesPerResolution * levels;
}

function isZodError(e: unknown): e is Error & { issues: unknown } {
    return e instanceof Error && "issues" in e;
}

async function fetchSingleFileOmeTiffOffsets(url: string) {
    // No offsets for multifile OME-TIFFs
    if (url.includes("companion.ome")) {
        return undefined;
    }
    const offsetsUrl = url.replace(/ome\.tif(f?)/gi, "offsets.json");
    try {
        const res = await fetch(offsetsUrl); //, {mode: 'no-cors'}); todo 2024-08-13: sort out CORS issues
        return res.status === 200 ? await res.json() : undefined;
    } catch {
        console.warn(`Failed to fetch offsets for ${url}`);
        return undefined;
    }
}

/**
 * Given an image source, creates a PixelSource[] and returns XML-meta
 */
export async function createLoader(
    urlOrFile: UrlOrFiles,
    handleOffsetsNotFound: (arg0: boolean) => void,
    handleLoaderError: (msg: string | null) => void,
) {
    // If the loader fails to load, handle the error (show an error snackbar).
    // Otherwise load.
    try {
        // OME-TIFF
        if (!isArray(urlOrFile) && isOmeTiff(urlOrFile)) {
            if (urlOrFile instanceof File) {
                // TODO(2021-05-09): temporarily disable `pool` until inline worker module is fixed.
                const source = await loadOmeTiff(urlOrFile, {
                    images: "all",
                    // pool: false
                });
                return source;
            }

            const maybeOffsets = await fetchSingleFileOmeTiffOffsets(urlOrFile);

            // TODO(2021-05-06): temporarily disable `pool` until inline worker module is fixed.
            const source = await loadOmeTiff(urlOrFile, {
                offsets: maybeOffsets,
                images: "all",
                // todo 2024-08-13: sort out CORS issues
                // headers: {
                //     mode: 'no-cors'
                // }
                // pool: false
            });

            // Show a warning if the total number of channels/images exceeds a fixed amount.
            const totalImageCount = await getTotalImageCount(source);
            if (
                !maybeOffsets &&
                totalImageCount > MAX_CHANNELS_FOR_SNACKBAR_WARNING
            ) {
                handleOffsetsNotFound(true);
            }
            return source;
        }

        if (
            isArray(urlOrFile) &&
            typeof urlOrFile[0].arrayBuffer !== "function"
        ) {
            throw new UnsupportedBrowserError(
                "Cannot upload a local Zarr or flat TIFF files with this browser. Try using Chrome, Firefox, or Microsoft Edge.",
            );
        }

        // Multiple flat tiffs
        if (isMultiTiff(urlOrFile)) {
            const mutiTiffSources = await generateMultiTiffSources(urlOrFile);
            const source = await loadMultiTiff(mutiTiffSources, {
                // images: 'all',
                // pool: false
            });
            return source;
        }

        // Bio-Formats Zarr
        let source;
        try {
            //@ts-ignore !ruh-roh
            source = await loadBioformatsZarr(urlOrFile);
        } catch (e) {
            if (isZodError(e)) {
                // If the error is a ZodError, it means there was an OME-XML file
                // but it was invalid. We shouldn't try to load the file as a OME-Zarr.
                throw e;
            }

            // try ome-zarr
            //@ts-ignore !ruh-roh
            const res = await loadOmeZarr(urlOrFile, { type: "multiscales" });
            // extract metadata into OME-XML-like form
            const metadata = {
                Pixels: {
                    Channels: res.metadata.omero.channels.map((c) => ({
                        Name: c.label,
                        SamplesPerPixel: 1,
                    })),
                },
            };
            source = { data: res.data, metadata };
        }
        return source;
    } catch (e) {
        if (e instanceof UnsupportedBrowserError) {
            handleLoaderError(e.message);
        } else {
            console.error(e); // eslint-disable-line
            handleLoaderError(null);
        }
        throw e; // Re-throw the error to allow the calling component to catch it
        // return { data: null };
    }
}

// Get the last part of a url (minus query parameters) to be used
// as a display name for avivator.
export function getNameFromUrl(url: string) {
    return url.split("?")[0].split("/").slice(-1)[0];
}
//type DIMENSION_FIELDx = typeof GLOBAL_SLIDER_DIMENSION_FIELDS[number];
// type DIMENSION_FIELD = "x" | "y" | "z" | "t" | "c"; //?
type DIMENSION_FIELD = "z" | "t";
/**
 * Return the midpoint of the global dimensions as a default selection.
 * 
 * n.b. some of the more abstract way in which dimension could be addressed has been
 * removed in favor of a more direct approach - would be good to verify that the premise
 * is still valid.
 */
function getDefaultGlobalSelection(dimensions: {name: string, size: number}[]) {
    const globalSelectableDimensions = dimensions.filter((d) =>
        GLOBAL_SLIDER_DIMENSION_FIELDS.includes(d.name.toLowerCase()),
    ) as unknown as DIMENSION_FIELD[];

    const selection: Partial<VivSelection> = {};
    for (const dim of globalSelectableDimensions) {
        //@ts-ignore
        selection[dim.name] = Math.floor(dim.size / 2);
    }

    return selection as VivSelection; //!probably not partial at this point?
}

function isGlobalOrXYDimension(name: string): name is "x" | "y" {
    // normalize name to lowercase
    name = name.toLowerCase();
    return (
        name === "x" ||
        name === "y" ||
        GLOBAL_SLIDER_DIMENSION_FIELDS.includes(name)
    );
}

/**
 * @param shape loader shape
 */
export function isInterleaved(shape: number[]) {
    const lastDimSize = shape[shape.length - 1];
    return lastDimSize === 3 || lastDimSize === 4;
}

function zip<A, B>(a: A[], b: B[]): [A, B][] {
    if (a.length !== b.length) {
        // todo consider using the min length of the two arrays
        throw new Error("Array lengths must be equal");
    }
    return a.map((val, i) => [val, b[i]]);
}
/**
 * Create a default selection using the midpoint of the available global dimensions,
 * and then the first four available selections from the first selectable channel.
 * @param pixelSource
 */
export function buildDefaultSelection({ labels, shape }: { labels: string[], shape: number[] }): VivSelection[] {
    const selection: VivSelection[] = [];

    const dimensions = zip(labels, shape).map(([name, size]) => ({
        name,
        size,
    }));

    const globalSelection = getDefaultGlobalSelection(dimensions);

    // First non-global dimension with some sort of selectable values.
    const firstNonGlobalSelectableDimension = dimensions.find(
        (dim) => !isGlobalOrXYDimension(dim.name),
    );

    // If there are no additional selectable dimensions, return the global selection.
    if (!firstNonGlobalSelectableDimension) {
        return [globalSelection];
    }

    for (
        let i = 0;
        i < Math.min(4, firstNonGlobalSelectableDimension.size);
        i += 1
    ) {
        //@ts-ignore
        selection.push({
            [firstNonGlobalSelectableDimension.name]: i,
            ...globalSelection,
        });
    }

    if (isInterleaved(shape)) {
        return [{ ...selection[0], c: 0 }];
    }

    return selection;
}

export function hexToRgb(hex: string) {
    // https://stackoverflow.com/questions/5623838/rgb-to-hex-and-hex-to-rgb
    const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    if (!result) {
        console.error(`Invalid hex color: ${hex}`);
        return [255, 0, 255];
    }
    return result.map((d) => Number.parseInt(d, 16)).slice(1);
}

export function range(length: number) {
    return [...Array(length).keys()];
}

// consider a version of this that uses outerContainer in mdv
export function useWindowSize(scaleWidth = 1, scaleHeight = 1) {
    function getSize() {
        return {
            width: window.innerWidth * scaleWidth,
            height: window.innerHeight * scaleHeight,
        };
    }
    const [windowSize, setWindowSize] = useState(getSize());
    useEffect(() => {
        const handleResize = () => {
            setWindowSize(getSize());
        };
        window.addEventListener("resize", handleResize);
        return () => {
            window.removeEventListener("resize", handleResize);
        };
    });
    return windowSize;
}
type LOADER = OME_TIFF | any; //idk
function narrowLimits(limits: number[]): limits is [number, number] {
    return limits.length === 2;
}
function narrowStats(stats: { domain: number[]; contrastLimits: number[] }): stats is { domain: [number, number]; contrastLimits: [number, number] } {
    return narrowLimits(stats.domain) && narrowLimits(stats.contrastLimits);
}
export async function getSingleSelectionStats2D({ loader, selection }: { loader: LOADER, selection: VivSelection}) {
    const data = Array.isArray(loader) ? loader[loader.length - 1] : loader;
    const raster = await data.getRaster({ selection });
    const selectionStats = getChannelStats(raster.data);
    if (!narrowStats(selectionStats)) {
        throw new Error("expected getChannelStats() to return domain and contrastLimits as [number, number]");
    }
    const { domain, contrastLimits } = selectionStats;
    // Edge case: if the contrast limits are the same, set them to the domain.
    if (contrastLimits[0] === contrastLimits[1]) {
        contrastLimits[0] = domain[0];
        contrastLimits[1] = domain[1];
    }
    return { domain, contrastLimits };
}

export async function getSingleSelectionStats3D({ loader, selection }: { loader: LOADER, selection: VivSelection }) {
    const lowResSource = loader[loader.length - 1];
    const { shape, labels } = lowResSource;
    const sizeZ = shape[labels.indexOf("z")];
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
        ] satisfies [number, number],
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
        ] satisfies [number, number],
    };
}

export const getSingleSelectionStats = async ({ loader, selection, use3d }: { loader: LOADER, selection: VivSelection, use3d: boolean }) => {
    const getStats = use3d
        ? getSingleSelectionStats3D
        : getSingleSelectionStats2D;
    return getStats({ loader, selection });
};

export const getMultiSelectionStats = async ({ loader, selections, use3d }: { loader: LOADER, selections: VivSelection[], use3d: boolean }) => {
    const stats = await Promise.all(
        selections.map((selection) =>
            getSingleSelectionStats({ loader, selection, use3d }),
        ),
    );
    const domains = stats.map((stat) => stat.domain);
    const contrastLimits = stats.map((stat) => stat.contrastLimits);
    return { domains, contrastLimits };
};

/* eslint-disable no-useless-escape */
// https://stackoverflow.com/a/11381730
export function isMobileOrTablet() {
    let check = false;
    // eslint-disable-next-line func-names
    ((a) => {
        if (
            /(android|bb\d+|meego).+mobile|avantgo|bada\/|blackberry|blazer|compal|elaine|fennec|hiptop|iemobile|ip(hone|od)|iris|kindle|lge |maemo|midp|mmp|mobile.+firefox|netfront|opera m(ob|in)i|palm( os)?|phone|p(ixi|re)\/|plucker|pocket|psp|series(4|6)0|symbian|treo|up\.(browser|link)|vodafone|wap|windows ce|xda|xiino|android|ipad|playbook|silk/i.test(
                a,
            ) ||
            /1207|6310|6590|3gso|4thp|50[1-6]i|770s|802s|a wa|abac|ac(er|oo|s\-)|ai(ko|rn)|al(av|ca|co)|amoi|an(ex|ny|yw)|aptu|ar(ch|go)|as(te|us)|attw|au(di|\-m|r |s )|avan|be(ck|ll|nq)|bi(lb|rd)|bl(ac|az)|br(e|v)w|bumb|bw\-(n|u)|c55\/|capi|ccwa|cdm\-|cell|chtm|cldc|cmd\-|co(mp|nd)|craw|da(it|ll|ng)|dbte|dc\-s|devi|dica|dmob|do(c|p)o|ds(12|\-d)|el(49|ai)|em(l2|ul)|er(ic|k0)|esl8|ez([4-7]0|os|wa|ze)|fetc|fly(\-|_)|g1 u|g560|gene|gf\-5|g\-mo|go(\.w|od)|gr(ad|un)|haie|hcit|hd\-(m|p|t)|hei\-|hi(pt|ta)|hp( i|ip)|hs\-c|ht(c(\-| |_|a|g|p|s|t)|tp)|hu(aw|tc)|i\-(20|go|ma)|i230|iac( |\-|\/)|ibro|idea|ig01|ikom|im1k|inno|ipaq|iris|ja(t|v)a|jbro|jemu|jigs|kddi|keji|kgt( |\/)|klon|kpt |kwc\-|kyo(c|k)|le(no|xi)|lg( g|\/(k|l|u)|50|54|\-[a-w])|libw|lynx|m1\-w|m3ga|m50\/|ma(te|ui|xo)|mc(01|21|ca)|m\-cr|me(rc|ri)|mi(o8|oa|ts)|mmef|mo(01|02|bi|de|do|t(\-| |o|v)|zz)|mt(50|p1|v )|mwbp|mywa|n10[0-2]|n20[2-3]|n30(0|2)|n50(0|2|5)|n7(0(0|1)|10)|ne((c|m)\-|on|tf|wf|wg|wt)|nok(6|i)|nzph|o2im|op(ti|wv)|oran|owg1|p800|pan(a|d|t)|pdxg|pg(13|\-([1-8]|c))|phil|pire|pl(ay|uc)|pn\-2|po(ck|rt|se)|prox|psio|pt\-g|qa\-a|qc(07|12|21|32|60|\-[2-7]|i\-)|qtek|r380|r600|raks|rim9|ro(ve|zo)|s55\/|sa(ge|ma|mm|ms|ny|va)|sc(01|h\-|oo|p\-)|sdk\/|se(c(\-|0|1)|47|mc|nd|ri)|sgh\-|shar|sie(\-|m)|sk\-0|sl(45|id)|sm(al|ar|b3|it|t5)|so(ft|ny)|sp(01|h\-|v\-|v )|sy(01|mb)|t2(18|50)|t6(00|10|18)|ta(gt|lk)|tcl\-|tdg\-|tel(i|m)|tim\-|t\-mo|to(pl|sh)|ts(70|m\-|m3|m5)|tx\-9|up(\.b|g1|si)|utst|v400|v750|veri|vi(rg|te)|vk(40|5[0-3]|\-v)|vm40|voda|vulc|vx(52|53|60|61|70|80|81|83|85|98)|w3c(\-| )|webc|whit|wi(g |nc|nw)|wmlb|wonu|x700|yas\-|your|zeto|zte\-/i.test(
                a.substr(0, 4),
            )
        )
            check = true;
    })(navigator.userAgent || navigator.vendor); // || window.opera);
    return check;
}
/* eslint-disable no-useless-escape */

/**
 * @param { import('../../src/loaders/omexml').OMEXML[0] } imgMeta
 */
export function guessRgb({ Pixels }: any) {
    const numChannels = Pixels.Channels.length;
    const { SamplesPerPixel } = Pixels.Channels[0];

    const is3Channel8Bit = numChannels === 3 && Pixels.Type === "uint8";
    const interleavedRgb =
        Pixels.SizeC === 3 && numChannels === 1 && Pixels.Interleaved;

    return SamplesPerPixel === 3 || is3Channel8Bit || interleavedRgb;
}
export function truncateDecimalNumber(value: number, maxLength: number) {
    if (!value && value !== 0) return "";
    const stringValue = value.toString();
    return stringValue.length > maxLength
        ? stringValue.substring(0, maxLength).replace(/\.$/, "")
        : stringValue;
}

/**
 * Get physical size scaling Matrix4
 * @param {Object} loader PixelSource - todo: clarify type
 */
export function getPhysicalSizeScalingMatrix(loader: PixelSource | any) {
    const { x, y, z } = loader?.meta?.physicalSizes ?? {};
    if (x?.size && y?.size && z?.size) {
        const min = Math.min(z.size, x.size, y.size);
        const ratio = [x.size / min, y.size / min, z.size / min];
        return new Matrix4().scale(ratio);
    }
    return new Matrix4().identity();
}

export function getBoundingCube(loader: PixelSource) {
    const source = Array.isArray(loader) ? loader[0] : loader;
    const { shape, labels } = source;
    const physicalSizeScalingMatrix = getPhysicalSizeScalingMatrix(source);
    const xSlice: [number, number] = [
        0,
        physicalSizeScalingMatrix[0] * shape[labels.indexOf("x")],
    ];
    const ySlice: [number, number] = [
        0,
        physicalSizeScalingMatrix[5] * shape[labels.indexOf("y")],
    ];
    const zSlice: [number, number] = [
        0,
        physicalSizeScalingMatrix[10] * shape[labels.indexOf("z")],
    ];
    return [xSlice, ySlice, zSlice];
}

export type RenderingMode = RENDERING_MODES[keyof RENDERING_MODES];
/**
 * Return an appropriate 3D extension for a given combination of `colormap` and `renderingMode`
 * @param colormap - supposedly a string but only used as a boolean
 * @param renderingMode
 */
export function get3DExtension(colormap: boolean, renderingMode: RenderingMode) {
    const extensions = colormap
        ? AdditiveColormap3DExtensions
        : ColorPalette3DExtensions;
    if (renderingMode === RENDERING_MODES.MAX_INTENSITY_PROJECTION) {
        return new extensions.MaximumIntensityProjectionExtension(undefined);
    }
    if (renderingMode === RENDERING_MODES.MIN_INTENSITY_PROJECTION) {
        return new extensions.MinimumIntensityProjectionExtension(undefined);
    }
    if (renderingMode === RENDERING_MODES.ADDITIVE) {
        return new extensions.AdditiveBlendExtension(undefined);
    }
    throw new Error(`${renderingMode} rendering mode not supported`);
}

const SI_PREFIXES = [
    { symbol: "Y", exponent: 24 },
    { symbol: "Z", exponent: 21 },
    { symbol: "E", exponent: 18 },
    { symbol: "P", exponent: 15 },
    { symbol: "T", exponent: 12 },
    { symbol: "G", exponent: 9 },
    { symbol: "M", exponent: 6 },
    { symbol: "k", exponent: 3 },
    { symbol: "h", exponent: 2 },
    { symbol: "da", exponent: 1 },
    { symbol: "", exponent: 0 },
    { symbol: "d", exponent: -1 },
    { symbol: "c", exponent: -2 },
    { symbol: "m", exponent: -3 },
    { symbol: "µ", exponent: -6 },
    { symbol: "n", exponent: -9 },
    { symbol: "p", exponent: -12 },
    { symbol: "f", exponent: -15 },
    { symbol: "a", exponent: -18 },
    { symbol: "z", exponent: -21 },
    { symbol: "y", exponent: -24 },
] as const;
type SizeUnit = typeof SI_PREFIXES[number]["symbol"];
/**
 * Convert a size value to meters.
 * @param size Size in original units.
 * @param unit String like 'mm', 'cm', 'dam', 'm', 'km', etc.
 * @returns Size in meters.
 */
export function sizeToMeters(size: number, unit: SizeUnit): number {
    if (!unit || unit === "m") {
        // Already in meters.
        return size;
    }
    if (unit.length > 1) {
        // We remove the trailing 'm' from the unit, so 'cm' becomes 'c' and 'dam' becomes 'da'.
        let unitPrefix = unit.substring(0, unit.length - 1);
        // Support 'u' as a prefix for micrometers.
        if (unitPrefix === "u") {
            unitPrefix = "µ";
        }
        const unitObj = SI_PREFIXES.find((p) => p.symbol === unitPrefix);
        if (unitObj) {
            return size * 10 ** unitObj.exponent;
        }
    }
    throw new Error("Received unknown unit");
}
// export function sizeToNiceUnits(size, unit) {
//     if (unit === 'm') {
//         return size;
//     }
//     if (unit.length > 1) {
//         const unitPrefix = unit.substring(0, unit.length - 1);
//         const unitObj = SI_PREFIXES.find(p => p.symbol === unitPrefix);
//         if (unitObj) {
//             return size / 10 ** unitObj.exponent + unitObj.symbol;
//         }
//     }
//     throw new Error('Received unknown unit');
// }
