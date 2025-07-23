import { BaseDecoder, type TypedArray, addDecoder } from "geotiff";

import openJpegFactory from "./chafey/openjpegjs.js";
// - nb, comments from https://github.com/cornerstonejs/cornerstone3D/blob/014f4c4cc2b973b200ec9af2e16783464b9a2a0d/packages/dicomImageLoader/src/shared/decoders/decodeJPEG2000.ts#L5
// Webpack asset/resource copies this to our output folder

// TODO: At some point maybe we can use this instead.
// This is closer to what Webpack 5 wants but it doesn't seem to work now
// const wasm = new URL('./blah.wasm', import.meta.url)
// import openjpegWasm from '@cornerstonejs/codec-openjpeg/decodewasm';
const openjpegWasm = new URL("./chafey/openjpegjs.wasm", import.meta.url);

// --- START OF TYPES ---
// These are based on the emscripten bindings and should ideally live in a .d.ts file
// but for now, we'll define them here to ensure they are applied.
interface Size {
    width: number;
    height: number;
}
interface Point {
    x: number;
    y: number;
}

//  OPJ_CLRSPC_UNKNOWN = -1,    /**< not supported by the library */
//  OPJ_CLRSPC_UNSPECIFIED = 0, /**< not specified in the codestream */
//  OPJ_CLRSPC_SRGB = 1,        /**< sRGB */
//  OPJ_CLRSPC_GRAY = 2,        /**< grayscale */
//  OPJ_CLRSPC_SYCC = 3,        /**< YUV */
//  OPJ_CLRSPC_EYCC = 4,        /**< e-YCC */
//  OPJ_CLRSPC_CMYK = 5         /**< CMYK */
type ColorSpace = -1 | 0 | 1 | 2 | 3 | 4 | 5;

interface FrameInfo {
    width: number;
    height: number;
    bitsPerSample: number;
    componentCount: number;
    isSigned: boolean;
}

interface J2KDecoder {
    /**
     * Decodes the encoded HTJ2K bitstream.  The caller must have copied the
     * HTJ2K encoded bitstream into the encoded buffer before calling this
     * method, see getEncodedBuffer() and getEncodedBytes() (latter not exported to js)
     */
    decode: () => void;
    readHeader: () => void;
    calculateSizeAtDecompositionLevel: (level: number) => Size;
    decodeSubResolution: (resolution: number, layer: number) => void;
    getBlockDimensions: () => Size;
    getColorSpace: () => ColorSpace;
    getDecodedBuffer: () => TypedArray;
    getEncodedBuffer: (length: number) => TypedArray;
    getFrameInfo: () => FrameInfo;
    getImageOffset: () => Point;
    getIsReversible: () => boolean;
    getNumDecompositions: () => number;
    getNumLayers: () => number;
    getProgressionOrder: () => number;
    getTileOffset: () => Point;
    getTileSize: () => Size;
}

interface OpenJpegModule {
    J2KDecoder: {
        new(): J2KDecoder;
    };
}

// --- END OF TYPES ---


// passed as argument to BaseDecoder constructor / decode
/* e.g.
{
    "NewSubfileType": 1,
    "ImageWidth": 361,
    "ImageLength": 220,
    "BitsPerSample": {
        "0": 16
    },
    "Compression": 34712,
    "PhotometricInterpretation": 1,
    "SamplesPerPixel": 1,
    "XResolution": {
        "0": 1,
        "1": 1
    },
    "YResolution": {
        "0": 1,
        "1": 1
    },
    "ResolutionUnit": 1,
    "TileWidth": 1024,
    "TileLength": 1024,
    "TileOffsets": [
        967136
    ],
    "TileByteCounts": [
        2193
    ]
}
*/

interface FileDirectory {
    TileWidth?: number;
    TileLength?: number;
    ImageWidth: number;
    ImageLength: number;
    BitsPerSample: Record<string, number>;
    Compression: number;
}

type CompressedImageFrame = {
    length: number;
    // offset: number;
    // tileWidth: number;
    // tileLength: number;
    // imageWidth: number;
    // imageLength: number;
};


export default class Jpeg2000Decoder extends BaseDecoder {
    private openjpeg: OpenJpegModule | null = null;

    // constructor(fileDirectory: FileDirectory) {
    //     super();
    // }

    private async getOpenJPEG(): Promise<OpenJpegModule> {
        if (!this.openjpeg) {
            try {
                // Try WASM version first, fall back to JS version
                this.openjpeg = (await openJpegFactory({
                    locateFile: (file: string) => {
                        if (file.endsWith(".wasm")) {
                            return openjpegWasm.href;
                        } else {
                            return file;
                        }
                    },
                })) as any as OpenJpegModule;
            } catch (error) {
                console.warn("WASM version failed, JS version fallback not attempted:", error);
                throw new Error("Failed to initialize OpenJPEG codec");
            }
        }
        return this.openjpeg;
    }

    async decodeBlock(compressedImageFrame: ArrayBuffer) {
        // this is more-or-less a copy of the code in conerstone3D's decodeJPEG2000.ts
        try {
            const openjpeg = await this.getOpenJPEG();
            const decoder = new openjpeg.J2KDecoder();
            const buffer = new Uint8Array(compressedImageFrame);
            const encodedBufferInWASM = decoder.getEncodedBuffer(buffer.length);
            encodedBufferInWASM.set(buffer);
            decoder.decode();
            // get information about the decoded image
            const frameInfo = decoder.getFrameInfo();
            // we could probably look at the frameInfo to see if the decode method worked properly
            // it won't throw an error otherwise...
            if (frameInfo.width === 0 || frameInfo.height === 0) {
                throw new Error("Failed to decode JPEG2000 image");
            }
            // console.log("Frame info:", frameInfo);
            // get the decoded pixels
            const decodedBufferInWASM = decoder.getDecodedBuffer();
            const imageFrame = new Uint8Array(decodedBufferInWASM.length);
            imageFrame.set(decodedBufferInWASM);
            const pixelData = getPixelData(frameInfo, decodedBufferInWASM);

            return pixelData;
        } catch (error) {
            console.error("JPEG2000 decoding failed:", error);

            // If the error is related to invalid field types, provide a more helpful message
            if (error instanceof RangeError && error.message.includes("Invalid field type")) {
                throw new Error(
                    `TIFF file appears to be corrupted with invalid field type. This may be due to file corruption or an unsupported TIFF variant. Original error: ${error.message}`,
                );
            }

            throw error;
        }
    }
}
function getPixelData(frameInfo: FrameInfo, decodedBuffer: TypedArray) {
    let pixelData: TypedArray;
    if (frameInfo.bitsPerSample > 8) {
        pixelData = frameInfo.isSigned
            ? new Int16Array(decodedBuffer.buffer, decodedBuffer.byteOffset, decodedBuffer.byteLength / 2)
            : new Uint16Array(decodedBuffer.buffer, decodedBuffer.byteOffset, decodedBuffer.byteLength / 2);
    } else {
        pixelData = frameInfo.isSigned
            ? new Int8Array(decodedBuffer.buffer, decodedBuffer.byteOffset, decodedBuffer.byteLength)
            : new Uint8Array(decodedBuffer.buffer, decodedBuffer.byteOffset, decodedBuffer.byteLength);
    }

    // Return a buffer containing only the pixel data
    return pixelData.buffer.slice(pixelData.byteOffset, pixelData.byteOffset + pixelData.byteLength);
}

// Register the decoder with geotiff
// JPEG2000 compression code in TIFF
const JPEG2000_COMPRESSION = 34712; // this is what we observe `getDecoder` using internally (slightly sidetracked by another image having 8)
try {
    addDecoder(JPEG2000_COMPRESSION, () => {
        return Promise.resolve(Jpeg2000Decoder);
    });
} catch (error) {
    console.error("Failed to register JPEG2000 decoder:", error);
}
