import { BaseDecoder, type TypedArray, addDecoder } from "geotiff";

import openJpegFactory from "@cornerstonejs/codec-openjpeg/decodewasmjs";
// - nb, comments from https://github.com/cornerstonejs/cornerstone3D/blob/014f4c4cc2b973b200ec9af2e16783464b9a2a0d/packages/dicomImageLoader/src/shared/decoders/decodeJPEG2000.ts#L5
// Webpack asset/resource copies this to our output folder

// TODO: At some point maybe we can use this instead.
// This is closer to what Webpack 5 wants but it doesn't seem to work now
// const wasm = new URL('./blah.wasm', import.meta.url)
// import openjpegWasm from '@cornerstonejs/codec-openjpeg/decodewasm';
const openjpegWasm = new URL("@cornerstonejs/codec-openjpeg/decodewasm", import.meta.url);

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

    constructor(fileDirectory: FileDirectory) {
        super();
    }

    private async getOpenJPEG(): Promise<OpenJpegModule> {
        if (!this.openjpeg) {
            console.log(">>> Initializing OpenJPEG");
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

    // async decode(fileDirectory: FileDirectory, buffer: TypedArray) {
    //     const compressedImageFrame = await this.decodeBlock(buffer);
    //     return compressedImageFrame;
    // }
    
    async decodeBlock(compressedImageFrame: TypedArray) {
        // this is more-or-less a copy of the code in conerstone3D's decodeJPEG2000.ts
        try {
            const openjpeg = await this.getOpenJPEG();
            // console.log("Inspecting openjpeg object:", Object.keys(openjpeg));
            // the thing we have as openjpeg here doesn't have a decode method...
            // but it might have J2KDecoder...
            const decoder = new openjpeg.J2KDecoder();
            const encodedBufferInWASM = decoder.getEncodedBuffer(compressedImageFrame.length);
            encodedBufferInWASM.set(compressedImageFrame);
            // this won't throw - but it will complain...
            // openjpegwasm_decode.js:9 [INFO] Stream reached its end !
            // openjpegwasm_decode.js:9 [ERROR] JP2H box missing. Required.
            // openjpegwasm_decode.js: 9[ERROR] opj_decompress: failed to read the header
            decoder.decode();
            // get information about the decoded image
            const frameInfo = decoder.getFrameInfo();
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
function getPixelData(frameInfo: any, decodedBuffer: any) {
    if (frameInfo.bitsPerSample > 8) {
        if (frameInfo.isSigned) {
            return new Int16Array(decodedBuffer.buffer, decodedBuffer.byteOffset, decodedBuffer.byteLength / 2);
        }

        return new Uint16Array(decodedBuffer.buffer, decodedBuffer.byteOffset, decodedBuffer.byteLength / 2);
    }

    if (frameInfo.isSigned) {
        return new Int8Array(decodedBuffer.buffer, decodedBuffer.byteOffset, decodedBuffer.byteLength);
    }

    return new Uint8Array(decodedBuffer.buffer, decodedBuffer.byteOffset, decodedBuffer.byteLength);
}

// Register the decoder with geotiff
// JPEG2000 compression code in TIFF
const JPEG2000_COMPRESSION = 34712; // this is what we observe `getDecoder` using internally (slightly sidetracked by another image having 8)
try {
    addDecoder(JPEG2000_COMPRESSION, () => {
        console.log(">>> Using JPEG2000 decoder");
        return Promise.resolve(Jpeg2000Decoder);
    });
    console.log("JPEG2000 decoder registered successfully");
} catch (error) {
    console.error("Failed to register JPEG2000 decoder:", error);
}
