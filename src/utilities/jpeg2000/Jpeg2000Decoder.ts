import { BaseDecoder, type TypedArray, addDecoder } from "geotiff";
import type { FrameInfo, MainModule as OpenJpegModule } from "./chafey/openjpegjs";
import openJpegFactory from "./chafey/openjpegjs.js";

// - nb, some comments from https://github.com/cornerstonejs/cornerstone3D/blob/014f4c4cc2b973b200ec9af2e16783464b9a2a0d/packages/dicomImageLoader/src/shared/decoders/decodeJPEG2000.ts#L5

// todo - publish a clean version of this to npm
const openjpegWasm = new URL("./chafey/openjpegjs.wasm", import.meta.url);

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

*/

export default class Jpeg2000Decoder extends BaseDecoder {
    private static openjpeg?: Promise<OpenJpegModule>;

    private async getOpenJPEG(): Promise<OpenJpegModule> {
        if (!Jpeg2000Decoder.openjpeg) {
            try {
                Jpeg2000Decoder.openjpeg = openJpegFactory({
                    locateFile: (file: string) => {
                        if (file.endsWith(".wasm")) {
                            return openjpegWasm.href;
                        } else {
                            return file;
                        }
                    },
                });
            } catch (error) {
                console.warn("WASM version failed, JS version fallback not attempted:", error);
                throw new Error("Failed to initialize OpenJPEG codec");
            }
        }
        return Jpeg2000Decoder.openjpeg;
    }

    async decodeBlock(compressedImageFrame: ArrayBuffer) {
        // this is more-or-less a copy of the code in conerstone3D's decodeJPEG2000.ts
        // (or it was originally, now various changes to make it do what we need)
        let decoder: InstanceType<OpenJpegModule["J2KDecoder"]> | null = null;
        try {
            const openjpeg = await this.getOpenJPEG();
            decoder = new openjpeg.J2KDecoder();
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
            // get the decoded pixels
            const decodedBufferInWASM = decoder.getDecodedBuffer();
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
        } finally {
            try {
                decoder?.delete();
            } catch (error) {}
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
const JPEG2000_COMPRESSION = [34712, 33004, 33003]; // this is what we observe `getDecoder` using internally (slightly sidetracked by another image having 8)
try {
    addDecoder(JPEG2000_COMPRESSION, () => {
        return Promise.resolve(Jpeg2000Decoder);
    });
} catch (error) {
    console.error("Failed to register JPEG2000 decoder:", error);
}
