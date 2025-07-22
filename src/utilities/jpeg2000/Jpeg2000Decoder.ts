import { BaseDecoder, type TypedArray, addDecoder } from "geotiff";

// @ts-ignore
import openJpegFactory from "@cornerstonejs/codec-openjpeg/decodewasmjs";
// - nb, comments from https://github.com/cornerstonejs/cornerstone3D/blob/014f4c4cc2b973b200ec9af2e16783464b9a2a0d/packages/dicomImageLoader/src/shared/decoders/decodeJPEG2000.ts#L5
// Webpack asset/resource copies this to our output folder

// TODO: At some point maybe we can use this instead.
// This is closer to what Webpack 5 wants but it doesn't seem to work now
// const wasm = new URL('./blah.wasm', import.meta.url)
// @ts-ignore
// import openjpegWasm from '@cornerstonejs/codec-openjpeg/decodewasm';
const openjpegWasm = new URL("@cornerstonejs/codec-openjpeg/decodewasm", import.meta.url);

// passed as argument to BaseDecoder constructor
// interface FileDirectory {
//     TileWidth?: number;
//     TileLength?: number;
//     ImageWidth: number;
//     ImageLength: number;
//     BitsPerSample: number[];
//     Compression: number;
// }

type CompressedImageFrame = {
    length: number;
    // offset: number;
    // tileWidth: number;
    // tileLength: number;
    // imageWidth: number;
    // imageLength: number;
};

// JPEG2000 compression code in TIFF
const JPEG2000_COMPRESSION = 34712; // 0x8790

export default class Jpeg2000Decoder extends BaseDecoder {
    private openjpeg: any;

    private async getOpenJPEG() {
        if (!this.openjpeg) {
            try {
                // Try WASM version first, fall back to JS version
                this.openjpeg = await openJpegFactory({
                    locateFile: (file: string) => {
                        if (file.endsWith(".wasm")) {
                            return openjpegWasm.toString();
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
        return this.openjpeg;
    }

    async decodeBlock(compressedImageFrame: TypedArray) {
        // this is more-or-less a copy of the code in conerstone3D's decodeJPEG2000.ts
        try {
            const openjpeg = await this.getOpenJPEG();
            // the thing we have as openjpeg here doesn't have a decode method...
            // but it might have J2KDecoder...
            const decoder = new openjpeg.J2KDecoder();
            const encodedBufferInWASM = decoder.getEncodedBuffer(compressedImageFrame.length);
            encodedBufferInWASM.set(compressedImageFrame);
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
try {
    addDecoder(JPEG2000_COMPRESSION, () => Promise.resolve(Jpeg2000Decoder));
    console.log("JPEG2000 decoder registered successfully");
} catch (error) {
    console.error("Failed to register JPEG2000 decoder:", error);
}
