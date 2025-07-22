/// <reference types="emscripten" />

declare module "@cornerstonejs/codec-openjpeg/decodewasmjs" {
    export class J2KDecoder {
        decode: () => any;
        readHeader: () => any;
        calculateSizeAtDecompositionLevel: (level: number) => any;
        /**
         * Decodes the encoded HTJ2K[sic] bitstream to the requested decomposition level.
         * The caller must have copied the HTJ2K[sic] encoded bitstream into the encoded 
         * buffer before calling this method, see getEncodedBuffer() and
         * getEncodedBytes() above.
         */
        decodeSubResolution: (resolution: number, layer: number) => void;
        getBlockDimensions: () => any;
        getColorSpace: () => any;
        getDecodedBuffer: () => any;
        getEncodedBuffer: (length: number) => any;
        getFrameInfo: () => any;
        getImageOffset: () => any;
        getIsReversible: () => any;
        getNumDecompositions: () => any;
        getNumLayers: () => any;
        getProgressionOrder: () => number;
        getTileOffset: () => any;
        getTileSize: () => any;
    }
    export interface OpenJpegModule extends EmscriptenModule {
        J2KDecoder: typeof J2KDecoder;
    }
    declare const Module: EmscriptenModuleFactory<OpenJpegModule>;
    export default Module;
} 