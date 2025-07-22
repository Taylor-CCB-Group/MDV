/// <reference types="emscripten" />
// see https://github.com/cornerstonejs/cornerstone3D/blob/014f4c4cc2b973b200ec9af2e16783464b9a2a0d/packages/dicomImageLoader/src/types/codec-openjpeg.d.ts
declare module "@cornerstonejs/codec-openjpeg/dist/openjpegwasm_decode" {
    export class J2KDecoder {
        decode: () => any;
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
