/// <reference types="emscripten" />
// see https://github.com/cornerstonejs/cornerstone3D/blob/014f4c4cc2b973b200ec9af2e16783464b9a2a0d/packages/dicomImageLoader/src/types/codec-openjpeg.d.ts
// & the emscripten source

/*
```cpp

EMSCRIPTEN_BINDINGS(J2KDecoder) {
  class_<J2KDecoder>("J2KDecoder")
    .constructor<>()
    .function("getEncodedBuffer", &J2KDecoder::getEncodedBuffer)
    .function("getDecodedBuffer", &J2KDecoder::getDecodedBuffer)
    .function("readHeader", &J2KDecoder::readHeader)
    .function("calculateSizeAtDecompositionLevel", &J2KDecoder::calculateSizeAtDecompositionLevel)
    .function("decode", &J2KDecoder::decode)
    .function("decodeSubResolution", &J2KDecoder::decodeSubResolution)
    .function("getFrameInfo", &J2KDecoder::getFrameInfo)
    .function("getNumDecompositions", &J2KDecoder::getNumDecompositions)
    .function("getIsReversible", &J2KDecoder::getIsReversible)
    .function("getProgressionOrder", &J2KDecoder::getProgressionOrder)
    .function("getImageOffset", &J2KDecoder::getImageOffset)
    .function("getTileSize", &J2KDecoder::getTileSize)
    .function("getTileOffset", &J2KDecoder::getTileOffset)
    .function("getBlockDimensions", &J2KDecoder::getBlockDimensions)
    .function("getNumLayers", &J2KDecoder::getNumLayers)
    .function("getColorSpace", &J2KDecoder::getColorSpace)
   ;
}
```
 */


declare module "@cornerstonejs/codec-openjpeg/dist/openjpegwasm_decode" {
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
