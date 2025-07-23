export interface Size {
    width: number;
    height: number;
}
export interface Point {
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

export interface FrameInfo {
    width: number;
    height: number;
    bitsPerSample: number;
    componentCount: number;
    isSigned: boolean;
}

export interface J2KDecoder {
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

export interface OpenJpegModule {
    J2KDecoder: {
        new(): J2KDecoder;
    };
}