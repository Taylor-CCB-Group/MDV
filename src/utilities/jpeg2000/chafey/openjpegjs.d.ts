// TypeScript bindings for emscripten-generated code.  Automatically generated at compile time.
declare namespace RuntimeExports {
    /**
     * @param {string|null=} returnType
     * @param {Array=} argTypes
     * @param {Arguments|Array=} args
     * @param {Object=} opts
     */
    function ccall(ident: any, returnType?: (string | null) | undefined, argTypes?: any[] | undefined, args?: (Arguments | any[]) | undefined, opts?: any | undefined): any;
    let HEAPF32: any;
    let HEAPF64: any;
    let HEAP_DATA_VIEW: any;
    let HEAP8: any;
    let HEAPU8: any;
    let HEAP16: any;
    let HEAPU16: any;
    let HEAP32: any;
    let HEAPU32: any;
    let HEAP64: any;
    let HEAPU64: any;
}
interface WasmModule {
}

export interface ClassHandle {
  isAliasOf(other: ClassHandle): boolean;
  delete(): void;
  deleteLater(): this;
  isDeleted(): boolean;
  clone(): this;
}
export type FrameInfo = {
  width: number,
  height: number,
  bitsPerSample: number,
  componentCount: number,
  isSigned: boolean
};

export type Point = {
  x: number,
  y: number
};

export type Size = {
  width: number,
  height: number
};

export interface J2KDecoder extends ClassHandle {
  getEncodedBuffer(_0: number): any;
  getDecodedBuffer(): any;
  readHeader(): void;
  calculateSizeAtDecompositionLevel(_0: number): Size;
  decode(): void;
  decodeSubResolution(_0: number, _1: number): void;
  getFrameInfo(): FrameInfo;
  getNumDecompositions(): number;
  getIsReversible(): boolean;
  getProgressionOrder(): number;
  getImageOffset(): Point;
  getTileSize(): Size;
  getTileOffset(): Point;
  getBlockDimensions(): Size;
  getNumLayers(): number;
  getColorSpace(): number;
}

export interface J2KEncoder extends ClassHandle {
  getDecodedBuffer(_0: FrameInfo): any;
  getEncodedBuffer(): any;
  encode(): void;
  setDecompositions(_0: number): void;
  setQuality(_0: boolean, _1: number): void;
  setProgressionOrder(_0: number): void;
  setDownSample(_0: number, _1: Point): void;
  setImageOffset(_0: Point): void;
  setTileSize(_0: Size): void;
  setTileOffset(_0: Point): void;
  setBlockDimensions(_0: Size): void;
  setNumPrecincts(_0: number): void;
  setPrecinct(_0: number, _1: Size): void;
  setCompressionRatio(_0: number, _1: number): void;
}

interface EmbindModule {
  getVersion(): string;
  J2KDecoder: {
    new(): J2KDecoder;
  };
  J2KEncoder: {
    new(): J2KEncoder;
  };
}

export type MainModule = WasmModule & typeof RuntimeExports & EmbindModule;
export default function MainModuleFactory (options?: unknown): Promise<MainModule>;
