import type { UploadFileType, UploadMethod } from "./fileUploadTypes";

export interface FileTypeConfig {
    type: UploadFileType;
    extensions: string[];
    mimeTypes: string[];
    maxSize?: number;
    processingConfig: {
        defaultWidth: number;
        defaultHeight: number;
        requiresMetadata?: boolean;
        endpoint: string;
    };
}

const USE_SOCKETIO_UPLOAD = true;

const FILE_TYPES: Record<string, FileTypeConfig> = {
    CSV: {
        type: "csv",
        extensions: [".csv"],
        mimeTypes: ["text/csv"],
        maxSize: 10000 * 1024 * 1024,
        processingConfig: {
            defaultWidth: 800,
            defaultHeight: 745,
            requiresMetadata: false,
            endpoint: "add_datasource",
        },
    },
    TSV: {
        type: "tsv",
        extensions: [".tsv", ".tab", ".txt"],
        mimeTypes: ["text/tab-separated-values", "text/plain"],
        maxSize: 10000 * 1024 * 1024,
        processingConfig: {
            defaultWidth: 800,
            defaultHeight: 745,
            requiresMetadata: false,
            endpoint: "add_datasource",
        },
    },
    TIFF: {
        type: "tiff",
        extensions: [".tiff", ".tif"],
        mimeTypes: ["image/tiff"],
        maxSize: 10000 * 1024 * 1024,
        processingConfig: {
            defaultWidth: 1032,
            defaultHeight: 580,
            requiresMetadata: true,
            endpoint: "add_or_update_image_datasource",
        },
    },
    H5: {
        type: "h5",
        extensions: [".h5", ".h5ad"],
        mimeTypes: ["application/x-hdf5", "application/x-hdf"],
        maxSize: 10000 * 1024 * 1024,
        processingConfig: {
            defaultWidth: 1000,
            defaultHeight: 800,
            requiresMetadata: true,
            endpoint: "add_anndata",
        },
    },
};

const IS_PRODUCTION = import.meta.env.PROD;

export const ALLOWED_FILE_TYPES: Record<string, FileTypeConfig> = IS_PRODUCTION
    ? {
          CSV: FILE_TYPES.CSV,
          TSV: FILE_TYPES.TSV,
      }
    : FILE_TYPES;

export const DEFAULT_UPLOAD_METHOD: UploadMethod = USE_SOCKETIO_UPLOAD
    ? "socketio"
    : "http";

export function getFileTypeFromExtension(fileName: string): FileTypeConfig | null {
    const extension = `.${fileName.split(".").pop()?.toLowerCase()}`;
    return (
        Object.values(ALLOWED_FILE_TYPES).find((config) =>
            config.extensions.includes(extension),
        ) || null
    );
}

export function generateDropzoneAccept(): Record<string, string[]> {
    return Object.values(ALLOWED_FILE_TYPES).reduce(
        (acc, config) => {
            config.mimeTypes.forEach((mimeType) => {
                acc[mimeType] = config.extensions;
            });
            return acc;
        },
        {} as Record<string, string[]>,
    );
}

export function getSupportedFileExtensions(): string[] {
    return Array.from(
        new Set(
            Object.values(ALLOWED_FILE_TYPES).flatMap(
                (config) => config.extensions,
            ),
        ),
    );
}
