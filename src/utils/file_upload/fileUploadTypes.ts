export type UploadMethod = "http" | "socketio";

export type UploadFileType = "csv" | "tiff" | "tsv" | "h5";

export interface UploadErrorPayload {
    message: string;
    status: number;
    traceback?: string;
}

export interface UploadValidationResult {
    columnNames: string[];
    columnTypes: string[];
}

export interface ConflictData {
    temp_folder: string;
}

export type UploadStage =
    | { kind: "idle" }
    | { kind: "validating" }
    | { kind: "ready" }
    | { kind: "uploading" }
    | { kind: "processing" }
    | { kind: "conflict" }
    | { kind: "redirecting"; viewName: string }
    | { kind: "success" }
    | { kind: "error"; error: UploadErrorPayload };

export interface DatasourceSummary {
    datasourceName: string;
    fileName: string;
    fileSize: string;
    rowCount: number;
    columnCount: number;
}

export interface TiffSummary {
    fileName: string;
    fileSize: string;
}
