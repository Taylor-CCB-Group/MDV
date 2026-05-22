import type { SocketIOUploadClient } from "@/charts/dialogs/SocketIOUploadClient";
import { DEFAULT_UPLOAD_METHOD } from "./fileUploadConfig";
import type {
    ConflictData,
    UploadFileType,
    UploadMethod,
    UploadStage,
    UploadValidationResult,
} from "./fileUploadTypes";

export interface UploadDialogState {
    selectedFiles: File[];
    stage: UploadStage;
    validationResult: UploadValidationResult | null;
    fileType: UploadFileType | null;
    tiffMetadata: unknown;
    h5Metadata: unknown;
    conflictData: ConflictData | null;
    uploadMethod: UploadMethod;
    socketioClient: SocketIOUploadClient | null;
}

export type UploadDialogAction =
    | { type: "SET_SELECTED_FILES"; payload: File[] }
    | { type: "SET_STAGE"; payload: UploadStage }
    | { type: "SET_VALIDATION_RESULT"; payload: UploadValidationResult | null }
    | { type: "SET_FILE_TYPE"; payload: UploadFileType | null }
    | { type: "SET_TIFF_METADATA"; payload: unknown }
    | { type: "SET_H5_METADATA"; payload: unknown }
    | { type: "SET_CONFLICT_DATA"; payload: ConflictData | null }
    | { type: "SET_UPLOAD_METHOD"; payload: UploadMethod }
    | { type: "SET_SOCKETIO_CLIENT"; payload: SocketIOUploadClient | null }
    | { type: "RESET_DIALOG" };

export function createInitialUploadDialogState(): UploadDialogState {
    return {
        selectedFiles: [],
        stage: { kind: "idle" },
        validationResult: null,
        fileType: null,
        tiffMetadata: null,
        h5Metadata: null,
        conflictData: null,
        uploadMethod: DEFAULT_UPLOAD_METHOD,
        socketioClient: null,
    };
}

export function fileUploadReducer(
    state: UploadDialogState,
    action: UploadDialogAction,
): UploadDialogState {
    switch (action.type) {
        case "SET_SELECTED_FILES":
            return { ...state, selectedFiles: action.payload };
        case "SET_STAGE":
            return { ...state, stage: action.payload };
        case "SET_VALIDATION_RESULT":
            return { ...state, validationResult: action.payload };
        case "SET_FILE_TYPE":
            return { ...state, fileType: action.payload };
        case "SET_TIFF_METADATA":
            return { ...state, tiffMetadata: action.payload };
        case "SET_H5_METADATA":
            return { ...state, h5Metadata: action.payload };
        case "SET_CONFLICT_DATA":
            return { ...state, conflictData: action.payload };
        case "SET_UPLOAD_METHOD":
            return { ...state, uploadMethod: action.payload };
        case "SET_SOCKETIO_CLIENT":
            return { ...state, socketioClient: action.payload };
        case "RESET_DIALOG":
            return {
                ...createInitialUploadDialogState(),
                uploadMethod: state.uploadMethod,
            };
        default:
            return state;
    }
}
