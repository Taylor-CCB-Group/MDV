import { describe, expect, test } from "vitest";

import {
    createInitialUploadDialogState,
    fileUploadReducer,
} from "@/utils/file_upload/fileUploadReducer";

describe("fileUploadReducer", () => {
    test("tracks stage transitions explicitly", () => {
        const validatingState = fileUploadReducer(
            createInitialUploadDialogState(),
            { type: "SET_STAGE", payload: { kind: "validating" } },
        );

        const processingState = fileUploadReducer(
            validatingState,
            { type: "SET_STAGE", payload: { kind: "processing" } },
        );

        expect(validatingState.stage).toEqual({ kind: "validating" });
        expect(processingState.stage).toEqual({ kind: "processing" });
    });

    test("stores structured error stage payloads", () => {
        const state = fileUploadReducer(createInitialUploadDialogState(), {
            type: "SET_STAGE",
            payload: {
                kind: "error",
                error: {
                    message: "Upload failed",
                    status: 500,
                    traceback: "boom",
                },
            },
        });

        expect(state.stage).toEqual({
            kind: "error",
            error: {
                message: "Upload failed",
                status: 500,
                traceback: "boom",
            },
        });
    });

    test("reset clears workflow state but preserves chosen upload method", () => {
        let state = createInitialUploadDialogState();
        state = fileUploadReducer(state, {
            type: "SET_UPLOAD_METHOD",
            payload: "http",
        });
        state = fileUploadReducer(state, {
            type: "SET_SELECTED_FILES",
            payload: [new File(["hello"], "data.csv", { type: "text/csv" })],
        });
        state = fileUploadReducer(state, {
            type: "SET_STAGE",
            payload: { kind: "ready" },
        });

        const resetState = fileUploadReducer(state, { type: "RESET_DIALOG" });

        expect(resetState.uploadMethod).toBe("http");
        expect(resetState.stage).toEqual({ kind: "idle" });
        expect(resetState.selectedFiles).toHaveLength(0);
        expect(resetState.validationResult).toBeNull();
    });
});
