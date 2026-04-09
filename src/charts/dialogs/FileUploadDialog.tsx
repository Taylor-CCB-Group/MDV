import type React from "react";
import {
    useState,
    useCallback,
    useReducer,
    useEffect,
    useMemo
} from "react";
import { useDropzone } from "react-dropzone";
import { observer } from "mobx-react-lite";

import axios, { type AxiosError, type AxiosProgressEvent } from "axios";
import { useProject } from "../../modules/ProjectContext";
import {
    useViewerStoreApi,
    useChannelsStoreApi,
} from "../../react/components/avivatorish/state";

import AnndataConflictDialog from "./AnndataConflictDialog";
import ReusableDialog from "./ReusableDialog";

import { createSocketIOUpload, type SocketIOUploadClient } from "./SocketIOUploadClient";
import {
    Container,
    DropzoneContainer,
    DynamicText,
    FileInputLabel,
} from "./file_upload/FileUploadUI";
import {
    FileDropzoneView,
    H5ReviewView,
    ProcessingView,
    RedirectingView,
    TabularReviewView,
    TiffReviewView,
    UploadErrorView,
    UploadingView,
    UploadSuccessView,
    ValidatingView,
} from "./file_upload/FileUploadViews";
import {
    ALLOWED_FILE_TYPES,
    generateDropzoneAccept,
    getFileTypeFromExtension,
    getSupportedFileExtensions,
} from "@/utils/file_upload/fileUploadConfig";
import {
    buildViewRedirectUrl,
    formatFileSizeMb,
    getProjectIdFromRoot,
    getSocketPath,
} from "@/utils/file_upload/fileUploadHelpers";
import {
    createInitialUploadDialogState,
    fileUploadReducer,
} from "@/utils/file_upload/fileUploadReducer";
import type {
    DatasourceSummary,
    TiffSummary,
    UploadErrorPayload,
    UploadValidationResult,
} from "@/utils/file_upload/fileUploadTypes";

export {
    Container,
    DropzoneContainer,
    DynamicText,
    FileInputLabel,
};

// Use dynamic import for the worker
const DatasourceWorker = new Worker(
    new URL("./datasourceWorker.ts", import.meta.url),
    {
        type: "module",
    },
);

const EMPTY_VALIDATION_RESULT: UploadValidationResult = {
    columnNames: [],
    columnTypes: [],
};

const EMPTY_DATASOURCE_SUMMARY: DatasourceSummary = {
    datasourceName: "",
    fileName: "",
    fileSize: "",
    rowCount: 0,
    columnCount: 0,
};

const EMPTY_TIFF_SUMMARY: TiffSummary = {
    fileName: "",
    fileSize: "",
};

export interface FileUploadDialogComponentProps {
    onClose: () => void;
    onResize: (width: number, height: number) => void;
    onLoadingStateChange: (isLoading: boolean) => void;
    socketioClientRef?: (client: SocketIOUploadClient | null) => void;
}

// Custom hook for file upload progress
const useFileUploadProgress = () => {
    const [progress, setProgress] = useState(0);

    const resetProgress = useCallback(() => {
        setProgress(0);
    }, []);

    return { progress, setProgress, resetProgress };
};

//todo figure out why type inference is failing with observer if we take out the `React.FC<...>` here
//also - are we actually observing anything here? not that it particularly matters (apart from type inference, apparently)
const FileUploadDialogComponent: React.FC<FileUploadDialogComponentProps> =
    observer(({ onClose, onResize, onLoadingStateChange, socketioClientRef }: FileUploadDialogComponentProps) => {
        const { root, chartManager, mainApiRoute } = useProject();

        const [selectedOption, setSelectedOption] = useState<string | null>(
            null,
        );

        const [updatedNamesArray, setUpdatedNamesArray] = useState<string[]>(
            [],
        );

        const handleSelect = useCallback((value: string) => {
            setSelectedOption(value);
            setDatasourceName(value);
            console.log("Selected:", value);
        }, []);

        const [datasourceName, setDatasourceName] = useState("");

        const [showMetadata, setShowMetadata] = useState(true);

        const toggleView = () => {
            setShowMetadata((prevState) => !prevState);
        };

        const [datasourceSummary, setDatasourceSummary] = useState<DatasourceSummary>(
            EMPTY_DATASOURCE_SUMMARY,
        );

        const [tiffSummary, setTiffSummary] = useState<TiffSummary>(
            EMPTY_TIFF_SUMMARY,
        );

        const [columnNames, setColumnNames] = useState<string[]>([]);
        const [columnTypes, setColumnTypes] = useState<string[]>([]);
        const [secondRowValues, setSecondRowValues] = useState<string[]>([]);

        // TIFF
        const channelsStore = useChannelsStoreApi();

        const [state, dispatch] = useReducer(
            fileUploadReducer,
            undefined,
            createInitialUploadDialogState,
        );

        const isBusy =
            state.stage.kind === "uploading" ||
            state.stage.kind === "processing";

        useEffect(() => {
            onLoadingStateChange(isBusy);
        }, [isBusy, onLoadingStateChange]);

        useEffect(() => {
            if (socketioClientRef) {
                socketioClientRef(state.socketioClient);
            }
        }, [state.socketioClient, socketioClientRef]);

        const viewerStore = useViewerStoreApi();

        const { progress, resetProgress, setProgress } =
            useFileUploadProgress();

        const resetLocalState = useCallback(() => {
            setDatasourceName("");
            setShowMetadata(true);
            setDatasourceSummary(EMPTY_DATASOURCE_SUMMARY);
            setTiffSummary(EMPTY_TIFF_SUMMARY);
            resetProgress();
        }, [resetProgress]);

        const setErrorStage = useCallback((error: UploadErrorPayload) => {
            dispatch({
                type: "SET_STAGE",
                payload: { kind: "error", error },
            });
        }, []);

        const onDrop = useCallback(
            async (acceptedFiles: File[]) => {
                if (acceptedFiles.length === 0) {
                    return;
                }

                const file = acceptedFiles[0];
                const fileConfig = getFileTypeFromExtension(file.name);

                if (!fileConfig) {
                    setErrorStage({
                        message: "Unsupported file type",
                        status: 400,
                    });
                    return;
                }

                dispatch({ type: "RESET_DIALOG" });
                resetLocalState();
                dispatch({
                    type: "SET_SELECTED_FILES",
                    payload: acceptedFiles,
                });
                dispatch({
                    type: "SET_FILE_TYPE",
                    payload: fileConfig.type,
                });
                dispatch({
                    type: "SET_STAGE",
                    payload: { kind: "validating" },
                });

                switch (fileConfig.type) {
                    case "tsv":
                    case "csv": {
                        const newDatasourceName = file.name;
                        setDatasourceName(newDatasourceName);
                        setDatasourceSummary({
                            datasourceName: newDatasourceName,
                            fileName: file.name,
                            fileSize: formatFileSizeMb(file.size),
                            rowCount: 0,
                            columnCount: 0,
                        });
                        // // Process CSV with Web Worker
                        // DatasourceWorker.postMessage({
                        //     file: file,
                        //     fileType: fileConfig.type,
                        // });
                        // DatasourceWorker.onmessage = (
                        //     event: MessageEvent,
                        // ) => {
                        //     const {
                        //         columnNames,
                        //         columnTypes,
                        //         secondRowValues,
                        //         rowCount,
                        //         columnCount,
                        //         error,
                        //     } = event.data;

                        //     if (error) {
                        //         setErrorStage({
                        //             message: "Validation failed.",
                        //             status: 400,
                        //             traceback: error,
                        //         });
                        //     } else {
                        //         setColumnNames(columnNames);
                        //         setColumnTypes(columnTypes);
                        //         setSecondRowValues(secondRowValues);
                        //         setDatasourceSummary(
                        //             (prevDatasrouceSummary) => ({
                        //                 ...prevDatasrouceSummary,
                        //                 rowCount,
                        //                 columnCount,
                        //             }),
                        //         );

                        //         const totalWidth = calculateTotalWidth(
                        //             columnNames,
                        //             columnTypes,
                        //             secondRowValues,
                        //         );
                        //         onResize(totalWidth, 745);
                        //     }
                        // };
                        dispatch({
                            type: "SET_VALIDATION_RESULT",
                            payload: EMPTY_VALIDATION_RESULT,
                        });
                        dispatch({
                            type: "SET_STAGE",
                            payload: { kind: "ready" },
                        });
                        break;
                    }
                    case "tiff": {
                        const dataSources =
                            window.mdv.chartManager?.dataSources ?? [];
                        const namesArray = dataSources.map(
                            (dataSource) => dataSource.name,
                        );
                        setUpdatedNamesArray([
                            ...namesArray,
                            "new datasource",
                        ]);
                        setTiffSummary({
                            fileName: file.name,
                            fileSize: formatFileSizeMb(file.size),
                        });
                        // Process TIFF file
                        // handleSubmitFile(acceptedFiles);
                        onResize(
                            fileConfig.processingConfig.defaultWidth,
                            fileConfig.processingConfig.defaultHeight,
                        );
                        dispatch({
                            type: "SET_VALIDATION_RESULT",
                            payload: EMPTY_VALIDATION_RESULT,
                        });
                        dispatch({
                            type: "SET_STAGE",
                            payload: { kind: "ready" },
                        });
                        break;
                    }
                    case "h5": {
                        try {
                            setDatasourceName(file.name);
                            onResize(
                                fileConfig.processingConfig.defaultWidth,
                                fileConfig.processingConfig.defaultHeight,
                            );
                            dispatch({
                                type: "SET_VALIDATION_RESULT",
                                payload: EMPTY_VALIDATION_RESULT,
                            });
                            dispatch({
                                type: "SET_STAGE",
                                payload: { kind: "ready" },
                            });
                        } catch (error) {
                            console.error("H5 processing error:", error);
                            setErrorStage({
                                message: "Error processing H5 file",
                                status: 500,
                                traceback:
                                    error instanceof Error
                                        ? error.message
                                        : String(error),
                            });
                        }
                        break;
                    }
                    default:
                        console.warn("Unhandled file type:", fileConfig.type);
                        return;
                }
            },
            [onResize, resetLocalState, setErrorStage],
        );

        const handleDatasourceNameChange = useCallback((event: any) => {
            const { value } = event.target;
            setDatasourceName(value); // Update the state with the new value
            setDatasourceSummary((prevDatasrouceSummary) => ({
                ...prevDatasrouceSummary,
                datasourceName: value, // Update datasourceName in the summary object
            }));
        }, []);

        const { getRootProps, getInputProps, isDragActive, fileRejections } =
            useDropzone({
                onDrop,
                accept: generateDropzoneAccept(),
                multiple: false,
                maxSize: Math.max(
                    ...Object.values(ALLOWED_FILE_TYPES).map(
                        (config) => config.maxSize || 0,
                    ),
                ),
            });

        const rejectionMessage = useMemo(() => {
            let message = "";
            if (fileRejections.length > 0) {
                const fileRejectionMessages: string[] = [];
                fileRejections.forEach((fileRejection) => {
                    fileRejection.errors.forEach((error) => {
                        fileRejectionMessages.push(error.message);
                    })
                })

                message = fileRejectionMessages.join(", ");
            } else {
                message = "Drag and drop files here or click the button below to upload";
            }
            return message;
        }, [fileRejections])

        const rejectionMessageStyle =
            fileRejections.length > 0 ? "text-red-500" : "";

        const supportedFileTypes = useMemo(
            () => getSupportedFileExtensions(),
            [],
        );

        const getTextWidth = (
            canvas: HTMLCanvasElement,
            context: CanvasRenderingContext2D,
            text: string,
        ) => {
            context.font = `${getComputedStyle(document.body).fontSize} Arial`;
            return context.measureText(text).width;
        };

        // Function to calculate the maximum total width needed for the ColumnPreview component
        const calculateTotalWidth = (
            columnNames: string[],
            columnTypes: string[],
            secondRowValues: string[],
        ) => {
            const canvas = document.createElement("canvas");
            const context = canvas.getContext("2d");
            if (!context) {
                throw new Error("Could not get 2D context from canvas");
            }

            let maxColumnNameWidth = 0;
            let maxColumnTypeWidth = 0;
            let maxColumnSecondRowWidth = 0;

            // Calculate the maximum width of column names
            maxColumnNameWidth = Math.max(
                ...columnNames.map((name) =>
                    getTextWidth(canvas, context, name),
                ),
            );

            // Calculate the maximum width of column types
            maxColumnTypeWidth = Math.max(
                ...columnTypes.map((type) =>
                    getTextWidth(canvas, context, type),
                ),
            );

            // Calculate the maximum width of second row values
            maxColumnSecondRowWidth = Math.max(
                ...secondRowValues.map((value) =>
                    getTextWidth(canvas, context, value),
                ),
            );

            // Calculate the total width needed for the ColumnPreview component
            const totalWidth =
                maxColumnNameWidth +
                maxColumnTypeWidth +
                maxColumnSecondRowWidth +
                32; // Add padding
            canvas.remove();
            return Math.max(800, totalWidth);
        };

        const handleUploadError = useCallback((error: AxiosError | unknown) => {
            const errorPayload = axios.isAxiosError(error)
                ? {
                      message:
                          error.response?.data ||
                          error.request?.responseText ||
                          "An error occurred while processing the upload.",
                      status: error.response?.status || 500,
                      traceback: error.message,
                  }
                : {
                      message: "An error occurred while processing the upload.",
                      status: 500,
                      traceback:
                          error instanceof Error
                              ? error.message
                              : String(error),
                  };

            setErrorStage(errorPayload);
        }, [setErrorStage]);

        // HTTP Upload function (original implementation)
        const handleHttpUpload = async () => {
            if (!state.selectedFiles.length) {
                setErrorStage({
                    message:
                        "No files selected. Please select a file before uploading.",
                    status: 400,
                });
                return;
            }

            const file = state.selectedFiles[0];
            const fileConfig = getFileTypeFromExtension(file.name);

            if (!fileConfig) {
                setErrorStage({
                    message: "Unsupported file type",
                    status: 400,
                });
                return;
            }

            dispatch({
                type: "SET_STAGE",
                payload: { kind: "uploading" },
            });
            resetProgress();

            const config = {
                headers: {
                    "Content-Type": "multipart/form-data",
                },
                onUploadProgress: (progressEvent: AxiosProgressEvent) => {
                    if (progressEvent.total) {
                        const percentComplete = Math.round(
                            (progressEvent.loaded * 100) / progressEvent.total,
                        );
                        setProgress(percentComplete);
                    } else {
                        // console.warn("Upload progress not computable");
                    }
                },
            };

            try {
                const formData = new FormData();
                formData.append("file", file);

                // Add type-specific form data
                switch (fileConfig.type) {
                    case "tsv":
                    case "csv":
                        formData.append("name", datasourceName);
                        formData.append("replace", "");
                        break;

                    case "tiff":
                        formData.append("datasourceName", datasourceName);
                        if (fileConfig.processingConfig.requiresMetadata) {
                            try {
                                formData.append(
                                    "tiffMetadata",
                                    JSON.stringify(state.tiffMetadata),
                                );
                            } catch (jsonError) {
                                throw new Error(
                                    "Invalid JSON format in tiffMetadata",
                                );
                            }
                        }
                        break;
                    case "h5":
                        // Payload content depends on backend implementation
                        break;
                }

                const response = await axios.post(
                    `${root}/${fileConfig.processingConfig.endpoint}`,
                    formData,
                    config,
                );
                if (response.status === 200) {
                    if (fileConfig.type === "tiff") {
                        chartManager.saveState();
                    }

                    const viewName = response.data?.created_view ?? response.data?.view;
                    if (viewName) {
                        dispatch({
                            type: "SET_STAGE",
                            payload: { kind: "redirecting", viewName },
                        });
                        window.location.href = buildViewRedirectUrl(root, viewName);
                    } else {
                        dispatch({
                            type: "SET_STAGE",
                            payload: { kind: "success" },
                        });
                    }
                } else {
                    throw new Error(
                        `Server responded with status ${response.status}`,
                    );
                }
            } catch (error) {
                if (
                    axios.isAxiosError(error) &&
                    error.response?.status === 409
                ) {
                    dispatch({
                        type: "SET_CONFLICT_DATA",
                        payload: {
                            temp_folder: error.response.data.temp_folder,
                        },
                    });
                    dispatch({
                        type: "SET_STAGE",
                        payload: { kind: "conflict" },
                    });
                } else {
                    console.error("Error uploading file:", error);
                    handleUploadError(error);
                }
            }
        };

        // SocketIO Upload function
        const handleSocketIOUpload = async () => {
            if (!state.selectedFiles.length) {
                setErrorStage({
                    message:
                        "No files selected. Please select a file before uploading.",
                    status: 400,
                });
                return;
            }

            const file = state.selectedFiles[0];
            const fileConfig = getFileTypeFromExtension(file.name);

            if (!fileConfig) {
                setErrorStage({
                    message: "Unsupported file type",
                    status: 400,
                });
                return;
            }

            dispatch({
                type: "SET_STAGE",
                payload: { kind: "uploading" },
            });
            resetProgress();

            try {
                const projectId = getProjectIdFromRoot(root);
                if (!projectId) {
                    throw new Error("Could not determine project ID from current route");
                }

                const socketPath = getSocketPath(String(mainApiRoute));
                console.log('Project ID:', projectId);
                console.log('Namespace:', `/project/${projectId}`);
                console.log('Socket Path:', socketPath);
                
                const isTabularForSocket =
                    fileConfig.type === "csv" || fileConfig.type === "tsv";

                // Create SocketIO upload client
                const uploadClient = createSocketIOUpload({
                    namespace: `/project/${projectId}`,
                    socketPath: socketPath, // socket path for path variable: '/test/socket.io'
                    file: file,
                    fileName: file.name,
                    // Tabular datasource options (CSV / TSV / .txt treated as TSV in UI)
                    datasourceName: isTabularForSocket ? datasourceName : undefined,
                    replace: isTabularForSocket ? true : undefined,
                    view: isTabularForSocket ? "default" : undefined,
                    suppliedOnly: isTabularForSocket ? false : undefined,
                    onProgress: (progressPercent: number, uploaded: number, total: number) => {
                        setProgress(progressPercent);
                    },
                    onStatusChange: (status: string, message?: string) => {
                        console.log(`Upload status: ${status}`, message);
                        if (status === 'processing') {
                            dispatch({
                                type: "SET_STAGE",
                                payload: { kind: "processing" },
                            });
                        }
                    },
                    onError: (error: any) => {
                        console.error('SocketIO upload error:', error);
                        setErrorStage({
                            message: error.message || 'Upload failed',
                            status: 500,
                        });
                    },
                    onSuccess: (result: any) => {
                        console.log('SocketIO upload success:', result);

                        if (fileConfig.type === "tiff") {
                            chartManager.saveState();
                        }

                        const viewName = result?.result?.created_view ?? result?.result?.view;
                        const tabularSuccess =
                            !isTabularForSocket ||
                            (result?.result?.success === true &&
                                Array.isArray(result?.result?.metadata?.columns) &&
                                result.result.metadata.columns.length > 0);

                        if (!tabularSuccess) {
                            setErrorStage({
                                message:
                                    "Upload failed: file format is invalid or could not be parsed. Please verify delimiter and header row.",
                                status: 400,
                            });
                            return;
                        }
                        if (viewName) {
                            dispatch({
                                type: "SET_STAGE",
                                payload: { kind: "redirecting", viewName },
                            });
                            window.location.href = buildViewRedirectUrl(root, viewName);
                        } else {
                            dispatch({
                                type: "SET_STAGE",
                                payload: { kind: "success" },
                            });
                        }

                    }
                });

                // Store client in state for potential cancellation
                dispatch({ type: "SET_SOCKETIO_CLIENT", payload: uploadClient });

                // Start the upload
                await uploadClient.upload();

            } catch (error) {
                console.error("SocketIO upload error:", error);
                setErrorStage({
                    message: error instanceof Error ? error.message : 'Upload failed',
                    status: 500,
                });
            }
        };

        // Main upload handler that chooses between HTTP and SocketIO
        const handleUploadClick = async () => {
            console.log(`Using ${state.uploadMethod} upload method`);
            
            if (state.uploadMethod === 'socketio') {
                await handleSocketIOUpload();
            } else {
                await handleHttpUpload();
            }
        };

        const handleClose = useCallback(async () => {
            if (state.socketioClient) {
                state.socketioClient.cancel();
            }

            dispatch({ type: "RESET_DIALOG" });
            resetLocalState();
            onResize(450, 320);
            onClose();
        }, [onClose, onResize, resetLocalState, state.socketioClient]);

        const handleCombineConfirm = async (label: string) => {
            if (!state.conflictData) return;

            const formData = new FormData();
            formData.append("temp_folder", state.conflictData.temp_folder);
            formData.append("combine", "true");
            formData.append("label", label);

            try {
                const response = await axios.patch(
                    `${root}/combine_anndata`,
                    formData,
                    {
                        headers: {
                            "Content-Type": "multipart/form-data",
                        },
                    },
                );

                if (response.status === 200) {
                    dispatch({
                        type: "SET_STAGE",
                        payload: { kind: "success" },
                    });
                    if (state.fileType === "h5") {
                        chartManager.saveState();
                    }
                }
            } catch (error) {
                console.error("Error replacing file:", error);
                handleUploadError(error);
            } finally {
                dispatch({ type: "SET_CONFLICT_DATA", payload: null });
            }
        };

        const displayUploadDialog = useCallback(() => {
            dispatch({ type: "RESET_DIALOG" });
            resetLocalState();
            onResize(450, 320);
        }, [onResize, resetLocalState]);

        const handleCombineError = (error: {
            message: string;
            status?: number;
        }) => {
            console.error("Conflict dialog error:", error);
            setErrorStage({
                message: error.message,
                status: error.status ?? 500,
            });
        };

        const handleCombineCancel = async () => {
            if (!state.conflictData) return;

            const formData = new FormData();
            formData.append("temp_folder", state.conflictData.temp_folder);
            formData.append("combine", "false");

            try {
                const fullPath = `${root}/combine_anndata`;

                const response = await axios.patch(fullPath, formData, {
                    headers: {
                        "Content-Type": "multipart/form-data",
                        Accept: "application/json, text/plain, */*",
                    },
                });
                console.log("Response:", response);
            } catch (error) {
                console.error(
                    "Full URL that failed:",
                    `${root}/combine_anndata`,
                );
                if (axios.isAxiosError(error)) {
                    console.error("Response data:", error.response?.data);
                    console.error("Response status:", error.response?.status);
                }
                throw error;
            } finally {
                dispatch({ type: "RESET_DIALOG" });
                resetLocalState();
                onResize(450, 320);
                onClose();
            }
        };

        // Method toggle function for testing/debugging
        const toggleUploadMethod = () => {
            const newMethod = state.uploadMethod === "http" ? "socketio" : "http";
            dispatch({ type: "SET_UPLOAD_METHOD", payload: newMethod });
            console.log(`Switched to ${newMethod} upload method`);
        };

        const isReadyForReview =
            state.stage.kind === "ready" && state.validationResult !== null;
        const currentError =
            state.stage.kind === "error" ? state.stage.error : null;

        return (
            <Container>
                {/* Upload Method Toggle (for development/testing), uncomment when required */}
                {/* {process.env.NODE_ENV === 'development' && (
                    <div className="mb-4 p-2 bg-gray-100 dark:bg-gray-700 rounded">
                        <button
                            type="button"
                            onClick={toggleUploadMethod}
                            className="text-sm px-3 py-1 bg-blue-500 text-white rounded hover:bg-blue-600"
                        >
                            Current: {state.uploadMethod.toUpperCase()} - Click to toggle
                        </button>
                    </div>
                )} */}
                {state.stage.kind === "redirecting" ? (
                    <RedirectingView />
                ) : state.stage.kind === "uploading" ? (
                    <UploadingView progress={progress} />
                ) : state.stage.kind === "processing" ? (
                    <ProcessingView />
                ) : state.stage.kind === "success" ? (
                    <UploadSuccessView onRefresh={() => window.location.reload()} />
                ) : currentError ? (
                    <UploadErrorView
                        error={currentError}
                        onUploadAnother={displayUploadDialog}
                        onCancel={handleClose}
                    />
                ) : state.stage.kind === "conflict" ? (
                    <AnndataConflictDialog
                        open={true}
                        onClose={handleCombineCancel}
                        onConfirm={(label) => handleCombineConfirm(label)}
                        onError={handleCombineError}
                    />
                ) : state.stage.kind === "validating" ? (
                    <ValidatingView />
                ) : isReadyForReview ? (
                    <>
                        {(state.fileType === "csv" ||
                            state.fileType === "tsv") && (
                            <TabularReviewView
                                datasourceName={datasourceName}
                                datasourceSummary={datasourceSummary}
                                onDatasourceNameChange={handleDatasourceNameChange}
                                isDatasourceDisabled={state.selectedFiles.length === 0}
                                onUpload={handleUploadClick}
                                onCancel={handleClose}
                                columnPreview={(
                                    <>
                                        {/* <ColumnPreview
                                            columnNames={columnNames}
                                            columnTypes={columnTypes}
                                            secondRowValues={secondRowValues}
                                        /> */}
                                    </>
                                )}
                            />
                        )}

                        {state.fileType === "tiff" && state.tiffMetadata && (
                            <TiffReviewView
                                tiffSummary={tiffSummary}
                                tiffMetadata={state.tiffMetadata}
                                selectedFile={state.selectedFiles[0]}
                                showMetadata={showMetadata}
                                onToggleView={toggleView}
                                onUpload={handleUploadClick}
                                onCancel={handleClose}
                                datasourceDropdown={(
                                    <>
                                        {/* <DatasourceDropdown
                                            options={
                                                updatedNamesArray
                                            }
                                            onSelect={handleSelect}
                                        /> */}
                                    </>
                                )}
                            />
                        )}

                        {state.fileType === "h5" && (
                            <H5ReviewView
                                selectedFile={state.selectedFiles[0]}
                                h5Metadata={state.h5Metadata}
                                onUpload={handleUploadClick}
                                onCancel={handleClose}
                            />
                        )}
                    </>
                ) : (
                    <FileDropzoneView
                        rootProps={getRootProps()}
                        inputProps={getInputProps()}
                        isDragActive={isDragActive}
                        selectedFileName={state.selectedFiles[0]?.name ?? null}
                        rejectionMessage={rejectionMessage}
                        rejectionMessageStyle={rejectionMessageStyle}
                        supportedFileTypes={supportedFileTypes}
                    />
                )}
            </Container>
        );
    });

const Wrapper = (props: Omit<FileUploadDialogComponentProps, 'onLoadingStateChange' | 'socketioClientRef'>) => {
    const [open, setOpen] = useState(true);
    const [, setIsLoading] = useState(false);
    const [socketioClient, setSocketioClient] = useState<SocketIOUploadClient | null>(null);

    const handleClose = () => {
        if (socketioClient) {
            socketioClient.cancel();
        }
        setOpen(false);
        props.onClose();
    };
    
    return (
        <ReusableDialog 
            open={open} 
            handleClose={handleClose} 
            // Pass component with the new callback prop
            component={
                <FileUploadDialogComponent 
                    {...props} 
                    onClose={handleClose} 
                    onLoadingStateChange={setIsLoading}
                    socketioClientRef={setSocketioClient}
                />
            } 
        />
    );
};

export default Wrapper;
