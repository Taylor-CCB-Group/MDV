import type React from "react";
import {
    useState,
    useCallback,
    useReducer,
    type PropsWithChildren,
    forwardRef
} from "react";
import { useDropzone } from "react-dropzone";
import { observer } from "mobx-react-lite";

import axios, { type AxiosError, type AxiosProgressEvent } from "axios";
import { useProject } from "../../modules/ProjectContext";
import { ColumnPreview } from "./ColumnPreview";

import {
    useViewerStoreApi,
    useChannelsStoreApi,
    type PixelSource,
    type Metadata,
} from "../../react/components/avivatorish/state";
import { createLoader } from "../../react/components/avivatorish/utils";
import { unstable_batchedUpdates } from "react-dom";

import { TiffPreview } from "./TiffPreview";
import { TiffMetadataTable } from "./TiffMetadataTable";
import TiffVisualization from "./TiffVisualization";
import { DatasourceDropdown } from "./DatasourceDropdown";
import {
    CloudUpload as CloudUploadIcon
} from "@mui/icons-material";
import { isArray } from "@/lib/utils";
// import processH5File from "./utils/h5Processing";
import H5MetadataPreview from "./H5MetadataPreview";
import DebugErrorComponent from "./DebugErrorComponent";
import AnndataConflictDialog from "./AnndataConflictDialog";
import ReusableDialog from "./ReusableDialog";

// Use dynamic import for the worker
const DatasourceWorker = new Worker(
    new URL("./datasourceWorker.ts", import.meta.url),
    {
        type: "module",
    },
);

export const Container = ({ children }: PropsWithChildren) => {
    return (
        <div className="flex flex-col content-center items-center h-max dark:bg-black-800 dark:text-white">
            {children}
        </div>
    );
};

const StatusContainer = ({ children }: PropsWithChildren) => {
    return (
        <div className="flex flex-col justify-center items-center w-full h-150 dark:bg-#333">
            {children}
        </div>
    );
};

const SuccessContainer = ({ children }: PropsWithChildren) => (
    <div className="flex flex-col items-center justify-center bg-[#f0f8ff] shadow-md border border-[#e0e0e0] m-4 dark:bg-black dark:border-gray-600">
        {children}
    </div>
);

const SuccessHeading = ({ children }: PropsWithChildren) => (
    <h1 className="text-[#333] mb-1 dark:text-white">{children}</h1>
);

const SuccessText = ({ children }: PropsWithChildren) => (
    <p className="text-2xl text-[#555] mb-3 text-center dark:text-gray-300">
        {children}
    </p>
);

export const DropzoneContainer = forwardRef(
    ({ isDragOver, children, ...props }: any, ref) => (
        <div
            {...props}
            ref={ref}
            className={`p-4 mt-2 z-50 text-center border-2 border-dashed rounded-lg ${isDragOver ? "bg-gray-300 dark:bg-slate-800" : "bg-white dark:bg-black"} min-w-[90%]`}
        >
            {children}
        </div>
    ),
);

//! warning: className isn't being merged etc, I think it just gets overriden.
// not sure if there's a better type for `htmlFor`.
export const FileInputLabel = ({
    children,
    htmlFor,
    ...props
}: PropsWithChildren & { className: string; htmlFor: string }) => (
    <label
        htmlFor={htmlFor}
        {...props}
        className="mt-8 px-5 py-2.5 border bg-stone-200 hover:bg-stone-300 rounded cursor-pointer inline-block my-2.5 dark:bg-stone-600 dark:hover:bg-stone-500"
    >
        {children}
    </label>
);

const Spinner = () => {
    return (
        <div
            className="w-16 h-16 border-8 mt-10 border-blue-500 border-dashed rounded-full animate-spin"
            style={{
                borderColor: "blue transparent blue transparent",
            }}
        />
    );
};

const colorStyles = {
    blue: {
        bgColor: "bg-blue-600",
        hoverColor: "hover:bg-blue-700",
        darkBgColor: "dark:bg-blue-800",
        darkHoverColor: "dark:bg-blue-900",
    },
    red: {
        bgColor: "bg-red-600",
        hoverColor: "hover:bg-red-700",
        darkBgColor: "dark:bg-red-800",
        darkHoverColor: "dark:bg-red-900",
    },
    green: {
        bgColor: "bg-green-600",
        hoverColor: "hover:bg-green-700",
        darkBgColor: "dark:bg-green-800",
        darkHoverColor: "dark:bg-green-900",
    },
    gray: {
        bgColor: "bg-gray-600",
        hoverColor: "hover:bg-gray-700",
        darkBgColor: "dark:bg-gray-800",
        darkHoverColor: "dark:bg-gray-900",
    },
};
type ButtonProps = {
    onClick: () => void;
    color?: "blue" | "red" | "green" | "gray";
    disabled?: boolean;
    size?: string;
    marginTop?: string;
} & PropsWithChildren;
export const Button = ({
    onClick,
    color = "blue",
    disabled = false,
    size = "px-5 py-2.5",
    marginTop = "mt-2.5",
    children,
}: ButtonProps) => {
    const { bgColor, hoverColor, darkBgColor, darkHoverColor } =
        colorStyles[color] || colorStyles.blue;

    return (
        <button
            type="button"
            onClick={onClick}
            className={`${size} ${marginTop} ${bgColor} ${hoverColor} ${darkBgColor} ${darkHoverColor} text-white rounded self-center cursor-pointer disabled:cursor-not-allowed disabled:bg-gray-500 disabled:opacity-70`}
            disabled={disabled}
        >
            {children}
        </button>
    );
};

//@ts-ignore CBA
const ProgressBar = ({ value, max }) => (
    <progress
        className="w-full h-8 mb-5 mt-10 bg-gray-200 dark:bg-white-200 border border-gray-300 rounded"
        value={value}
        max={max}
    />
);

const Message = ({ children }: PropsWithChildren) => (
    <p className="text-lg font-bold text-[#333] dark:text-white text-center">
        {children}
    </p>
);

const FileSummary = ({ children }: PropsWithChildren) => (
    <div className="max-w-[800px] w-[100%] text-center mb-5">{children}</div>
);

const FileSummaryHeading = ({ children }: PropsWithChildren) => (
    <h1 className="text-gray-800 dark:text-white mb-3">{children}</h1>
);

const FileSummaryText = ({ children }: PropsWithChildren) => (
    <>
        {typeof children === "string" ? (
            <p className="text-lg text-gray-700 dark:text-white my-1">
                {children}
            </p>
        ) : (
            <div className="text-lg text-gray-700 dark:text-white my-1">
                {children}
            </div>
        )}
    </>
);

const ErrorContainer = ({ children }: PropsWithChildren) => (
    <div className="bg-[#fff0f0] text-[#d8000c] border border-[#ffbaba] p-5 my-5 rounded shadow-sm text-left w-[90%] max-w-[600px]">
        {children}
    </div>
);

export const DynamicText = ({
    text,
    className = "",
}: { text: string; className: string }) => (
    <div className="w-96 h-20 overflow-hidden flex items-center justify-center">
        <p
            //consider using `cn()` from lib/utils
            className={`text-center m-0 font-bold text-sm sm:text-lg md:text-m ${className}`}
        >
            {text}
        </p>
    </div>
);

const ErrorHeading = ({ children }: PropsWithChildren) => (
    <h4 className="mt-0 text-lg font-bold">{children}</h4>
);
//@ts-ignore CBA
const DatasourceNameInput = ({ value, onChange, isDisabled }) => (
    <div className="flex-left items-center space-x-2 pr-4">
        <label className="text-lg text-gray-700 dark:text-white my-1" htmlFor="datasourceName">
            <strong>Datasouce Name:</strong>
        </label>
        <input
            id="datasourceName"
            type="text"
            value={value}
            onChange={onChange}
            className="p-2 border rounded focus:outline-none focus:ring focus:border-blue-300 dark:text-gray-300 dark:bg-gray-800 flex-grow"
            placeholder="Enter datasource name"
            disabled={isDisabled}
        />
    </div>
);

type UploadActionType =
    | "SET_SELECTED_FILES"
    | "SET_IS_UPLOADING"
    | "SET_IS_INSERTING"
    | "SET_SUCCESS"
    | "SET_ERROR"
    | "SET_IS_VALIDATING"
    | "SET_VALIDATION_RESULT"
    | "SET_FILE_TYPE"
    | "SET_TIFF_METADATA"
    | "SET_H5_METADATA"
    | "SET_FILE_SUMMARY"
    | "SET_CONFLICT_DATA"
    | "SET_SHOW_CONFLICT_DIALOG";

// Reducer function
const DEFAULT_REDUCER_STATE = {
    selectedFiles: [] as File[],
    isUploading: false,
    isInserting: false,
    isValidating: false,
    validationResult: null as unknown,
    success: false,
    error: null as unknown,
    fileType: null as "csv" | "tiff" | "tsv" | "h5" | null,
    tiffMetadata: null as unknown,
    h5Metadata: null as H5Metadata | unknown,
    showReplaceDialog: false,
    conflictData: null as { tempFolder: string } | null,
} as const;
type ReducerState = typeof DEFAULT_REDUCER_STATE;
// TODO - would be good to type this, this is not how I should be spending my weekend.
//                                    ^^ only myself to blame for this - I should get a life
type ReducerPayload<T extends UploadActionType> = T extends "SET_SELECTED_FILES"
    ? File[]
    : T extends "SET_IS_UPLOADING"
      ? boolean
      : T extends "SET_IS_INSERTING"
        ? boolean
        : T extends "SET_SUCCESS"
          ? boolean
          : T extends "SET_ERROR"
            ? { message: string; status: number; traceback?: string }
            : T extends "SET_IS_VALIDATING"
              ? boolean
              : T extends "SET_VALIDATION_RESULT"
                ? { columnNames: string[]; columnTypes: string[] }
                : T extends "SET_FILE_TYPE"
                  ? "csv" | "tiff" | "tsv" | null
                  : T extends "SET_TIFF_METADATA"
                    ? unknown
                    : T extends "SET_CONFLICT_DATA"
                      ? { temp_folder: string } | null
                      : T extends "SET_SHOW_CONFLICT_DIALOG"
                        ? boolean
                        : never;
type ReducerAction<T extends UploadActionType> = {
    type: T;
    payload: ReducerPayload<T>;
};
const reducer = <T extends UploadActionType>(
    state: ReducerState,
    action: { type: T; payload: any },
) => {
    switch (action.type) {
        case "SET_SELECTED_FILES":
            return { ...state, selectedFiles: action.payload };
        case "SET_IS_UPLOADING":
            return { ...state, isUploading: action.payload };
        case "SET_IS_INSERTING":
            return { ...state, isInserting: action.payload };
        case "SET_SUCCESS":
            return { ...state, success: action.payload };
        case "SET_ERROR":
            return { ...state, error: action.payload };
        case "SET_IS_VALIDATING":
            return { ...state, isValidating: action.payload };
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
        case "SET_SHOW_CONFLICT_DIALOG":
            return { ...state, showReplaceDialog: action.payload };
        default:
            return state;
    }
};

// Constants
const UPLOAD_INTERVAL = 30;
const UPLOAD_DURATION = 3000;
const UPLOAD_STEP = 100 / (UPLOAD_DURATION / UPLOAD_INTERVAL);

export interface FileUploadDialogComponentProps {
    onClose: () => void;
    onResize: (width: number, height: number) => void; // Add this prop
}

// Define supported file types and their configurations
export interface FileTypeConfig {
    type: string;
    extensions: string[];
    mimeTypes: string[];
    maxSize?: number; // in bytes
    processingConfig: {
        defaultWidth: number;
        defaultHeight: number;
        requiresMetadata?: boolean;
        endpoint: string;
    };
}

const FILE_TYPES: Record<string, FileTypeConfig> = {
    CSV: {
        type: "csv",
        extensions: [".csv"],
        mimeTypes: ["text/csv"],
        maxSize: 10000 * 1024 * 1024, // 10 GB
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
        mimeTypes: ["text/tab-separated-values"],
        maxSize: 10000 * 1024 * 1024, // 10 GB
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
        maxSize: 10000 * 1024 * 1024, // 10 GB
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

interface H5Metadata {
    uns: Record<string, any>;
    obs: Record<string, any>;
    var: Record<string, any>;
}

// Helper functions for file type checking
const getFileTypeFromExtension = (fileName: string): FileTypeConfig | null => {
    const extension = `.${fileName.split(".").pop()?.toLowerCase()}`;
    return (
        Object.values(FILE_TYPES).find((config) =>
            config.extensions.includes(extension),
        ) || null
    );
};

const generateDropzoneAccept = () => {
    return Object.values(FILE_TYPES).reduce(
        (acc, config) => {
            config.mimeTypes.forEach((mimeType) => {
                acc[mimeType] = config.extensions;
            });
            return acc;
        },
        {} as Record<string, string[]>,
    );
};

// Custom hook for file upload progress
const useFileUploadProgress = () => {
    const [progress, setProgress] = useState(0);

    const startProgress = useCallback(() => {
        setProgress(0);
        const interval = setInterval(() => {
            setProgress((prevProgress) => {
                const updatedProgress = prevProgress + UPLOAD_STEP;
                if (updatedProgress >= 100) {
                    clearInterval(interval);
                    return 100;
                }
                return updatedProgress;
            });
        }, UPLOAD_INTERVAL);
    }, []);

    const resetProgress = useCallback(() => {
        setProgress(0);
    }, []);

    return { progress, setProgress, startProgress, resetProgress };
};

//todo figure out why type inference is failing with observer if we take out the `React.FC<...>` here
//also - are we actually observing anything here? not that it particularly matters (apart from type inference, apparently)
const FileUploadDialogComponent: React.FC<FileUploadDialogComponentProps> =
    observer(({ onClose, onResize }: FileUploadDialogComponentProps) => {
        const { root, projectName, chartManager } = useProject();

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

        const [datasourceSummary, setDatasourceSummary] = useState({
            datasourceName: "",
            fileName: "",
            fileSize: "",
            rowCount: 0,
            columnCount: 0,
        });

        const [tiffSummary, settiffSummary] = useState({
            fileName: "",
            fileSize: "",
        });

        const [columnNames, setColumnNames] = useState<string[]>([]);
        const [columnTypes, setColumnTypes] = useState<string[]>([]);
        const [secondRowValues, setSecondRowValues] = useState<string[]>([]);

        // TIFF
        const channelsStore = useChannelsStoreApi();

        const [state, dispatch] = useReducer(reducer, DEFAULT_REDUCER_STATE);

        const viewerStore = useViewerStoreApi();

        const { progress, resetProgress, setProgress } =
            useFileUploadProgress();

        const onDrop = useCallback(
            async (acceptedFiles: File[]) => {
                if (acceptedFiles.length > 0) {
                    const file = acceptedFiles[0];
                    const fileConfig = getFileTypeFromExtension(file.name);

                    if (!fileConfig) {
                        dispatch({
                            type: "SET_ERROR",
                            payload: {
                                message: "Unsupported file type",
                                status: 400,
                            },
                        });
                        return;
                    }

                    dispatch({
                        type: "SET_SELECTED_FILES",
                        payload: acceptedFiles,
                    });
                    dispatch({ type: "SET_IS_VALIDATING", payload: true });

                    // Handle file based on its type
                    switch (fileConfig.type) {
                        case "tsv":
                        case "csv": {
                            const newDatasourceName = file.name;
                            setDatasourceName(newDatasourceName);
                            setDatasourceSummary({
                                datasourceName: newDatasourceName,
                                fileName: file.name,
                                fileSize: (file.size / (1024 * 1024)).toFixed(
                                    2,
                                ),
                                rowCount: 0,
                                columnCount: 0,
                            });

                            // Process CSV with Web Worker
                            DatasourceWorker.postMessage({
                                file: file,
                                fileType: fileConfig.type,
                            });
                            DatasourceWorker.onmessage = (
                                event: MessageEvent,
                            ) => {
                                const {
                                    columnNames,
                                    columnTypes,
                                    secondRowValues,
                                    rowCount,
                                    columnCount,
                                    error,
                                } = event.data;

                                if (error) {
                                    dispatch({
                                        type: "SET_ERROR",
                                        payload: {
                                            message: "Validation failed.",
                                            traceback: error,
                                        },
                                    });
                                    dispatch({
                                        type: "SET_IS_VALIDATING",
                                        payload: false,
                                    });
                                } else {
                                    setColumnNames(columnNames);
                                    setColumnTypes(columnTypes);
                                    setSecondRowValues(secondRowValues);
                                    setDatasourceSummary(
                                        (prevDatasrouceSummary) => ({
                                            ...prevDatasrouceSummary,
                                            rowCount,
                                            columnCount,
                                        }),
                                    );

                                    const totalWidth = calculateTotalWidth(
                                        columnNames,
                                        columnTypes,
                                        secondRowValues,
                                    );
                                    onResize(totalWidth, 745);
                                    dispatch({
                                        type: "SET_IS_VALIDATING",
                                        payload: false,
                                    });
                                    dispatch({
                                        type: "SET_VALIDATION_RESULT",
                                        payload: { columnNames, columnTypes },
                                    });
                                }
                            };
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
                            settiffSummary({
                                fileName: file.name,
                                fileSize: (file.size / (1024 * 1024)).toFixed(
                                    2,
                                ),
                            });

                            // Process TIFF file
                            handleSubmitFile(acceptedFiles);
                            onResize(
                                fileConfig.processingConfig.defaultWidth,
                                fileConfig.processingConfig.defaultHeight,
                            );
                            dispatch({
                                type: "SET_IS_VALIDATING",
                                payload: false,
                            });
                            dispatch({
                                type: "SET_VALIDATION_RESULT",
                                payload: { columnNames, columnTypes },
                            });
                            break;
                        }
                        case "h5": {
                            try {
                                const processH5File = (await import(
                                    "./utils/h5Processing"
                                )).default;
                                const h5Metadata = await processH5File(file);
                                dispatch({
                                    type: "SET_H5_METADATA",
                                    payload: h5Metadata,
                                });

                                const newDatasourceName = file.name;
                                setDatasourceName(newDatasourceName);

                                onResize(
                                    fileConfig.processingConfig.defaultWidth,
                                    fileConfig.processingConfig.defaultHeight,
                                );
                                dispatch({
                                    type: "SET_IS_VALIDATING",
                                    payload: false,
                                });
                                dispatch({
                                    type: "SET_VALIDATION_RESULT",
                                    payload: { metadata: h5Metadata },
                                });
                            } catch (error) {
                                console.error("H5 processing error:", error);
                                dispatch({
                                    type: "SET_ERROR",
                                    payload: {
                                        message: "Error processing H5 file",
                                        traceback:
                                            error instanceof Error
                                                ? error.message
                                                : String(error),
                                    },
                                });
                                dispatch({
                                    type: "SET_IS_VALIDATING",
                                    payload: false,
                                });
                            }
                            break;
                        }
                        default:
                            console.warn(
                                "Unhandled file type:",
                                fileConfig.type,
                            );
                            return;
                    }

                    // Set file type in state
                    dispatch({
                        type: "SET_FILE_TYPE",
                        payload: fileConfig.type,
                    });
                }
            },
            [onResize, columnNames, columnTypes],
        );

        const handleSubmitFile = useCallback(
            async (files: File[]) => {
                let newSource;
                if (files.length === 1) {
                    newSource = {
                        urlOrFile: files[0],
                        description: files[0].name,
                    };
                } else {
                    newSource = {
                        urlOrFile: files,
                        description: "data.zarr",
                    };
                }

                viewerStore.setState({ isChannelLoading: [true] });
                viewerStore.setState({ isViewerLoading: true });

                try {
                    const newLoader = await createLoader(
                        newSource.urlOrFile,
                        () => {}, // placeholder for toggleIsOffsetsSnackbarOn
                        (message) => {
                            viewerStore.setState({
                                loaderErrorSnackbar: { on: true, message },
                            });
                        },
                    );
                    // types are still somewhat on a wing and a prayer... we'll get there
                    let nextMeta: Metadata;
                    let nextLoader: PixelSource["data"] | PixelSource["data"][];
                    //let nextLoader: OME_TIFF['data'] | OME_TIFF['data'][];
                    if (isArray(newLoader)) {
                        if (newLoader.length > 1) {
                            //@ts-ignore flagging so we can get back to this
                            nextMeta = newLoader.map((l) => l.metadata);
                            nextLoader = newLoader.map((l) => l.data);
                        } else {
                            nextMeta = newLoader[0].metadata;
                            nextLoader = newLoader[0].data;
                        }
                    } else {
                        //@ts-ignore flagging so we can get back to this
                        nextMeta = newLoader.metadata;
                        nextLoader = newLoader.data;
                    }

                    if (nextLoader) {
                        console.log(
                            "Metadata (in JSON-like form) for current file being viewed: ",
                            nextMeta,
                        );

                        unstable_batchedUpdates(() => {
                            channelsStore.setState({ loader: nextLoader });
                            viewerStore.setState({ metadata: nextMeta });
                        });

                        dispatch({
                            type: "SET_TIFF_METADATA",
                            payload: nextMeta,
                        });
                    }
                } catch (error: any) {
                    console.error("Error loading file:", error);
                    dispatch({
                        type: "SET_ERROR",
                        payload: {
                            message: "Error loading TIFF file.",
                            traceback: error.message,
                        },
                    });
                } finally {
                    viewerStore.setState({ isChannelLoading: [false] });
                    viewerStore.setState({ isViewerLoading: false });
                }
            },
            [viewerStore, channelsStore],
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
                maxSize: Math.max(
                    ...Object.values(FILE_TYPES).map(
                        (config) => config.maxSize || 0,
                    ),
                ),
            });

        const rejectionMessage =
            fileRejections.length > 0
                ? "Only CSV, TSV, and TIFF files can be selected"
                : "Drag and drop files here or click the button below to upload";

        const rejectionMessageStyle =
            fileRejections.length > 0 ? "text-red-500" : "";

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
            if (!context)
                throw new Error("Could not get 2D context from canvas");

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

        const handleUploadClick = async () => {
            if (!state.selectedFiles.length) {
                dispatch({
                    type: "SET_ERROR",
                    payload: {
                        message:
                            "No files selected. Please select a file before uploading.",
                        status: 400,
                    },
                });
                return;
            }

            const file = state.selectedFiles[0];
            const fileConfig = getFileTypeFromExtension(file.name);

            if (!fileConfig) {
                dispatch({
                    type: "SET_ERROR",
                    payload: {
                        message: "Unsupported file type",
                        status: 400,
                    },
                });
                return;
            }

            dispatch({ type: "SET_IS_UPLOADING", payload: true });
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
                    dispatch({ type: "SET_IS_UPLOADING", payload: false });
                    dispatch({ type: "SET_SUCCESS", payload: true });

                    if (fileConfig.type === "tiff") {
                        chartManager.saveState();
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
                    dispatch({ type: "SET_IS_UPLOADING", payload: false });
                    dispatch({ type: "SET_SUCCESS", payload: false });
                    dispatch({
                        type: "SET_CONFLICT_DATA",
                        payload: {
                            AnndataConflictDialog,
                            temp_folder: error.response.data.temp_folder,
                        },
                    });
                    dispatch({
                        type: "SET_SHOW_CONFLICT_DIALOG",
                        payload: true,
                    });
                } else {
                    console.error("Error uploading file:", error);
                    handleUploadError(error);
                }
            }
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
                      traceback: (error as Error).message,
                  };

            dispatch({
                type: "SET_ERROR",
                payload: errorPayload,
            });
            dispatch({ type: "SET_IS_UPLOADING", payload: false });
        }, []);

        const handleClose = useCallback(async () => {
            dispatch({ type: "SET_FILE_SUMMARY", payload: null });
            onResize(450, 320);
            onClose();
        }, [onClose, onResize]);

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
                    dispatch({ type: "SET_SUCCESS", payload: true });
                    if (state.fileType === "h5") {
                        chartManager.saveState();
                    }
                }
            } catch (error) {
                console.error("Error replacing file:", error);
                handleUploadError(error);
            } finally {
                dispatch({ type: "SET_SHOW_CONFLICT_DIALOG", payload: false });
                dispatch({ type: "SET_CONFLICT_DATA", payload: null });
            }
        };

        const displayUploadDialog = async () => {
            if (state.error) {
                dispatch({ type: "SET_SELECTED_FILES", payload: [] });
                dispatch({ type: "SET_IS_UPLOADING", payload: false });
                dispatch({ type: "SET_IS_INSERTING", payload: false });
                dispatch({ type: "SET_SUCCESS", payload: false });
                dispatch({ type: "SET_ERROR", payload: null });
                dispatch({ type: "SET_IS_VALIDATING", payload: false });
                dispatch({ type: "SET_VALIDATION_RESULT", payload: null });
                dispatch({ type: "SET_FILE_TYPE", payload: null });
                dispatch({ type: "SET_TIFF_METADATA", payload: null });
                dispatch({ type: "SET_H5_METADATA", payload: null });
                dispatch({ type: "SET_CONFLICT_DATA", payload: null });
                dispatch({ type: "SET_SHOW_CONFLICT_DIALOG", payload: false });
                setDatasourceName("");
                setDatasourceSummary({
                    datasourceName: "",
                    fileName: "",
                    fileSize: "",
                    rowCount: 0,
                    columnCount: 0,
                });
                return;
            }
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
                dispatch({ type: "SET_SHOW_CONFLICT_DIALOG", payload: false });
                dispatch({ type: "SET_CONFLICT_DATA", payload: null });
                dispatch({ type: "SET_VALIDATION_RESULT", payload: null });
                dispatch({ type: "SET_H5_METADATA", payload: null });
                dispatch({ type: "SET_SELECTED_FILES", payload: [] });
                dispatch({ type: "SET_FILE_TYPE", payload: null });
                dispatch({ type: "SET_FILE_SUMMARY", payload: null });
                onResize(450, 320);
                onClose();
            }
        };

        const handleCombineError = (error: {
            message: string;
            status?: number;
        }) => {
            console.error("Conflict dialog error:", error);
            dispatch({ type: "SET_ERROR", payload: error });
        };

        return (
            <Container>
                {state.isUploading ? (
                    <StatusContainer>
                        <Message>
                            {"Your file is being uploaded, please wait..."}
                        </Message>
                        <ProgressBar value={progress} max="100" />
                    </StatusContainer>
                ) : state.isInserting ? (
                    <StatusContainer>
                        <Message>
                            {"Your file is being processed, please wait..."}
                        </Message>
                        <Spinner />
                    </StatusContainer>
                ) : state.success ? (
                    <>
                        <SuccessContainer>
                            <SuccessHeading>Success!</SuccessHeading>
                            <SuccessText>
                                The file was uploaded successfully to the
                                database.
                            </SuccessText>
                        </SuccessContainer>
                        <Button
                            color="green"
                            onClick={() => window.location.reload()}
                        >
                            Refresh Page
                        </Button>
                    </>
                ) : state.error ? (
                    <>
                        <DebugErrorComponent error={state.error} />
                        <div className="flex justify-center items-center gap-6">
                            <Button
                                marginTop="mt-1"
                                onClick={displayUploadDialog}
                            >
                                Upload Another File
                            </Button>
                            <Button
                                color="red"
                                size="px-14 py-2.5"
                                marginTop="mt-1"
                                onClick={handleClose}
                            >
                                Cancel
                            </Button>
                        </div>
                    </>
                ) : state.showReplaceDialog ? (
                    <AnndataConflictDialog
                        open={state.showReplaceDialog}
                        onClose={handleCombineCancel}
                        onConfirm={(label) => handleCombineConfirm(label)}
                        onError={handleCombineError}
                    />
                ) : state.isValidating ? (
                    <StatusContainer>
                        <Message>{"Validating data, please wait..."}</Message>
                        <Spinner />
                    </StatusContainer>
                ) : state.validationResult ? (
                    <>
                        {(state.fileType === "csv" ||
                            state.fileType === "tsv") && (
                            <>
                                <FileSummary>
                                    <FileSummaryHeading>
                                        {"Uploaded File Summary"}
                                    </FileSummaryHeading>
                                    <FileSummaryText>
                                        <DatasourceNameInput
                                            value={datasourceName}
                                            onChange={
                                                handleDatasourceNameChange
                                            }
                                            isDisabled={
                                                state.selectedFiles.length === 0
                                            }
                                        />
                                    </FileSummaryText>
                                    <FileSummaryText>
                                        <strong>{"File name"}</strong>{" "}
                                        {datasourceSummary.fileName}
                                    </FileSummaryText>
                                    <FileSummaryText>
                                        <strong>{"Number of rows"}</strong>{" "}
                                        {datasourceSummary.rowCount}
                                    </FileSummaryText>
                                    <FileSummaryText>
                                        <strong>{"Number of columns"}</strong>{" "}
                                        {datasourceSummary.columnCount}
                                    </FileSummaryText>
                                    <FileSummaryText>
                                        <strong>{"File size"}</strong>{" "}
                                        {datasourceSummary.fileSize} MB
                                    </FileSummaryText>
                                </FileSummary>
                                <ColumnPreview
                                    columnNames={columnNames}
                                    columnTypes={columnTypes}
                                    secondRowValues={secondRowValues}
                                />

                                <div className="flex justify-center items-center gap-6 mt-4">
                                    <Button
                                        marginTop="mt-1"
                                        onClick={handleUploadClick}
                                    >
                                        {"Upload"}
                                    </Button>
                                    <Button
                                        color="red"
                                        size="px-6 py-2.5"
                                        marginTop="mt-1"
                                        onClick={handleClose}
                                    >
                                        {"Cancel"}
                                    </Button>
                                </div>
                            </>
                        )}

                        {state.fileType === "tiff" && state.tiffMetadata && (
                            <>
                                <div className="h-full w-full flex">
                                    <div className="w-1/3 flex flex-col">
                                        <div className="flex items-center justify-center">
                                            <div className="flex items-start justify-start pt-5 pl-4">
                                                <FileSummary>
                                                    <FileSummaryHeading>
                                                        {
                                                            "Uploaded File Summary"
                                                        }
                                                    </FileSummaryHeading>
                                                    <FileSummaryText>
                                                        <strong>
                                                            {"File name:"}
                                                        </strong>{" "}
                                                        {tiffSummary.fileName}
                                                    </FileSummaryText>
                                                    <FileSummaryText>
                                                        <strong>
                                                            {"File size:"}
                                                        </strong>{" "}
                                                        {tiffSummary.fileSize}{" "}
                                                        MB
                                                    </FileSummaryText>
                                                    <DatasourceDropdown
                                                        options={
                                                            updatedNamesArray
                                                        }
                                                        onSelect={handleSelect}
                                                    />
                                                </FileSummary>
                                            </div>
                                        </div>
                                        <div className="flex items-center justify-center">
                                            <div className="flex flex-col items-center justify-center ">
                                                <FileSummaryText>
                                                    <strong>
                                                        Image Preview:
                                                    </strong>
                                                </FileSummaryText>
                                                <TiffVisualization
                                                    metadata={
                                                        state.tiffMetadata
                                                    }
                                                    file={
                                                        state.selectedFiles[0]
                                                    }
                                                />
                                            </div>
                                        </div>
                                    </div>
                                    <div className="w-2/3">
                                        <div className="flex items-start justify-start h-[70%] min-h-[420px]">
                                            {showMetadata ? (
                                                <TiffMetadataTable
                                                    metadata={
                                                        state.tiffMetadata
                                                    }
                                                />
                                            ) : (
                                                <TiffPreview
                                                    metadata={
                                                        state.tiffMetadata
                                                    }
                                                />
                                            )}
                                        </div>
                                        <div className="flex justify-between items-center ml-4 mr-2 mt-12">
                                            <Button
                                                color="gray"
                                                size="px-6 py-2.5"
                                                marginTop="mt-1"
                                                onClick={toggleView}
                                            >
                                                {showMetadata
                                                    ? "View Metadata as XML"
                                                    : "View Metadata as Table"}
                                            </Button>
                                            <div className="flex space-x-4 ml-auto">
                                                <Button
                                                    marginTop="mt-6"
                                                    onClick={handleUploadClick}
                                                >
                                                    {"Upload"}
                                                </Button>
                                                <Button
                                                    color="red"
                                                    size="px-6 py-2.5"
                                                    marginTop="mt-6"
                                                    onClick={handleClose}
                                                >
                                                    {"Cancel"}
                                                </Button>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </>
                        )}

                        {state.fileType === "h5" && state.h5Metadata && (
                            <>
                                <FileSummary>
                                    <FileSummaryHeading>
                                        H5 File Summary
                                    </FileSummaryHeading>
                                    <FileSummaryText>
                                        <strong>File name:</strong>{" "}
                                        {state.selectedFiles[0].name}
                                    </FileSummaryText>
                                    <FileSummaryText>
                                        <strong>File size:</strong>{" "}
                                        {(
                                            state.selectedFiles[0].size /
                                            (1024 * 1024)
                                        ).toFixed(2)}{" "}
                                        MB
                                    </FileSummaryText>
                                </FileSummary>

                                {/*<H5MetadataPreview
                                    metadata={state.h5Metadata}
                                /> */}

                                <div className="flex justify-center items-center gap-6 mt-4">
                                    <Button
                                        marginTop="mt-1"
                                        onClick={handleUploadClick}
                                    >
                                        Upload
                                    </Button>
                                    <Button
                                        color="red"
                                        size="px-6 py-2.5"
                                        marginTop="mt-1"
                                        onClick={handleClose}
                                    >
                                        Cancel
                                    </Button>
                                </div>
                            </>
                        )}
                    </>
                ) : (
                    <>
                        <DropzoneContainer
                            {...getRootProps()}
                            isDragOver={isDragActive}
                            aria-label={"dropzoneLabel"}
                        >
                            <input {...getInputProps()} />
                            <div className="flex flex-col items-center justify-center space-y-0">
                                <CloudUploadIcon
                                    className="text-gray-400"
                                    style={{ fontSize: "5rem" }}
                                />
                                {isDragActive ? (
                                    <DynamicText
                                        text={"Drop files here..."}
                                        className="text-sm"
                                    />
                                ) : (
                                    <DynamicText
                                        text={
                                            state.selectedFiles.length > 0
                                                ? `Selected file: ${state.selectedFiles[0].name}`
                                                : rejectionMessage
                                        }
                                        className={`${rejectionMessageStyle} text-sm`}
                                    />
                                )}
                                <FileInputLabel
                                    htmlFor="fileInput"
                                    className="text-sm"
                                >
                                    {"Choose File"}
                                </FileInputLabel>
                            </div>
                        </DropzoneContainer>
                    </>
                )}
            </Container>
        );
    });

const Wrapper = (props: FileUploadDialogComponentProps) => {
    const [open, setOpen] = useState(true);

    const handleClose = () => {
        setOpen(false);
    };
    return (
        <ReusableDialog open={open} handleClose={handleClose} component={<FileUploadDialogComponent {...props} />} />
    );
};

export default Wrapper;
