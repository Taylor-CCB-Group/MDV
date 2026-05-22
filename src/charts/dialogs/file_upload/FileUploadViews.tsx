import type { ReactNode } from "react";
import CloudUploadIcon from "@mui/icons-material/CloudUpload";

import DebugErrorComponent from "../DebugErrorComponent";
import type {
    DatasourceSummary,
    TiffSummary,
    UploadErrorPayload,
} from "@/utils/file_upload/fileUploadTypes";
import {
    Button,
    DatasourceNameInput,
    DynamicText,
    FileInputLabel,
    FileSummary,
    FileSummaryHeading,
    FileSummaryText,
    ProgressBar,
    Message,
    Spinner,
    StatusContainer,
    SuccessContainer,
    SuccessHeading,
    SuccessText,
    DropzoneContainer,
} from "./FileUploadUI";

export function RedirectingView() {
    return (
        <StatusContainer>
            <Message>Redirecting to your new view…</Message>
            <Spinner />
        </StatusContainer>
    );
}

export function UploadingView({ progress }: { progress: number }) {
    return (
        <StatusContainer>
            <Message>
                {"Your file is being uploaded, please wait..."}
            </Message>
            <ProgressBar value={progress} max="100" />
        </StatusContainer>
    );
}

export function ProcessingView() {
    return (
        <StatusContainer>
            <Message>
                {"Your file is being processed, please wait..."}
            </Message>
            <Spinner />
        </StatusContainer>
    );
}

export function UploadSuccessView({ onRefresh }: { onRefresh: () => void }) {
    return (
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
                onClick={onRefresh}
            >
                Refresh Page
            </Button>
        </>
    );
}

export function UploadErrorView({
    error,
    onUploadAnother,
    onCancel,
}: {
    error: UploadErrorPayload;
    onUploadAnother: () => void;
    onCancel: () => void;
}) {
    return (
        <>
            <DebugErrorComponent error={error} />
            <div className="flex justify-center items-center gap-6">
                <Button
                    marginTop="mt-1"
                    onClick={onUploadAnother}
                >
                    Upload Another File
                </Button>
                <Button
                    color="red"
                    size="px-14 py-2.5"
                    marginTop="mt-1"
                    onClick={onCancel}
                >
                    Cancel
                </Button>
            </div>
        </>
    );
}

export function ValidatingView() {
    return (
        <StatusContainer>
            <Message>{"Validating data, please wait..."}</Message>
            <Spinner />
        </StatusContainer>
    );
}

export function TabularReviewView({
    datasourceName,
    datasourceSummary,
    onDatasourceNameChange,
    isDatasourceDisabled,
    onUpload,
    onCancel,
    columnPreview,
}: {
    datasourceName: string;
    datasourceSummary: DatasourceSummary;
    onDatasourceNameChange: (event: unknown) => void;
    isDatasourceDisabled: boolean;
    onUpload: () => void;
    onCancel: () => void;
    columnPreview?: ReactNode;
}) {
    return (
        <>
            <FileSummary>
                <FileSummaryHeading>
                    {"Uploaded File Summary"}
                </FileSummaryHeading>
                <FileSummaryText>
                    <DatasourceNameInput
                        value={datasourceName}
                        onChange={onDatasourceNameChange}
                        isDisabled={isDatasourceDisabled}
                    />
                </FileSummaryText>
                <FileSummaryText>
                    <strong>{"File name"}</strong>{" "}
                    {datasourceSummary.fileName}
                </FileSummaryText>
                {/* <FileSummaryText>
                    <strong>{"Number of rows"}</strong>{" "}
                    {datasourceSummary.rowCount}
                </FileSummaryText>
                <FileSummaryText>
                    <strong>{"Number of columns"}</strong>{" "}
                    {datasourceSummary.columnCount}
                </FileSummaryText> */}
                <FileSummaryText>
                    <strong>{"File size"}</strong>{" "}
                    {datasourceSummary.fileSize} MB
                </FileSummaryText>
            </FileSummary>
            {columnPreview}
            <div className="flex justify-center items-center gap-6 mt-4">
                <Button
                    marginTop="mt-1"
                    onClick={onUpload}
                >
                    {"Upload"}
                </Button>
                <Button
                    color="red"
                    size="px-6 py-2.5"
                    marginTop="mt-1"
                    onClick={onCancel}
                >
                    {"Cancel"}
                </Button>
            </div>
        </>
    );
}

export function TiffReviewView({
    tiffSummary,
    tiffMetadata,
    selectedFile,
    showMetadata,
    onToggleView,
    onUpload,
    onCancel,
    datasourceDropdown,
}: {
    tiffSummary: TiffSummary;
    tiffMetadata: unknown;
    selectedFile: File | undefined;
    showMetadata: boolean;
    onToggleView: () => void;
    onUpload: () => void;
    onCancel: () => void;
    datasourceDropdown?: ReactNode;
}) {
    return (
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
                            {datasourceDropdown}
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
                        {/* <TiffVisualization
                            metadata={
                                tiffMetadata
                            }
                            file={
                                selectedFile
                            }
                        /> */}
                    </div>
                </div>
            </div>
            <div className="w-2/3">
                <div className="flex items-start justify-start h-[70%] min-h-[420px]">
                    {/* {showMetadata ? (
                        <TiffMetadataTable
                            metadata={
                                tiffMetadata
                            }
                        />
                    ) : (
                        <TiffPreview
                            metadata={
                                tiffMetadata
                            }
                        />
                    )} */}
                </div>
                <div className="flex justify-between items-center ml-4 mr-2 mt-12">
                    <Button
                        color="gray"
                        size="px-6 py-2.5"
                        marginTop="mt-1"
                        onClick={onToggleView}
                    >
                        {showMetadata
                            ? "View Metadata as XML"
                            : "View Metadata as Table"}
                    </Button>
                    <div className="flex space-x-4 ml-auto">
                        <Button
                            marginTop="mt-6"
                            onClick={onUpload}
                        >
                            {"Upload"}
                        </Button>
                        <Button
                            color="red"
                            size="px-6 py-2.5"
                            marginTop="mt-6"
                            onClick={onCancel}
                        >
                            {"Cancel"}
                        </Button>
                    </div>
                </div>
            </div>
        </div>
    );
}

export function H5ReviewView({
    selectedFile,
    h5Metadata,
    onUpload,
    onCancel,
}: {
    selectedFile: File | undefined;
    h5Metadata: unknown;
    onUpload: () => void;
    onCancel: () => void;
}) {
    if (!selectedFile) {
        return null;
    }

    return (
        <>
            <FileSummary>
                <FileSummaryHeading>
                    H5 File Summary
                </FileSummaryHeading>
                <FileSummaryText>
                    <strong>File name:</strong>{" "}
                    {selectedFile.name}
                </FileSummaryText>
                <FileSummaryText>
                    <strong>File size:</strong>{" "}
                    {(
                        selectedFile.size /
                        (1024 * 1024)
                    ).toFixed(2)}{" "}
                    MB
                </FileSummaryText>
            </FileSummary>

            {/*<H5MetadataPreview
                metadata={h5Metadata}
            /> */}

            <div className="flex justify-center items-center gap-6 mt-4">
                <Button
                    marginTop="mt-1"
                    onClick={onUpload}
                >
                    Upload
                </Button>
                <Button
                    color="red"
                    size="px-6 py-2.5"
                    marginTop="mt-1"
                    onClick={onCancel}
                >
                    Cancel
                </Button>
            </div>
        </>
    );
}

export function FileDropzoneView({
    rootProps,
    inputProps,
    isDragActive,
    selectedFileName,
    rejectionMessage,
    rejectionMessageStyle,
    supportedFileTypes,
}: {
    rootProps: Record<string, unknown>;
    inputProps: Record<string, unknown>;
    isDragActive: boolean;
    selectedFileName: string | null;
    rejectionMessage: string;
    rejectionMessageStyle: string;
    supportedFileTypes: string[];
}) {
    return (
        <DropzoneContainer
            {...rootProps}
            isDragOver={isDragActive}
            aria-label={"dropzoneLabel"}
        >
            <input {...inputProps} />
            <div className="flex flex-col items-center justify-center">
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
                            selectedFileName
                                ? `Selected file: ${selectedFileName}`
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
                <div className="mt-6 w-full max-w-md rounded-md bg-gray-50 px-3 py-2 dark:bg-gray-900/50">
                    <p className="mb-3 text-center text-xs font-semibold uppercase tracking-wide text-gray-500 dark:text-gray-400">
                        Supported file types
                    </p>
                    <div className="flex flex-wrap items-center justify-center gap-2">
                    {supportedFileTypes.map((extension) => (
                        <span
                            key={extension}
                            className="rounded-full border border-gray-200 bg-white px-2.5 py-1 text-xs font-medium text-gray-700 shadow-sm dark:border-gray-600 dark:bg-gray-800 dark:text-gray-200"
                        >
                            {extension}
                        </span>
                    ))}
                    </div>
                </div>
            </div>
        </DropzoneContainer>
    );
}
