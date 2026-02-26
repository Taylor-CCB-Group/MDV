import useProjects from "@/catalog/hooks/useProjects";
import AlertErrorComponent from "@/charts/dialogs/AlertErrorComponent";
import DebugErrorComponent from "@/charts/dialogs/DebugErrorComponent";
import {
    Container,
    DropzoneContainer,
    DynamicText,
    FileInputLabel,
} from "@/charts/dialogs/FileUploadDialog";
import {
    Box,
    Button,
    CircularProgress,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
} from "@mui/material";
import { CloudUpload as CloudUploadIcon } from "@mui/icons-material";
import axios, { AxiosError, type AxiosProgressEvent } from "axios";
import type React from "react";
import { useCallback, useState } from "react";
import { useDropzone } from "react-dropzone";
import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import { createSocketIOUpload, type SocketIOUploadClient } from "../../charts/dialogs/SocketIOUploadClient";
import { ZipReader, BlobReader } from "@zip.js/zip.js";
import { useApiRoot } from "@/catalog/hooks/useApiRoot";
import ReusableAlertDialog from "@/charts/dialogs/ReusableAlertDialog";

// Constants moved to module level to avoid recreation
const REQUIRED_FILES = new Set(["views.json", "state.json", "datasources.json"]);

// Helper function moved outside component for better organization
const findRootPrefix = (names: string[]): string | null => {
    // Check if all required files are in the root of archive
    const rootFiles = new Set(names.filter(n => !n.includes("/")).map(n => n.split("/").pop() || ""));
    if (REQUIRED_FILES.size === [...REQUIRED_FILES].filter(file => rootFiles.has(file)).length) {
        return "";
    }
    
    // Check one level below if required files exist
    const dirs = new Set(names.filter(n => n.includes("/")).map(n => n.split("/")[0]));
    for (const dir of dirs) {
        const filesInDir = new Set(
            names
                .filter(n => n.startsWith(`${dir}/`))
                .map(n => n.split("/").pop() || "")
        );
        if (REQUIRED_FILES.size === [...REQUIRED_FILES].filter(file => filesInDir.has(file)).length) {
            return `${dir}/`;
        }
    }
    return null;
};

export const Loader = () => {
    return (
        <Box
            sx={{
                display: "flex",
                justifyContent: "center",
                alignItems: "center",
                height: "100%",
                width: "100%",
            }}
        >
            <CircularProgress />
        </Box>
    );
};

const ProgressBar = ({ value, max }: { value: number; max: number }) => (
    <progress
        className="w-full h-8 mb-5 mt-10 bg-gray-200 dark:bg-white-200 border border-gray-300 rounded"
        value={value}
        max={max}
    />
);

export type ImportProjectDialogProps = {
    open: boolean;
    setOpen: React.Dispatch<React.SetStateAction<boolean>>;
};

// Read the environment variable to determine the upload method
const USE_SOCKETIO_UPLOAD = true; // Change to false to use HTTP upload

const ImportProjectDialog = ({ open, setOpen }: ImportProjectDialogProps) => {
    const [file, setFile] = useState<File | null>(null);
    const [projectName, setProjectName] = useState("");
    const [isLoading, setIsLoading] = useState(false);
    const [progress, setProgress] = useState(0);
    const [error, setError] = useState<{ message: string; stack?: string }>();
    const [errorOpen, setErrorOpen] = useState(false);
    const [validationError, setValidationError] = useState<string | null>(null);
    const { fetchProjects } = useProjects();
    
    // State to hold the chosen upload method and the client instance for cancellation
    const [uploadMethod] = useState<'http' | 'socketio'>(
        USE_SOCKETIO_UPLOAD ? 'socketio' : 'http'
    );
    const [uploadClient, setUploadClient] = useState<SocketIOUploadClient | null>(null);

    // Get the mdvApiRoot from the hook
    const mdvApiRoot = useApiRoot();
    
    const validateMDVProject = useCallback(async (file: File): Promise<{ isValid: boolean; error?: string }> => {
        let zipReader: ZipReader<Blob> | null = null;
        
        try {
            zipReader = new ZipReader<Blob>(new BlobReader(file));
            const entries = await zipReader.getEntries();
            
            // Get just the filenames without loading file contents
            const fileNames = entries.map(entry => entry.filename);
            
            // Reject entries with absolute paths or ".."
            const badEntries = fileNames.filter(n => n.startsWith("/") || n.startsWith("\\") || n.includes(".."));
            if (badEntries.length > 0) {
                return { isValid: false, error: "Invalid ZIP file: unsafe paths detected" };
            }

            const root = findRootPrefix(fileNames);
            
            if (root === null) {
                return { isValid: false, error: "Not a valid MDV project: missing required files (views.json, state.json, datasources.json)" };
            }

            return { isValid: true };
        } catch (error) {
            console.error('ZIP validation error:', error);
            return { isValid: false, error: "Failed to read ZIP file contents" };
        } finally {
            if (zipReader) {
                try {
                    await zipReader.close();
                } catch (closeError) {
                    console.warn('Failed to close ZIP reader:', closeError);
                }
            }
        }
    }, []);

    const onDrop = useCallback(async (acceptedFiles: File[]) => {
        if (acceptedFiles[0]) {
            setValidationError(null);
            const file = acceptedFiles[0];
            
            const validation = await validateMDVProject(file);
            if (!validation.isValid) {
                setValidationError(validation.error || "Invalid MDV project file");
                return;
            }
            
            setFile(file);
            const name = file.name.match(/^(.+?)(?:\.mdv)?\.zip$/)?.[1];
            setProjectName(name ?? "");
        }
    }, [validateMDVProject]);

    const { getRootProps, getInputProps, isDragActive, fileRejections } = useDropzone({
        onDrop,
        multiple: false,
        accept: { "application/zip": [".zip"] },
        maxSize: 20000 * 1024 * 1024,
    });

    const rejectionMessage =
        fileRejections.length > 0
            ? "Only ZIP files are allowed under 20 GB. Please try again."
            : validationError
            ? validationError
            : "Drag and drop a file here or click the button below (only MDV project archive *.zip files are allowed)";

    const rejectionMessageStyle = fileRejections.length > 0 || validationError ? "text-red-500" : "";

    const onClose = useCallback(() => {
        if (uploadClient) {
            uploadClient.cancel();
        }
        setOpen(false);
    }, [setOpen, uploadClient]);

    // Original axios logic
    const handleHttpUpload = useCallback(async () => {
        setError(undefined);
        setErrorOpen(false);
        setIsLoading(true);
        setProgress(0);
        try {
            if (!file) {
                setError({ message: "No file selected. Please try again." });
                return;
            }

            const form = new FormData();
            form.append("file", file);

            if (projectName.trim()) {
                form.append("name", projectName);
            }

            // const res = await axios.post("import_project", form);
            const res = await axios.post("import_project", form, {
                onUploadProgress: (progressEvent: AxiosProgressEvent) => {
                    if (progressEvent.total) {
                        const percentCompleted = Math.round(
                            (progressEvent.loaded * 100) / progressEvent.total
                        );
                        setProgress(percentCompleted);
                    }
                },
            });
            if (res.status === 200) {
                if (res.data?.status === "success") {
                    await fetchProjects();
                    const base = import.meta.env.DEV ? "http://localhost:5170?dir=/" : "";
                    window.location.href = `${base}project/${res.data.id}`;
                } else {
                    throw new Error("An error occurred. Please try again.");
                }
            } else {
                throw new Error(res.data?.error || "An error occurred. Please try again.");
            }
        } catch (error) {
            let err;
            if (error instanceof AxiosError) {
                err = { message: error.response?.data?.error, stack: error.response?.status === 500 ? error?.stack : undefined };
            } else if (error instanceof Error) {
                err = { message: error.message, stack: error?.stack };
            } else {
                err = { message: "An error occurred while trying to upload the file. Please try again." };
            }
            setError(err);
            setErrorOpen(true);
        } finally {
            setIsLoading(false);
        }
    }, [file, fetchProjects, projectName]);

    const handleSocketIOUpload = useCallback(async () => {
        if (!file) {
            setError({ message: "No file selected. Please try again." });
            return;
        }
        setIsLoading(true);
        setProgress(0);
        setError(undefined);
        setErrorOpen(false);

        try {
            // Use mdvApiRoot for socket path
            const socketPath = `${mdvApiRoot}socket.io`;
            const client = createSocketIOUpload({
                namespace: "/",
                file: file,
                socketPath, // socket path for path variable: '/test/socket.io'
                onProgress: (percent) => {
                    setProgress(percent);
                },
                onStatusChange: (status, message) => {
                    console.log(`Project Upload Status: ${status}`, message);
                },
                onSuccess: async (result) => {
                    setIsLoading(false);
                    setUploadClient(null);
                    if (result?.result?.project_id) {
                        await fetchProjects();
                        const base = import.meta.env.DEV ? "http://localhost:5170?dir=/" : "";
                        window.location.href = `${base}project/${result.result.project_id}`;
                    } else {
                        throw new Error("Upload succeeded but did not return a project ID.");
                    }
                },
                onError: (error) => {
                    setIsLoading(false);
                    setUploadClient(null);
                    const err = { message: error.message || "An error occurred during upload.", stack: error.stack };
                    setError(err);
                    setErrorOpen(true);
                },
            });
            setUploadClient(client);
            await client.upload();

        } catch (error) {
            setIsLoading(false);
            setUploadClient(null);
            let err;
            if (error instanceof Error) {
                err = { message: error.message, stack: error?.stack };
            } else {
                err = { message: "An unexpected error occurred. Please try again." };
            }
            setError(err);
            setErrorOpen(true);
        }
    }, [file, fetchProjects, mdvApiRoot]);

    // Main handler to delegate to the correct upload function
    const handleUpload = useCallback(() => {
        console.log(`Using ${uploadMethod} upload method for project import.`);
        if (uploadMethod === 'socketio') {
            handleSocketIOUpload();
        } else {
            handleHttpUpload();
        }
    }, [uploadMethod, handleSocketIOUpload, handleHttpUpload]);

    return (
        <>
            <Dialog open={open} onClose={onClose} fullWidth>
                <DialogTitle>
                    Import Project
                    <DialogCloseIconButton onClose={onClose} />
                </DialogTitle>
                <DialogContent dividers>
                    {isLoading ? (
                        <div className="flex flex-col justify-center items-center w-full h-40 dark:bg-[#333]">
                            <p className="text-lg font-bold text-[#333] dark:text-white text-center">
                                Your project is being uploaded, please wait...
                            </p>
                            <ProgressBar value={progress} max={100} />
                        </div>
                    ) : (
                        <Container>
                            <DropzoneContainer
                                {...getRootProps()}
                                isDragOver={isDragActive}
                                aria-label={"dropzoneLabel"}
                                style={{ border: "none" }}
                            >
                                {!file && <input {...getInputProps()} />}
                                <div className="flex flex-col items-center justify-center space-y-0">
                                    <CloudUploadIcon className="text-gray-400" style={{ fontSize: "5rem" }} />
                                    {isDragActive ? (
                                        <DynamicText text={"Drop files here..."} className="text-sm" />
                                    ) : (
                                        <DynamicText
                                            text={file ? `Selected file: ${file.name} (${projectName})` : rejectionMessage}
                                            className={`${rejectionMessageStyle} text-sm`}
                                        />
                                    )}
                                    {!file && (
                                        <FileInputLabel htmlFor="fileInput" className="text-sm">
                                            Choose File
                                        </FileInputLabel>
                                    )}
                                </div>
                            </DropzoneContainer>
                            {file ? (
                                <Button onClick={handleUpload} sx={{ mt: 2, minWidth: "50%" }} variant="contained">
                                    Upload file
                                </Button>
                            ) : (
                                <></>
                            )}
                        </Container>
                    )}
                </DialogContent>
                <DialogActions>
                    <Button onClick={onClose}>Close</Button>
                </DialogActions>
            </Dialog>
            {error && (
                <ReusableAlertDialog
                    open={errorOpen}
                    handleClose={() => setErrorOpen(false)}
                    isAlertErrorComponent={!error?.stack}
                    component={
                        error?.stack ? (
                            <DebugErrorComponent
                                error={{ message: error.message, stack: error?.stack }}
                                extraMetadata={{ error: error.message }}
                            />
                        ) : (
                            <AlertErrorComponent message={error.message} />
                        )
                    }
                />
            )}
        </>
    );
};

export default ImportProjectDialog;