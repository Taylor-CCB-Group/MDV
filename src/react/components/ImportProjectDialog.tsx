import useProjects from "@/catalog/hooks/useProjects";
import AlertErrorComponent from "@/charts/dialogs/AlertErrorComponent";
import DebugErrorComponent from "@/charts/dialogs/DebugErrorComponent";
import {
    Container,
    DropzoneContainer,
    DynamicText,
    FileInputLabel,
    Button as FileUploadButton,
} from "@/charts/dialogs/FileUploadDialog";
import ReusableDialog from "@/charts/dialogs/ReusableDialog";
import {
    Box,
    Button,
    CircularProgress,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    TextField,
    Typography,
} from "@mui/material";
import { CloudUpload as CloudUploadIcon } from "@mui/icons-material";
import axios, { AxiosError } from "axios";
import type React from "react";
import { useCallback, useState } from "react";
import { useDropzone } from "react-dropzone";
import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";

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

export type ImportProjectDialogProps = {
    open: boolean;
    setOpen: React.Dispatch<React.SetStateAction<boolean>>;
};

const ImportProjectDialog = ({ open, setOpen }: ImportProjectDialogProps) => {
    const [file, setFile] = useState<File | null>(null);
    const [projectName, setProjectName] = useState("");
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<{ message: string; stack?: string }>();
    const [errorOpen, setErrorOpen] = useState(false);
    const { fetchProjects } = useProjects();

    // todo: add additional checks
    const onDrop = useCallback((acceptedFiles: File[]) => {
        if (acceptedFiles[0]) {
            setFile(acceptedFiles[0]);
            const name = acceptedFiles[0].name.split(".")[0];
            setProjectName(name ?? "");
        }
    }, []);

    const { getRootProps, getInputProps, isDragActive, fileRejections } = useDropzone({
        onDrop,
        multiple: false,
        accept: {
            "application/zip": [".zip"],
        },
        // todo: update later for larger files support
        maxSize: 2000 * 1024 * 1024, // 2 GB
    });

    const rejectionMessage =
        fileRejections.length > 0
            ? "Only ZIP files are allowed under 2 GB. Please try again."
            : "Drag and drop a file here or click the button below (only *.zip files are allowed)";

    const rejectionMessageStyle = fileRejections.length > 0 ? "text-red-500" : "";

    const onClose = useCallback(() => {
        setOpen(false);
    }, [setOpen]);

    const onUpload = useCallback(async () => {
        setError(undefined);
        setErrorOpen(false);
        setIsLoading(true);
        try {
            if (!file) {
                setError({ message: "No file selected. Please try again." });
                return;
            }

            // Append file to form data
            const form = new FormData();
            form.append("file", file);

            // Append name to form data if the name exists
            if (projectName.trim()) {
                form.append("name", projectName);
            }

            const res = await axios.post("import_project", form);
            if (res.status === 200) {
                // Navigate to the newly created project and fetch the projects
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
            // Handling errors based on the error type and status code
            if (error instanceof AxiosError) {
                err = {
                    message: error.response?.data?.error,
                    stack: error.response?.status === 500 ? error?.stack : undefined,
                };
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

    return (
        <>
            <Dialog open={open} onClose={!isLoading ? onClose : undefined} fullWidth>
                <DialogTitle>
                    Import Project
                    <DialogCloseIconButton onClose={onClose} disabled={isLoading} />
                </DialogTitle>
                <DialogContent dividers>
                    {!isLoading ? (
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
                                            text={file ? `Selected file: ${file.name}` : rejectionMessage}
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
                                <Button onClick={onUpload} sx={{ mt: 2, minWidth: "50%" }} variant="contained">
                                    Upload file
                                </Button>
                            ) : (
                                <></>
                            )}
                        </Container>
                    ) : (
                        <Loader />
                    )}
                </DialogContent>
                <DialogActions>
                    <Button onClick={onClose} disabled={isLoading}>Close</Button>
                </DialogActions>
            </Dialog>
            {error && (
                <ReusableDialog
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
