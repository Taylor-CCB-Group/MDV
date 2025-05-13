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

const Loader = () => {
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
    const onDrop = useCallback((acceptedFiles: any[]) => {
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
        maxSize: 2000 * 1024 * 1024, // 2 GB
    });

    const rejectionMessage =
        fileRejections.length > 0
            ? "Only ZIP files are allowed under 200 MB. Please try again."
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

            const form = new FormData();
            form.append("file", file);
            form.append("name", projectName);
            const res = await axios.post("import_project", form);
            if (res.status === 200) {
                // Navigate to the newly created project
                if (res.data?.status === "success") {
                    await fetchProjects();
                    const base = import.meta.env.DEV ? "http://localhost:5170?dir=/" : "";
                    window.location.href = `${base}project/${res.data.id}`;
                } else {
                    throw new Error("An error occurred. Please try again.");
                }
            } else {
                // Error Handling
                if (res.status === 403) {
                    setError({ message: res.data?.error || "Forbidden: You are not allowed to perform this action." });
                } else if (res.status === 500) {
                    setError({
                        message: "Internal Server Error. Please try again.",
                        stack: res.data?.error || "An unknown error occurred",
                    });
                } else {
                    setError({ message: res.data?.error || "An error occurred. Please try again." });
                }
                setErrorOpen(true);
            }
        } catch (error) {
            let err;
            if (error instanceof AxiosError) {
                err = { message: error.response?.data?.error, stack: error?.stack };
            } else if (error instanceof Error) {
                err = { message: error.message, stack: error?.stack };
            } else {
                err = { message: "An error occurred while trying to import the project. Please try again." };
            }
            setError(err);
            setErrorOpen(true);
        } finally {
            setIsLoading(false);
        }
    }, [file, fetchProjects, projectName]);

    return (
        <>
            <Dialog open={open} onClose={onClose} fullWidth>
                <DialogTitle>
                    Import Project
                    <DialogCloseIconButton onClose={onClose} />
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
                                    Import
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
                    <Button onClick={onClose}>Close</Button>
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
