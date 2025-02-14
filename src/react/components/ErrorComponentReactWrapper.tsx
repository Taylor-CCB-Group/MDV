import ErrorDisplay, { type ErrorMetadata } from "@/charts/dialogs/ErrorDisplay";
import ReusableDialog from "@/charts/dialogs/ReusableDialog";
import { Error as ErrorIcon } from "@mui/icons-material";
import { Alert, Button, Container, Paper } from "@mui/material";
import { useState } from "react";

export type ErrorComponentType = {
    error: {
        message: string,
        stack?: string
    },
    height?: number,
    width?: number,
    extraMetaData?: ErrorMetadata,
};

const ErrorComponent = ({error, extraMetaData}: ErrorComponentType) => {
    const [open, setOpen] = useState(false);
    const handleClose = () => setOpen(false);
    const handleOpen = () => setOpen(true);
    return (
        <Container>
            <Paper
                elevation={24}
                sx={{
                    bgcolor: "var(--background_color_error)",
                    color: "var(--text_color)",
                }}
            >
                <Alert
                    severity="error"
                    icon={
                        <ErrorIcon
                            sx={{
                                fontSize: 25,
                                color: "red",
                                marginTop: "1px",
                            }}
                        />
                    }
                    sx={{
                        "& .MuiAlert-message": {
                            width: "100%",
                        },
                        bgcolor: "transparent",
                        color: "var(--text_color)",
                        "& .MuiAlert-icon": {
                            color: "var(--icon_color_error)",
                        },
                        alignItems: "center",
                    }}
                >
                    <Button
                        onClick={handleOpen}
                        variant="text"
                        sx={{
                            fontSize: "0.9rem",
                            textTransform: "none",
                            color: "var(--text_color_error)",
                            ":hover": {
                                background: "none",
                            },
                        }}
                    >
                        ERROR: Click to view details
                    </Button>
                </Alert>
                <ReusableDialog
                    open={open}
                    handleClose={handleClose}
                    component={
                        <ErrorDisplay
                            error={{
                                message: error.message,
                                traceback: error.stack,
                            }}
                            extraMetadata={extraMetaData}
                        />
                    }
                />
            </Paper>
        </Container>
    );
};



const ErrorComponentReactWrapper = ({error, height, width, extraMetaData}: ErrorComponentType) => {
    return (
        <ErrorComponent error={error} height={height} width={width} extraMetaData={extraMetaData} />
    );
};

export default ErrorComponentReactWrapper;