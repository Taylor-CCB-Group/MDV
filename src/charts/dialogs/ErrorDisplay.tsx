import {
    Check,
    ContentCopy,
    Error as ErrorIcon,
    ExpandLess,
    ExpandMore,
} from "@mui/icons-material";
import {
    Alert,
    AlertTitle,
    Collapse,
    IconButton,
    Paper,
    Snackbar,
    Tooltip,
} from "@mui/material";
import React, { useState } from "react";

interface ErrorDisplayProps {
    error: {
        message: string;
        traceback?: string;
    };
    title?: string;
}

const ErrorDisplay = ({
    error,
    title = "Error Occurred",
}: ErrorDisplayProps) => {
    const [expanded, setExpanded] = useState(false);
    const [copied, setCopied] = useState(false);

    const handleCopy = async () => {
        try {
            await navigator.clipboard.writeText(
                `Error: ${error.message}\n\nTraceback:\n${error.traceback}`,
            );
            setCopied(true);
            setTimeout(() => setCopied(false), 2000);
        } catch (err) {
            console.error("Failed to copy error details:", err);
        }
    };

    return (
        <div style={{ maxWidth: 800, margin: "20px auto", width: "90%" }}>
            <Paper
                elevation={3}
                sx={{
                    bgcolor: "var(--background_color_error)",
                    color: "var(--text_color)",
                }}
            >
                <Alert
                    severity="error"
                    icon={<ErrorIcon />}
                    sx={{
                        "& .MuiAlert-message": {
                            width: "100%",
                        },
                        bgcolor: "transparent",
                        color: "var(--text_color)",
                        "& .MuiAlert-icon": {
                            color: "var(--icon_color_error)",
                        },
                    }}
                >
                    <AlertTitle
                        sx={{
                            fontSize: "1.1rem",
                            fontWeight: "bold",
                            display: "flex",
                            alignItems: "center",
                            justifyContent: "space-between",
                            color: "var(--text_color_error)",
                        }}
                    >
                        {title}
                        {error.traceback && (
                            <Tooltip title="Copy error details">
                                <IconButton
                                    size="small"
                                    onClick={handleCopy}
                                    sx={{
                                        ml: 1,
                                        color: "var(--text_color_error)",
                                        "&:hover": {
                                            backgroundColor:
                                                "var(--background_color)",
                                        },
                                    }}
                                >
                                    {copied ? <Check /> : <ContentCopy />}
                                </IconButton>
                            </Tooltip>
                        )}
                    </AlertTitle>

                    <div style={{ marginTop: "12px" }}>
                        <p
                            style={{
                                margin: "0 0 12px 0",
                                lineHeight: "1.5",
                                color: "var(--text_color_error)",
                                whiteSpace: "pre-wrap",
                                wordBreak: "break-word",
                            }}
                        >
                            {error.message}
                        </p>

                        {error.traceback && (
                            <div>
                                <IconButton
                                    size="small"
                                    onClick={() => setExpanded(!expanded)}
                                    sx={{
                                        mb: 1,
                                        fontSize: "0.875rem",
                                        display: "flex",
                                        alignItems: "center",
                                        gap: "4px",
                                        color: "var(--icon_color_error)",
                                        "&:hover": {
                                            backgroundColor:
                                                "var(--background_color)",
                                        },
                                    }}
                                >
                                    {expanded ? <ExpandLess /> : <ExpandMore />}
                                    {expanded ? "Hide Details" : "Show Details"}
                                </IconButton>

                                <Collapse in={expanded}>
                                    <pre
                                        style={{
                                            backgroundColor:
                                                "var(--background_color)",
                                            padding: "16px",
                                            borderRadius: "4px",
                                            margin: 0,
                                            maxHeight: "400px",
                                            overflowY: "auto",
                                            overflowX: "hidden",
                                            fontSize: "0.875rem",
                                            lineHeight: 1.5,
                                            fontFamily:
                                                "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
                                            border: "1px solid var(--input_border_color)",
                                            whiteSpace: "pre-wrap",
                                            wordBreak: "break-word",
                                            color: "var(--text_color)",
                                        }}
                                    >
                                        {error.traceback}
                                    </pre>
                                </Collapse>
                            </div>
                        )}
                    </div>
                </Alert>
            </Paper>

            <Snackbar
                open={copied}
                autoHideDuration={2000}
                onClose={() => setCopied(false)}
                message="Error details copied to clipboard"
                anchorOrigin={{ vertical: "bottom", horizontal: "center" }}
            />
        </div>
    );
};

export default ErrorDisplay;
