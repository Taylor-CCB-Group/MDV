import { Check, ContentCopy } from "@mui/icons-material";
import { Alert, AlertTitle, Box, Divider, IconButton, Link, Paper, Tooltip, Typography } from "@mui/material";
import { useState } from "react";
import { ALERT_ERROR_REPORT_MESSAGE, MDV_EMAIL } from "@/utilities/constants";

export type AlertType = "error" | "info" | "success" | "warning";

export type AlertErrorComponentProps = {
    message: string;
    title?: string;
    alertType?: AlertType;
};

const AlertErrorComponent = ({ message, title, alertType = "error" }: AlertErrorComponentProps) => {
    const [copied, setCopied] = useState(false);

    const handleCopy = async () => {
        try {
            await navigator.clipboard.writeText(message);
            setCopied(true);
            setTimeout(() => setCopied(false), 2000);
        } catch (err) {
            console.error("Failed to copy error message:", err);
        }
    };

    const backgroundColor = `var(--background_color_${alertType})`;
    const textColor = `var(--text_color_${alertType})`;

    return (
        <div
            style={{
                margin: "30px",
                maxWidth: "1000px",
                minHeight: "15vh",
                display: "flex",
                flexDirection: "column",
                justifyContent: "center",
                alignItems: "center",
            }}
        >
            <Paper
                elevation={3}
                sx={{
                    bgcolor: backgroundColor,
                    color: "var(--text_color)",
                    height: "100%",
                    width: "100%",
                }}
            >
                <>
                    <Alert
                        severity={alertType}
                        sx={{
                            "& .MuiAlert-message": {
                                width: "100%",
                            },
                            bgcolor: "transparent",
                            color: textColor,
                            padding: 2,
                            display: "flex",
                            alignItems: title ? "flex-start" : "center",
                            borderRadius: 0,
                            boxShadow: "none",
                        }}
                    >
                        {title && <AlertTitle variant="h6">{title}</AlertTitle>}
                        <Typography
                            variant="subtitle1"
                            sx={{
                                display: "flex",
                                justifyContent: "space-between",
                                alignItems: "center",
                                gap: 1,
                            }}
                        >
                            <span style={{ flex: 1, minWidth: 0, wordBreak: "break-word" }}>
                                {message}
                            </span>
                            {alertType === "error" && (
                                <Tooltip title="Copy error message">
                                    <IconButton
                                        size="small"
                                        onClick={handleCopy}
                                        sx={{
                                            ml: 1,
                                            color: "var(--text_color_error)",
                                        }}
                                    >
                                        {copied ? <Check /> : <ContentCopy />}
                                    </IconButton>
                                </Tooltip>
                            )}
                        </Typography>
                    </Alert>
                    <Divider />
                    {alertType === "error" ? (
                        <Box sx={{ p: 2, mt: 1 }}>
                            <Typography>
                                {ALERT_ERROR_REPORT_MESSAGE}{" "}
                                <Link 
                                    href={`mailto:${MDV_EMAIL}?subject=MDV%20Issue%20Report`}
                                    sx={{ 
                                        textDecoration: "none",
                                        "&.MuiLink-root": { color: "info.main" }
                                    }}
                                >
                                    {MDV_EMAIL}
                                </Link>
                            </Typography>
                        </Box>
                    ) : (
                        <></>
                    )}
                </>
            </Paper>
        </div>
    );
};

export default AlertErrorComponent;
