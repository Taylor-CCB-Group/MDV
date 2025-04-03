import useBuildInfo from "@/catalog/hooks/useBuildInfo";
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
    Box,
    Button,
    CircularProgress,
    Collapse,
    Container,
    Divider,
    IconButton,
    Paper,
    Snackbar,
    TextareaAutosize,
    Tooltip,
    Typography,
} from "@mui/material";
import axios from "axios";
import { useEffect, useState } from "react";
import JsonView from "react18-json-view";

/** Any extra information we may want to use
 *
 * @remarks Only an object is expected, large objects like data store will break the ErrorDisplay component
 */
//! Don't send large objects like DataStore which will break the ErrorDisplay component, also look for an alternative way to handle this issue
export type ErrorMetadata = object;
export interface ErrorDisplayProps {
    error: {
        message: string;
        traceback?: string;
        stack?: string;
    };
    title?: string;
    extraMetadata?: ErrorMetadata;
}

const ErrorDisplay = ({
    error,
    title = "Error Occurred",
    extraMetadata,
}: ErrorDisplayProps) => {
    const [expanded, setExpanded] = useState(true);
    const [copied, setCopied] = useState(false);
    const [userComments, setUserComments] = useState<string>();
    const info = useBuildInfo();
    const buildInfo = info?.buildInfo || {};
    const [isLoading, setIsLoading] = useState(false);
    const [emailSent, setEmailSent] = useState(false);
    const [emailNotSent, setEmailNotSent] = useState(false);
    const [metaData, setMetaData] = useState<Record<string, unknown> | null>(null);

    useEffect(() => {
        if (extraMetadata || error?.stack || error?.traceback)
            setMetaData({
                ...extraMetadata,
                // Create the key traceback or stack only if the value exists
                ...(error?.traceback != null && {
                    traceback: error?.traceback,
                }),
                ...(error?.stack != null && { stack: error?.stack }),
            });
    }, [extraMetadata, error]);

    // todo: Uncomment later
    //todo: Add the logic to send the error details to the email address and display the corresponding message to user
    // Send the error details and the user's comments (if any) to the support email address
    // const handleSend = async () => {
    //     setIsLoading(true);
    //     setEmailSent(false);
    //     setEmailNotSent(false);
    //     const errorDetails = {
    //         message: error.message,
    //         userComments: userComments ? userComments : null,
    //         extraMetadata: metaData ? metaData : null,
    //         // buildInfo
    //     };
    //     try {
    //         const res = await axios.post(
    //             "/send-error",
    //             JSON.stringify(errorDetails),
    //         );
    //         const data = res.data;
    //         setEmailSent(true);
    //     } catch (e) {
    //         console.log("error", e);
    //         setEmailNotSent(true);
    //     } finally {
    //         setIsLoading(false);
    //     }
    //     console.log("Send", errorDetails);
    // };

    const handleCopy = async () => {
        try {
            await navigator.clipboard.writeText(
                `Error: ${error.message}${metaData ? `\n\nMetaData:\n${JSON.stringify(metaData, null, 2)}` : ""}`,
            );
            setCopied(true);
            setTimeout(() => setCopied(false), 2000);
        } catch (err) {
            console.error("Failed to copy error details:", err);
        }
    };

    return (
        <div
            style={{
                // maxWidth: 800,
                // minWidth: 500,
                margin: "20px auto",
                // width: "90%",
                width: "45vw",
                maxWidth: "1000px",
                minHeight: "45vh",
                display: "flex",
                justifyContent: "center",
                alignItems: "center",
            }}
        >
            <Paper
                elevation={3}
                sx={{
                    bgcolor: "var(--background_color_error)",
                    color: "var(--text_color)",
                    height: "95%",
                    width: "95%",
                }}
            >
                {isLoading ? (
                    // Display loading icon when button is clicked
                    <Box
                        sx={{
                            height: "40vh",
                            width: "42vw",
                            display: "flex",
                            alignItems: "center",
                            justifyContent: "center",
                        }}
                    >
                        <CircularProgress color="inherit" />
                    </Box>
                ) : emailSent ? (
                    // Display success message when email is sent
                    <Box
                        sx={{
                            height: "40vh",
                            width: "42vw",
                            display: "flex",
                            alignItems: "center",
                            justifyContent: "center",
                        }}
                    >
                        <Typography>
                            Thank you for your feedback, we will do our best to
                            fix this issue as soon as possible.
                        </Typography>
                    </Box>
                ) : (
                    <>
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
                                {metaData && (
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
                                            {copied ? (
                                                <Check />
                                            ) : (
                                                <ContentCopy />
                                            )}
                                        </IconButton>
                                    </Tooltip>
                                )}
                            </AlertTitle>

                            <div style={{ marginTop: "12px" }}>
                                <p
                                    style={{
                                        margin: "0 0 12px 0",
                                        lineHeight: "1.6",
                                        padding: "4px 8px",
                                        color: "var(--text_color_error)",
                                        whiteSpace: "pre-wrap",
                                        wordBreak: "break-word",
                                    }}
                                >
                                    {error.message}
                                </p>

                                {metaData && (
                                    <div>
                                        <IconButton
                                            size="medium"
                                            onClick={() =>
                                                setExpanded(!expanded)
                                            }
                                            sx={{
                                                mb: 1,
                                                fontSize: "0.875rem",
                                                display: "flex",
                                                alignItems: "center",
                                                justifyContent: "center",
                                                padding: "6px 12px",
                                                borderRadius: "6px",
                                                color: "var(--icon_color_error)",
                                                "&:hover": {
                                                    backgroundColor:
                                                        "var(--background_color_hover)",
                                                },
                                            }}
                                        >
                                            {expanded ? (
                                                <ExpandLess />
                                            ) : (
                                                <ExpandMore />
                                            )}
                                            <span style={{ marginLeft: "6px" }}>
                                                {expanded
                                                    ? "Hide Details"
                                                    : "Show Details"}
                                            </span>
                                        </IconButton>

                                        <Collapse in={expanded}>
                                            <JsonView
                                                editable
                                                src={{
                                                    ...metaData, 
                                                    // buildInfo: buildInfo
                                                }}
                                                onAdd={(params) =>
                                                    setMetaData(params.src)
                                                }
                                                onDelete={(params) =>
                                                    setMetaData(params.src)
                                                }
                                                onEdit={(params) =>
                                                    setMetaData(params.src)
                                                }
                                                collapsed={1}
                                                customizeNode={(_params) => {
                                                    return { add: false };
                                                }}
                                            />
                                        </Collapse>
                                    </div>
                                )}
                            </div>
                        </Alert>
                        <Divider />
                        <Container
                            sx={{
                                pt: 2,
                            }}
                        >
                            <Typography variant="h6" sx={{ mb: 2 }}>
                                We are sorry for the inconvenience.
                            </Typography>
                            <Typography sx={{ mb: 2 }}>
                                To help us diagnose and fix the problem, please
                                copy the error information on the top right corner and send us.
                            </Typography>
                            {metaData && (
                                <Typography sx={{ mb: 2 }}>
                                    Please remove any sensitive information from
                                    the error details above before sending.
                                </Typography>
                            )}
                            {/* <TextareaAutosize
                                minRows={3}
                                style={{
                                    width: "100%",
                                    padding: "10px 10px",
                                    marginTop: 2,
                                    backgroundColor:
                                        "var(--input_background_color)",
                                    borderWidth: 2,
                                    borderColor: "var(--input_border_color)",
                                }}
                                placeholder="Please provide any additional comments (optional)"
                                value={userComments}
                                onChange={(e) =>
                                    setUserComments(e.target.value)
                                }
                            />
                            {emailNotSent && (
                                <Typography
                                    sx={{
                                        color: "var(--icon_color_error)",
                                        mt: 1,
                                    }}
                                >
                                    Something went wrong, please try again.
                                </Typography>
                            )}
                            <Button
                                variant="contained"
                                color="primary"
                                sx={{ mt: 2, mb: 3, justifySelf: "flex-start" }}
                                onClick={handleSend}
                            >
                                Send
                            </Button> */}
                        </Container>
                    </>
                )}
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
