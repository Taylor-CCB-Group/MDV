import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    IconButton,
    Paper,
    type PaperProps,
    TextField,
    Typography,
} from "@mui/material";
import { Close as CloseIcon } from "@mui/icons-material";
import { useCallback, useRef, useState } from "react";
import Draggable from "react-draggable";

export type FindAndReplaceDialogProps = {
    open: boolean;
    onClose: () => void;
    handleFind: (find: string) => void;
    handleFindPrev: () => void;
    handleFindNext: () => void;
    handleReplace: (find: string, replace: string) => void;
    handleReplaceAll: (find: string, replace: string) => void;
    foundMatches?: number | null;
    columnName?: string | null; // The column being searched/replaced
    disableFindNext?: boolean;
    disableFindPrev?: boolean;
    isColumnEditable?: boolean;
};

const PaperComponent = (props: PaperProps) => {
    const nodeRef = useRef<HTMLDivElement>(null);
    return (
        <Draggable
            nodeRef={nodeRef as React.RefObject<HTMLDivElement>}
            handle="#draggable-dialog-title"
            cancel={'[class*="MuiDialogContent-root"]'}
        >
            <Paper {...props} ref={nodeRef} />
        </Draggable>
    );
};

const FindAndReplaceDialog = ({
    open,
    onClose,
    handleFind,
    handleFindNext,
    handleFindPrev,
    handleReplace,
    handleReplaceAll,
    foundMatches,
    columnName,
    disableFindNext = false,
    disableFindPrev = false,
    isColumnEditable = true,
}: FindAndReplaceDialogProps) => {
    const [findText, setFindText] = useState("");
    const [replaceText, setReplaceText] = useState("");

    // Clear the find and replace text before closing
    const onDialogClose = useCallback(() => {
        setFindText("");
        setReplaceText("");
        onClose();
    }, [onClose]);

    return (
        <Dialog
            open={open}
            onClose={onDialogClose}
            PaperComponent={PaperComponent}
            disableAutoFocus
            disableEnforceFocus
            disableRestoreFocus
            fullWidth
            maxWidth="xs"
            hideBackdrop
            PaperProps={{
                style: {
                    boxShadow: "0 8px 32px rgba(0, 0 , 0, 0.6)",
                },
            }}
        >
            <DialogTitle style={{ cursor: "move" }} id="draggable-dialog-title">
                Find And Replace {columnName ? `in "${columnName}"` : ""}
            </DialogTitle>
            <IconButton
                aria-label="close"
                onClick={onDialogClose}
                sx={{
                    position: "absolute",
                    right: 8,
                    top: 8,
                    color: (theme) => theme.palette.grey[500],
                }}
            >
                <CloseIcon />
            </IconButton>
            <DialogContent dividers>
                <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
                    {/* Find Section */}
                    <Box sx={{ display: "flex", gap: 1, alignItems: "center" }}>
                        <TextField
                            label="Find Text"
                            variant="outlined"
                            size="small"
                            value={findText}
                            onChange={(e) => setFindText(e.target.value)}
                            onKeyDown={(e) => e.key === "Enter" && handleFind(findText)}
                            sx={{ flexGrow: 1 }}
                        />
                    </Box>
                    {/* foundMatches || typeof foundMatches === "number" - This condition is for zero */}
                    {findText && (foundMatches || typeof foundMatches === "number") ? (
                        <Typography>{`Found ${foundMatches} matches`}</Typography>
                    ) : (
                        <></>
                    )}

                    {/* Navigation Section */}
                    <Box sx={{ display: "flex", gap: 3 }}>
                        <Button
                            variant="outlined"
                            onClick={() => handleFind(findText)}
                            disabled={!findText.trim()}
                            color="inherit"
                            sx={{ textTransform: "capitalize", fontSize: "0.8rem" }}
                        >
                            Find
                        </Button>
                        <Button
                            variant="outlined"
                            onClick={() => handleFindPrev()}
                            size="small"
                            disabled={disableFindPrev}
                            color="inherit"
                            sx={{ textTransform: "capitalize", fontSize: "0.8rem" }}
                        >
                            Find Previous
                        </Button>
                        <Button
                            variant="outlined"
                            onClick={() => handleFindNext()}
                            size="small"
                            disabled={disableFindNext}
                            color="inherit"
                            sx={{ textTransform: "capitalize", fontSize: "0.8rem" }}
                        >
                            Find Next
                        </Button>
                    </Box>

                    {/* Replace Section */}
                    {isColumnEditable && (
                        <>
                            <Box sx={{ display: "flex", alignItems: "center", mt: 3 }}>
                                <TextField
                                    label="Replace Text"
                                    variant="outlined"
                                    size="small"
                                    value={replaceText}
                                    onChange={(e) => setReplaceText(e.target.value)}
                                    onKeyDown={(e) => e.key === "Enter" && handleReplace(findText, replaceText)}
                                    sx={{ flexGrow: 1 }}
                                />
                            </Box>

                            <Box sx={{ display: "flex", gap: 3 }}>
                                <Button
                                    variant="outlined"
                                    onClick={() => handleReplace(findText, replaceText)}
                                    disabled={!replaceText.trim()}
                                    color="inherit"
                                    sx={{ textTransform: "capitalize", fontSize: "0.8rem" }}
                                >
                                    Replace
                                </Button>
                                <Button
                                    variant="outlined"
                                    onClick={() => handleReplaceAll(findText, replaceText)}
                                    disabled={!replaceText.trim()}
                                    color="inherit"
                                    sx={{ textTransform: "capitalize", fontSize: "0.8rem" }}
                                >
                                    Replace All
                                </Button>
                            </Box>
                        </>
                    )}
                </Box>
            </DialogContent>
            <DialogActions>
                <Button color="error" onClick={onDialogClose}>
                    Close
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default FindAndReplaceDialog;
