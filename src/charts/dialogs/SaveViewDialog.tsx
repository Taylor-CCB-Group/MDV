import { Button, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from "@mui/material";
import { useEffect, useCallback } from "react";

const SaveViewDialogComponent = (props: {
    open: boolean;
    onClose: () => void;
    action?: () => void;
    content?: string;
}) => {
    const { viewManager } = window.mdv.chartManager;

    const onSave = useCallback(() => {
        viewManager.saveView();
    }, [viewManager]);

    useEffect(() => {
        const handleKeyDown = (event: KeyboardEvent) => {
            if (event.key === "y") {
                event.preventDefault();
                onSave();
                props?.action?.();
                props.onClose();
            }
            if (event.key === "n") {
                event.preventDefault();
                props?.action?.();
                props.onClose();
            }            
        };

        window.addEventListener("keydown", handleKeyDown);
        return () => {
            window.removeEventListener("keydown", handleKeyDown);
        };
    }, [props, onSave]);

    return (
        <Dialog open={props.open} onClose={props.onClose} fullWidth maxWidth="xs">
            <DialogTitle>Save View</DialogTitle>
            <DialogContent dividers>
                <Typography>
                    {props?.content
                        ? props?.content
                        : "You have unsaved changes, do you want to save the current view?"}
                </Typography>
            </DialogContent>
            <DialogActions>
                <Button
                    onClick={() => {
                        onSave();
                        props?.action?.();
                        props.onClose();
                    }}
                >
                    Yes
                </Button>
                <Button
                    color="error"
                    onClick={() => {
                        props?.action?.();
                        props.onClose();
                    }}
                >
                    No
                </Button>
                <Button
                    color="error"
                    onClick={() => {
                        props.onClose();
                    }}
                >
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default SaveViewDialogComponent;
