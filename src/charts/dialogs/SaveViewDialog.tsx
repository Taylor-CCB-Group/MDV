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
            if (!props.open) return; //probably doesn't happen, but no harm to check
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

        // seems ok to use window here since this is a modal dialog
        // saves faffing about with refs
        window.addEventListener("keydown", handleKeyDown);
        return () => {
            window.removeEventListener("keydown", handleKeyDown);
        };
    }, [props, onSave, props.open]);

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
