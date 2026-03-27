import { Button, Dialog, DialogActions, DialogContent, DialogTitle, TextField } from "@mui/material";
import { type FormEvent, useCallback, useState } from "react";

const SaveAsViewDialogComponent = (props: {
    open: boolean;
    onClose: () => void;
}) => {
    const { viewManager } = window.mdv.chartManager;
    const [viewName, setViewName] = useState("");

    const checkInvalidName = (name: string) => {
        return viewManager.all_views.includes(name);
    };

    const onSaveAs = useCallback(async () => {
        await viewManager.saveAsView(viewName);
        props.onClose();
    }, [viewManager, viewName, props]);

    const isSaveDisabled = !viewName || checkInvalidName(viewName);

    // Allow enter-to-submit behaviour
    const handleSubmit = (event: FormEvent<HTMLFormElement>) => {
        event.preventDefault();
        if (isSaveDisabled) return;
        onSaveAs();
    };

    return (
        <Dialog
            open={props.open}
            onClose={props.onClose}
            fullWidth
            maxWidth="xs"
            PaperProps={{
                component: "form",
                onSubmit: handleSubmit,
            }}
        >
            <DialogTitle>Save View As...</DialogTitle>
            <DialogContent dividers>
                {/* todo consider adding thumbnail here */}
                <TextField
                    label="Enter New View Name"
                    variant="outlined"
                    fullWidth
                    margin="normal"
                    value={viewName}
                    onChange={(e) => setViewName(e.target.value)}
                    error={checkInvalidName(viewName)}
                    helperText={checkInvalidName(viewName) && "Name already exists"}
                />
            </DialogContent>
            <DialogActions>
                <Button type="submit" disabled={isSaveDisabled}>
                    Save As
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

export default SaveAsViewDialogComponent;
