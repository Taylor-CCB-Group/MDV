import { Button, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from "@mui/material";
import type { FormEvent, KeyboardEvent } from "react";

const DeleteViewDialogComponent = (props: {
    open: boolean;
    onClose: () => void;
}) => {
    const { viewManager } = window.mdv.chartManager;

    const onDelete = () => {
        viewManager.deleteView();
        props.onClose();
    };

    const handleSubmit = (event: FormEvent<HTMLFormElement>) => {
        event.preventDefault();
        onDelete();
    };

    const handleKeyDown = (event: KeyboardEvent<HTMLDivElement>) => {
        if (event.key !== "Enter") return;
        event.preventDefault();
        onDelete();
    };

    return (
        <Dialog
            open={props.open}
            onClose={props.onClose}
            onKeyDown={handleKeyDown}
            fullWidth
            maxWidth="xs"
            PaperProps={{
                component: "form",
                onSubmit: handleSubmit,
            }}
        >
            <DialogTitle>Delete View</DialogTitle>
            <DialogContent dividers>
                <Typography>Do you want to delete the current view?</Typography>
            </DialogContent>
            <DialogActions>
                <Button color="primary" type="submit" autoFocus>
                    Yes
                </Button>
                <Button color="error" onClick={props.onClose}>
                    No
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default DeleteViewDialogComponent;
