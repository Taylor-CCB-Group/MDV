import { Button, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from "@mui/material";

const SaveViewDialogComponent = (props: {
    open: boolean;
    onClose: () => void;
    action?: () => void;
    content?: string;
}) => {
    const { viewManager } = window.mdv.chartManager;

    const onSave = async () => {
        await viewManager.saveView();
        props?.action?.();
        props.onClose();
    };

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
                <Button onClick={onSave}>
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
