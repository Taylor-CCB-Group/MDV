import {
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Typography,
} from "@mui/material";

const SaveViewDialogComponent = (props: {
    open: boolean;
    onClose: () => void;
    action?: () => void;
}) => {
    const { viewManager } = window.mdv.chartManager;

    const onSave = () => {
        viewManager.saveView();
    };

    return (
        <Dialog
            open={props.open}
            onClose={props.onClose}
            fullWidth
            maxWidth="xs"
        >
            <DialogTitle>Save View</DialogTitle>
            <DialogContent dividers>
                <Typography>Do you want to save the current view?</Typography>
            </DialogContent>
            <DialogActions>
                <Button
                    color="primary"
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
            </DialogActions>
        </Dialog>
    );
};

export default SaveViewDialogComponent;
