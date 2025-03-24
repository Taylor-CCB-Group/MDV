import {
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Typography,
} from "@mui/material";

const DeleteViewDialogComponent = (props: {
    open: boolean;
    onClose: () => void;
}) => {
    const { viewManager } = window.mdv.chartManager;

    const onDelete = () => {
        viewManager.deleteView();

        props.onClose();
    };
    return (
        <Dialog
            open={props.open}
            onClose={props.onClose}
            fullWidth
            maxWidth="xs"
        >
            <DialogTitle>Delete View</DialogTitle>
            <DialogContent dividers>
                <Typography>Do you want to delete the current view?</Typography>
            </DialogContent>
            <DialogActions>
                <Button color="primary" onClick={onDelete}>
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
