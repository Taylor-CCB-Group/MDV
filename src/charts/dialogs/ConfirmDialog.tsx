import { Button, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from "@mui/material";

export type ConfirmDialogType = {
    open: boolean;
    title: string;
    message: string;
    onClose: () => void;
    onConfirm: () => void;
    isDelete?: boolean;
};

const ConfirmDialog = ({ open, onClose, message, title, onConfirm, isDelete }: ConfirmDialogType) => {
    

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="xs">
            <DialogTitle>{title}</DialogTitle>
            <DialogContent dividers>
                <Typography>{message}</Typography>
            </DialogContent>
            <DialogActions>
                <Button color={isDelete ? "error" : "primary"} onClick={onConfirm}>
                    Yes
                </Button>
                <Button color={isDelete ? "primary" : "error"} onClick={onClose}>
                    No
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default ConfirmDialog;
