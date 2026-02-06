import { Button, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from "@mui/material";

export type ConfirmDialogType = {
    open: boolean;
    title: string;
    message: string;
    onClose: () => void;
    onConfirm: () => void;
};

const ConfirmDialog = ({ open, onClose, message, title, onConfirm }: ConfirmDialogType) => {
    const handleConfirm = () => {
        onConfirm();
        onClose();
    };
    
    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="xs">
            <DialogTitle>{title}</DialogTitle>
            <DialogContent dividers>
                <Typography>{message}</Typography>
            </DialogContent>
            <DialogActions>
                <Button color="primary" onClick={handleConfirm}>
                    Yes
                </Button>
                <Button color="error" onClick={onClose}>
                    No
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default ConfirmDialog;
