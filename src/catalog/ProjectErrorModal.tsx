import {
    Close as CloseIcon,
    ErrorOutline as ErrorIcon,
} from "@mui/icons-material";
import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    IconButton,
    Typography,
} from "@mui/material";
import type React from "react";
import { DialogCloseIconButton } from "./ProjectRenameModal";

export interface ErrorModalProps {
    open: boolean;
    message: string;
    onClose: () => void;
}

const ErrorModal: React.FC<ErrorModalProps> = ({ open, message, onClose }) => {
    return (
        <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
            <DialogTitle>
                <Box display="flex" alignItems="center">
                    <ErrorIcon color="error" sx={{ mr: 1 }} />
                    Error
                </Box>
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Typography variant="body1">{message}</Typography>
            </DialogContent>
            <DialogActions>
                <Box
                    sx={{
                        display: "flex",
                        justifyContent: "flex-end",
                        width: "100%",
                    }}
                >
                    <Button
                        onClick={onClose}
                        variant="contained"
                        color="primary"
                    >
                        Close
                    </Button>
                </Box>
            </DialogActions>
        </Dialog>
    );
};

export default ErrorModal;
