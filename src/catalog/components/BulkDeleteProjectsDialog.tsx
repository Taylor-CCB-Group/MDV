import {
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Typography,
} from "@mui/material";
import type React from "react";
import { DialogCloseIconButton } from "../ProjectRenameModal";

export interface BulkDeleteProjectsDialogProps {
    open: boolean;
    selectedCount: number;
    onClose: () => void;
    onConfirm: () => void | Promise<void>;
}

const BulkDeleteProjectsDialog: React.FC<BulkDeleteProjectsDialogProps> = ({
    open,
    selectedCount,
    onClose,
    onConfirm,
}) => {
    return (
        <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
            <DialogTitle>
                Move projects to recycle bin?
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Typography>
                    Move {selectedCount} selected project{selectedCount === 1 ? "" : "s"} to the
                    recycle bin?
                </Typography>
            </DialogContent>
            <DialogActions>
                <Button onClick={onClose}>Cancel</Button>
                <Button color="error" onClick={onConfirm}>
                    Move to recycle bin
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default BulkDeleteProjectsDialog;
