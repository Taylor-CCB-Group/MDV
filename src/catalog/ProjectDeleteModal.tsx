import { Close as CloseIcon, Edit as EditIcon } from "@mui/icons-material";
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

export interface ProjectDeleteModalProps {
    id: string;
    open: boolean;
    onDelete: (id: string) => Promise<void>;
    onClose: () => void;
}

const ProjectDeleteModal: React.FC<ProjectDeleteModalProps> = ({
    id,
    open,
    onDelete,
    onClose,
}) => {
    const handleDelete = async () => {
        await onDelete(id);
        onClose();
    };

    return (
        <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
            <DialogTitle>
                Move Project to Recycle Bin
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Typography variant="body1" align="center">
                    Are you sure you want to move this project to the recycle bin?
                </Typography>
            </DialogContent>
            <DialogActions>
                <Box
                    sx={{
                        display: "flex",
                        justifyContent: "flex-end",
                        width: "100%",
                    }}
                >
                    <Button onClick={onClose} color="primary" sx={{ mr: 2 }}>
                        Cancel
                    </Button>
                    <Button onClick={handleDelete} color="error">
                        Move to Recycle Bin
                    </Button>
                </Box>
            </DialogActions>
        </Dialog>
    );
};

export default ProjectDeleteModal;
