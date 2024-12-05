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

interface ProjectDeleteModalProps {
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
                Delete Project
                <IconButton
                    aria-label="close"
                    onClick={onClose}
                    sx={{
                        position: "absolute",
                        right: 8,
                        top: 8,
                        color: (theme) => theme.palette.grey[500],
                    }}
                >
                    <CloseIcon />
                </IconButton>
            </DialogTitle>
            <DialogContent dividers>
                <Typography variant="body1" align="center">
                    Are you sure you want to permanently{" "}
                    <strong>delete this project</strong>? This action is{" "}
                    <strong>irreversible</strong> and cannot be{" "}
                    <strong>undone</strong>.
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
                        Delete Project
                    </Button>
                </Box>
            </DialogActions>
        </Dialog>
    );
};

export default ProjectDeleteModal;
