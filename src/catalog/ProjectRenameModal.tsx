import { Close as CloseIcon } from "@mui/icons-material";
import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    IconButton,
    TextField,
} from "@mui/material";
import type React from "react";
import { useCallback, useState } from "react";

export interface ProjectRenameModalProps {
    id: string;
    name: string;
    open: boolean;
    onRename: (id: string, newName: string) => Promise<void>;
    onClose: () => void;
}

const ProjectRenameModal: React.FC<ProjectRenameModalProps> = ({
    id,
    name,
    open,
    onRename,
    onClose,
}) => {
    const [newName, setNewName] = useState(name);
    const handleSave = useCallback(() => {
        if (newName !== name) {
            onRename(id, newName);
        }
        onClose();
    }, [id, name, newName, onRename, onClose]);

    return (
        <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
            <DialogTitle>
                Rename Project
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
                <Box sx={{ mb: 2 }}>
                    <TextField
                        label="Project Name"
                        value={newName}
                        onChange={(e) => setNewName(e.target.value)}
                        fullWidth
                        variant="outlined"
                        margin="normal"
                        autoFocus
                        onKeyDown={(e) => {
                            if (e.key === "Enter") {
                                handleSave();
                            }
                        }}
                    />
                </Box>
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
                    <Button onClick={handleSave} color="primary">
                        Save Changes
                    </Button>
                </Box>
            </DialogActions>
        </Dialog>
    );
};

export default ProjectRenameModal;
