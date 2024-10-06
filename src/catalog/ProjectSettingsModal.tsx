import type React from "react";
import { useState } from "react";
import {
    Dialog,
    DialogTitle,
    DialogContent,
    DialogActions,
    TextField,
    Select,
    MenuItem,
    Button,
    IconButton,
    FormControl,
    InputLabel,
    Box,
} from "@mui/material";
import { Close as CloseIcon, Edit as EditIcon } from "@mui/icons-material";

interface ProjectSettingsModalProps {
    id: string;
    name: string;
    type: "Editable" | "Read-Only";
    open: boolean;
    onRename: (id: string, newName: string) => void;
    onChangeType: (id: string, newType: "Editable" | "Read-Only") => void;
    onDelete: (id: string) => void;
    onClose: () => void;
}

const ProjectSettingsModal: React.FC<ProjectSettingsModalProps> = ({
    id,
    name,
    type,
    open,
    onRename,
    onChangeType,
    onDelete,
    onClose,
}) => {
    const [newName, setNewName] = useState(name);
    const [newType, setNewType] = useState(type);
    const [isEditingName, setIsEditingName] = useState(false);

    const handleSave = () => {
        if (newName !== name) {
            onRename(id, newName);
        }
        if (newType !== type) {
            onChangeType(id, newType);
        }
        onClose();
    };

    const toggleEditName = () => {
        setIsEditingName((prev) => !prev);
    };

    const handleDelete = () => {
        if (window.confirm("Are you sure you want to delete this project?")) {
            onDelete(id);
            onClose();
        }
    };

    return (
        <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
            <DialogTitle>
                Project Settings
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
                        InputProps={{
                            readOnly: !isEditingName,
                            endAdornment: (
                                <IconButton onClick={toggleEditName} edge="end">
                                    <EditIcon />
                                </IconButton>
                            ),
                        }}
                    />
                </Box>
                <FormControl fullWidth variant="outlined" margin="normal">
                    <InputLabel>Project Type</InputLabel>
                    <Select
                        value={newType}
                        onChange={(e) =>
                            setNewType(
                                e.target.value as "Editable" | "Read-Only",
                            )
                        }
                        label="Project Type"
                    >
                        <MenuItem value="Editable">Editable</MenuItem>
                        <MenuItem value="Read-Only">Read-Only</MenuItem>
                    </Select>
                </FormControl>
            </DialogContent>
            <DialogActions>
                <Box
                    sx={{
                        display: "flex",
                        justifyContent: "space-between",
                        width: "100%",
                    }}
                >
                    <Button onClick={handleDelete} color="error">
                        Delete Project
                    </Button>
                    <Button onClick={handleSave} color="primary">
                        Save Changes
                    </Button>
                </Box>
            </DialogActions>
        </Dialog>
    );
};

export default ProjectSettingsModal;
