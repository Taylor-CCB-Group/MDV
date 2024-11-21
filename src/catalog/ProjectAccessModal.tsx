import { Close as CloseIcon, Edit as EditIcon } from "@mui/icons-material";
import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    FormControl,
    IconButton,
    InputLabel,
    MenuItem,
    Select,
} from "@mui/material";
import type React from "react";
import { useState } from "react";
import type { ProjectAccessType } from "./utils/projectUtils";

interface ProjectAccessModalProps {
    id: string;
    type: ProjectAccessType;
    open: boolean;
    onChangeType: (
        id: string,
        newType: ProjectAccessType,
    ) => Promise<void>;
    onClose: () => void;
}

const ProjectAccessModal: React.FC<ProjectAccessModalProps> = ({
    id,
    type,
    open,
    onChangeType,
    onClose,
}) => {
    const [newType, setNewType] = useState(type);

    const handleSave = () => {
        if (newType !== type) {
            onChangeType(id, newType);
        }
        onClose();
    };

    return (
        <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
            <DialogTitle>
                Edit Project Access Level
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
                <FormControl fullWidth variant="outlined" margin="normal">
                    <InputLabel>Project Type</InputLabel>
                    <Select
                        value={newType}
                        onChange={(e) =>
                            setNewType(
                                e.target.value as ProjectAccessType,
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

export default ProjectAccessModal;