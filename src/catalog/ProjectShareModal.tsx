import { Close as CloseIcon, ContentCopy as ContentCopyIcon, Email as EmailIcon } from "@mui/icons-material";
import {
    Autocomplete,
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    IconButton,
    MenuItem,
    Select,
    Snackbar,
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableRow,
    TextField,
    Typography,
    useTheme,
} from "@mui/material";
import type React from "react";
import { useState } from "react";
import AddIcon from "@mui/icons-material/Add";
import DeleteIcon from "@mui/icons-material/Delete";
import useProjectShare from "./hooks/useProjectShare";

export interface ProjectShareModalProps {
    open: boolean;
    onClose: () => void;
    projectId: string;
}

const ProjectShareModal: React.FC<ProjectShareModalProps> = ({ open, onClose, projectId }) => {
    const { username, setUsername, sharedUsers, setSharedUsers, addUser, updateSharedUsers } =
        useProjectShare(projectId);

    const theme = useTheme();

    const handleAddNewUser = async () => {
        if (username) {
            // todo: update api logic
            // await addUser(username);
            setUsername("");
        }
    };

    const handlePermissionChange = async (permission: string, index: number) => {
        setSharedUsers((prev) => {
            const updatedUsers = [...prev];
            updatedUsers[index].permission = permission;
            return updatedUsers;
        });
    };

    const handleDeleteUser = async (index: number) => {
        setSharedUsers((prev) => {
            const updatedUsers = [...prev];
            updatedUsers.splice(index, 1);
            return updatedUsers;
        });
    };

    const handleUpdate = async () => {
        // todo: update api logic
        // await updateSharedUsers(sharedUsers);
    };

    return (
        <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
            <DialogTitle>
                Manage Project Sharing
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
                <Box sx={{ mt: 2 }}>
                    {/* <Typography variant="subtitle1" gutterBottom>
                        Add New User
                    </Typography> */}
                    <Box sx={{ display: "flex", gap: 1, mb: 5 }}>
                        <Autocomplete
                            fullWidth
                            size="small"
                            options={[]}
                            sx={{ color: "inherit" }}
                            renderInput={(params) => (
                                <TextField {...params} placeholder="Enter username of the user you want to add" />
                            )}
                            value={username}
                            onChange={(_, value) => value && setUsername(value)}
                        />
                        <Button variant="contained" color="primary" onClick={handleAddNewUser} endIcon={<AddIcon />}>
                            Add
                        </Button>
                    </Box>
                    <Table
                        sx={{
                            mb: 3,
                            "& .MuiTableCell-head": {
                                fontWeight: 600,
                            },
                        }}
                    >
                        <TableHead>
                            <TableRow>
                                <TableCell>Shared Users</TableCell>
                                <TableCell>Permission</TableCell>
                                <TableCell align="right">Remove</TableCell>
                            </TableRow>
                        </TableHead>
                        <TableBody>
                            {sharedUsers.map((user, index) => (
                                <TableRow key={`${user.name}-${index}`}>
                                    <TableCell>{user.name}</TableCell>
                                    <TableCell>
                                        <Select
                                            value={user.permission}
                                            size="small"
                                            fullWidth
                                            onChange={(e) => handlePermissionChange(e.target.value, index)}
                                        >
                                            <MenuItem value="view">View</MenuItem>
                                            <MenuItem value="edit">Edit</MenuItem>
                                        </Select>
                                    </TableCell>
                                    <TableCell align="right">
                                        <IconButton onClick={() => handleDeleteUser(index)}>
                                            <DeleteIcon />
                                        </IconButton>
                                    </TableCell>
                                </TableRow>
                            ))}
                        </TableBody>
                    </Table>
                </Box>
            </DialogContent>
            <DialogActions sx={{ py: 2 }}>
                <Button onClick={handleUpdate} variant="outlined" color="primary" sx={{ fontWeight: "bold" }}>
                    Update
                </Button>
                <Button onClick={onClose} color="error" variant="outlined" sx={{ fontWeight: "bold" }}>
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default ProjectShareModal;
