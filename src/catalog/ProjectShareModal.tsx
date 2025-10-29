import { Close as CloseIcon } from "@mui/icons-material";
import {
    Autocomplete,
    Box,
    Button,
    CircularProgress,
    Dialog,
    DialogActions,
    DialogContent,
    DialogContentText,
    DialogTitle,
    IconButton,
    MenuItem,
    Select,
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableRow,
    TextField,
    Typography,
} from "@mui/material";
import type React from "react";
import DeleteIcon from "@mui/icons-material/Delete";
import useProjectShare, { type UserPermission, type RegisteredUser } from "./hooks/useProjectShare";
import { useState } from "react";
import { DialogCloseIconButton } from "./ProjectRenameModal";
import { matchEmail } from "@/lib/utils";
import usePermissions from "./PermissionsContext";

export interface ProjectShareModalProps {
    open: boolean;
    onClose: () => void;
    projectId: string;
}

const ProjectShareModal: React.FC<ProjectShareModalProps> = ({ open, onClose, projectId }) => {
    const [openTextDialog, setOpenTextDialog] = useState(false);
    const {
        email,
        setEmail,
        sharedUsers,
        addUser,
        isLoading,
        error,
        errorMsg,
        userList,
        newUser,
        setNewUser,
        deleteSharedUser,
        changeUserPermission,
    } = useProjectShare(projectId);
    const { permissions: operationPermissions } = usePermissions();

    const handleAddNewUser = async () => {
        if (newUser) {
            await addUser(newUser.id, "View");
            setEmail("");
            setNewUser(null);
        }
    };

    const handlePermissionChange = async (permission: UserPermission, userId: number) => {
        await changeUserPermission(userId, permission);
    };

    const handleDeleteUser = async (userId: number) => {
        await deleteSharedUser(userId);
    };

    const onInputChange = (value: string) => {
        setEmail(value);
    };

    const onChange = (user: RegisteredUser | null) => {
        if (user) {
            setNewUser(user);
        } else {
            setNewUser(null);
            setEmail("");
        }
    };

    return (
        <>
            <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
                <DialogTitle>
                    Manage Project Sharing
                    <DialogCloseIconButton onClose={onClose} />
                </DialogTitle>
                <DialogContent dividers>
                    <Box sx={{ mt: 2 }}>
                        {error ? (
                            <Box
                                sx={{
                                    display: "flex",
                                    justifyContent: "center",
                                    alignItems: "center",
                                    height: "20vh",
                                }}
                            >
                                <Typography variant="h6" color="error">
                                    {error}
                                </Typography>
                            </Box>
                        ) : (
                            <>
                                {errorMsg && (
                                    <Box sx={{ mb: 2 }}>
                                        <Typography color="error">{errorMsg}</Typography>
                                    </Box>
                                )}
                                <Box sx={{ display: "flex", gap: 1, mb: 5 }}>
                                    <Autocomplete
                                        fullWidth
                                        size="small"
                                        options={userList}
                                        sx={{ color: "inherit" }}
                                        renderInput={(params) => (
                                            <TextField {...params} placeholder="Enter email to search for the user" />
                                        )}
                                        inputValue={email}
                                        onChange={(_, value) => onChange(value)}
                                        onInputChange={(_, value) => onInputChange(value)}
                                        getOptionLabel={(option) => option.email}
                                        // Add isOptionEqualToValue here if it doesn't work properly
                                        noOptionsText={
                                            <Box
                                                sx={{
                                                    p: 1,
                                                    display: "flex",
                                                    justifyContent: "center",
                                                    alignItems: "center",
                                                }}
                                            >
                                                <Typography variant="body2" color="textPrimary">
                                                    If it's a new user:
                                                </Typography>
                                                <Button
                                                    variant="contained"
                                                    size="small"
                                                    sx={{ ml: 1 }}
                                                    onClick={() => setOpenTextDialog(true)}
                                                >
                                                    Add New User
                                                </Button>
                                            </Box>
                                        }
                                        filterOptions={(options, state) => {
                                            if (state.inputValue.length < 5) return [];
                                            const input = state.inputValue.toLowerCase();
                                            return options.filter((option) => {
                                                const email = option.email.toLowerCase();
                                                if (matchEmail(input, email)) return option;
                                            });
                                        }}
                                    />
                                    <Button
                                        variant="contained"
                                        color="primary"
                                        onClick={handleAddNewUser}
                                        disabled={isLoading || !newUser}
                                    >
                                        {isLoading ? <CircularProgress size="1.5rem" /> : "Add"}
                                    </Button>
                                </Box>
                                {sharedUsers?.length === 0 ? (
                                    <Box
                                        sx={{
                                            display: "flex",
                                            justifyContent: "center",
                                            my: 2,
                                        }}
                                    >
                                        <Typography variant="h6">No shared users</Typography>
                                    </Box>
                                ) : (
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
                                                {operationPermissions.editUserPermissions && <TableCell>Permission</TableCell>}
                                                {operationPermissions.removeUserFromProject && <TableCell align="right">Remove</TableCell>}
                                            </TableRow>
                                        </TableHead>
                                        <TableBody>
                                            {sharedUsers.map((user) => (
                                                <TableRow key={user.id}>
                                                    <TableCell>{user.email}</TableCell>
                                                    {operationPermissions.editUserPermissions && <TableCell>
                                                        <Select
                                                            value={user.permission}
                                                            size="small"
                                                            fullWidth
                                                            onChange={(e) =>
                                                                handlePermissionChange(
                                                                    e.target.value as UserPermission,
                                                                    user.id,
                                                                )
                                                            }
                                                        >
                                                            <MenuItem value="View">View</MenuItem>
                                                            <MenuItem value="Edit">Edit</MenuItem>
                                                            <MenuItem value="Owner">Owner</MenuItem>
                                                        </Select>
                                                    </TableCell>}
                                                    {operationPermissions.removeUserFromProject && <TableCell align="right">
                                                        <IconButton onClick={() => handleDeleteUser(user.id)}>
                                                            <DeleteIcon />
                                                        </IconButton>
                                                    </TableCell>}
                                                </TableRow>
                                            ))}
                                        </TableBody>
                                    </Table>
                                )}
                            </>
                        )}
                    </Box>
                </DialogContent>
                <DialogActions sx={{ py: 2 }}>
                    <Button onClick={onClose} color="error" variant="outlined" sx={{ fontWeight: "bold" }}>
                        Cancel
                    </Button>
                </DialogActions>
            </Dialog>
            {openTextDialog && (
                <Dialog open={openTextDialog} onClose={() => setOpenTextDialog(false)} fullWidth maxWidth="sm">
                    <DialogContent
                        dividers
                        sx={{ display: "flex", height: "20vh", justifyContent: "center", alignItems: "center" }}
                    >
                        <Typography>Please contact the administrator:</Typography>
                        <Typography sx={{ ml: 1 }} color="info">
                            stephen.taylor@well.ox.ac.uk
                        </Typography>
                    </DialogContent>
                    <DialogActions>
                        <Button onClick={() => setOpenTextDialog(false)} color="primary">
                            OK
                        </Button>
                    </DialogActions>
                </Dialog>
            )}
        </>
    );
};

export default ProjectShareModal;
