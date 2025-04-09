import { Close as CloseIcon } from "@mui/icons-material";
import {
    Autocomplete,
    Box,
    Button,
    CircularProgress,
    Dialog,
    DialogActions,
    DialogContent,
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
} from "@mui/material";
import type React from "react";
import DeleteIcon from "@mui/icons-material/Delete";
import useProjectShare, { type UserPermission, type RegisteredUser } from "./hooks/useProjectShare";
import ErrorDisplay from "@/charts/dialogs/ErrorDisplay";

export interface ProjectShareModalProps {
    open: boolean;
    onClose: () => void;
    projectId: string;
}

const ProjectShareModal: React.FC<ProjectShareModalProps> = ({ open, onClose, projectId }) => {
    const {
        email,
        setEmail,
        sharedUsers,
        setSharedUsers,
        addUser,
        // updateSharedUsers,
        isLoading,
        error,
        userList,
        newUser,
        setNewUser,
        deleteSharedUser,
        changeUserPermission
    } = useProjectShare(projectId);

    const handleAddNewUser = async () => {
        if (newUser) {
            // todo: update api logic
            // await addUser(newUser.id, "View");
            console.log("new user", newUser);
            setEmail("");
            setNewUser(undefined);
        }
    };

    const handlePermissionChange = async (permission: UserPermission, userId: number) => {
        await changeUserPermission(userId, permission);
    };

    const handleDeleteUser = async (userId: number) => {
       await deleteSharedUser(userId);
    };

    const onInputChange = (e: React.SyntheticEvent, value: string) => {
        console.log("oninputchange", value);
        setEmail(value);
    };

    const onChange = async (e: React.SyntheticEvent, user: RegisteredUser | null) => {
        console.log("onchange", user);
        if (user) {
            setNewUser(user);
        }
    }

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
                    {error ? (
                        <div>
                            <ErrorDisplay error={error} />
                        </div>
                    ) : (
                        <>
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
                                    onChange={onChange}
                                    onInputChange={onInputChange}
                                    getOptionLabel={(option) => option.email}
                                    getOptionKey={(option) => option.id}
                                />
                                <Button
                                    variant="contained"
                                    color="primary"
                                    onClick={handleAddNewUser}
                                    disabled={isLoading || !email}
                                >
                                    {isLoading ? <CircularProgress size="1.5rem" /> : "Add"}
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
                                        <TableRow key={user.id}>
                                            <TableCell>{user.email}</TableCell>
                                            <TableCell>
                                                <Select
                                                    value={user.permission}
                                                    size="small"
                                                    fullWidth
                                                    onChange={(e) => handlePermissionChange(e.target.value as UserPermission, user.id)}
                                                >
                                                    <MenuItem value="View">View</MenuItem>
                                                    <MenuItem value="Edit">Edit</MenuItem>
                                                    <MenuItem value="Owner">Owner</MenuItem>
                                                </Select>
                                            </TableCell>
                                            <TableCell align="right">
                                                <IconButton onClick={() => handleDeleteUser(user.id)}>
                                                    <DeleteIcon />
                                                </IconButton>
                                            </TableCell>
                                        </TableRow>
                                    ))}
                                </TableBody>
                            </Table>
                        </>
                    )}
                </Box>
            </DialogContent>
            <DialogActions sx={{ py: 2 }}>
                {/* <Button
                    onClick={handleUpdate}
                    variant="outlined"
                    color="primary"
                    sx={{ fontWeight: "bold" }}
                    disabled={isLoading}
                >
                    {isLoading ? <CircularProgress size="1.5rem" /> : "Update"}
                </Button> */}
                <Button onClick={onClose} color="error" variant="outlined" sx={{ fontWeight: "bold" }}>
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default ProjectShareModal;
