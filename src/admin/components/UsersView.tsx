import type { FormEvent } from "react";
import {
    Alert,
    Box,
    Button,
    Chip,
    MenuItem,
    Paper,
    Stack,
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableRow,
    TextField,
    Typography,
} from "@mui/material";
import type { AdminPermission, AdminProject, AdminUser } from "../api";
import {
    formatUserName,
    type InitialProjectAccess,
    type PermissionSetter,
} from "../adminUtils";

type UsersViewProps = {
    users: AdminUser[];
    projects: AdminProject[];
    email: string;
    firstName: string;
    lastName: string;
    projectId: string;
    permission: AdminPermission;
    initialProjectAccess: InitialProjectAccess[];
    saving: boolean;
    onEmailChange: (value: string) => void;
    onFirstNameChange: (value: string) => void;
    onLastNameChange: (value: string) => void;
    onProjectIdChange: (value: string) => void;
    onPermissionChange: (value: string, setter: PermissionSetter) => void;
    onPermissionStateChange: PermissionSetter;
    onAddInitialProjectAccess: () => void;
    onRemoveInitialProjectAccess: (projectId: number) => void;
    onCreateUser: (event: FormEvent<HTMLFormElement>) => void;
};

export function UsersView({
    users,
    projects,
    email,
    firstName,
    lastName,
    projectId,
    permission,
    initialProjectAccess,
    saving,
    onEmailChange,
    onFirstNameChange,
    onLastNameChange,
    onProjectIdChange,
    onPermissionChange,
    onPermissionStateChange,
    onAddInitialProjectAccess,
    onRemoveInitialProjectAccess,
    onCreateUser,
}: UsersViewProps) {
    return (
        <Stack spacing={2.5}>
            <Paper
                component="form"
                onSubmit={onCreateUser}
                className="admin-showcase-panel admin-compact-panel"
                variant="outlined"
            >
                <Stack spacing={2}>
                    <Box>
                        <Typography variant="h5">Create Deployment User</Typography>
                        <Typography variant="body2">
                            Create or resolve a deployment user and grant project access.
                        </Typography>
                    </Box>
                    <TextField
                        label="Email"
                        value={email}
                        onChange={(event) => onEmailChange(event.target.value)}
                        required
                        size="small"
                        type="email"
                    />
                    <Stack direction={{ xs: "column", sm: "row" }} spacing={1}>
                        <TextField
                            label="First name"
                            value={firstName}
                            onChange={(event) => onFirstNameChange(event.target.value)}
                            size="small"
                            fullWidth
                        />
                        <TextField
                            label="Last name"
                            value={lastName}
                            onChange={(event) => onLastNameChange(event.target.value)}
                            size="small"
                            fullWidth
                        />
                    </Stack>
                    <Stack direction={{ xs: "column", sm: "row" }} spacing={1}>
                        <TextField
                            select
                            label="Project"
                            value={projectId}
                            onChange={(event) => onProjectIdChange(event.target.value)}
                            size="small"
                            disabled={projects.length === 0}
                            fullWidth
                        >
                            {projects.map((project) => (
                                <MenuItem key={project.id} value={String(project.id)}>
                                    {project.name}
                                </MenuItem>
                            ))}
                        </TextField>
                        <TextField
                            select
                            label="Permission"
                            value={permission}
                            onChange={(event) => onPermissionChange(event.target.value, onPermissionStateChange)}
                            required
                            size="small"
                            sx={{ minWidth: { sm: 150 } }}
                        >
                            <MenuItem value="view">View</MenuItem>
                            <MenuItem value="edit">Edit</MenuItem>
                            <MenuItem value="owner">Owner</MenuItem>
                        </TextField>
                    </Stack>
                    <Button
                        type="button"
                        variant="outlined"
                        disabled={projects.length === 0}
                        onClick={onAddInitialProjectAccess}
                    >
                        Add initial access
                    </Button>
                    {initialProjectAccess.length === 0 ? (
                        <Alert severity="info">
                            User will be created without project access unless an assignment is added.
                        </Alert>
                    ) : (
                        <Stack spacing={1}>
                            {initialProjectAccess.map((access) => {
                                const project = projects.find((item) => item.id === access.projectId);
                                return (
                                    <Box key={access.projectId} className="admin-showcase-access-pill">
                                        <Typography variant="body2">
                                            {project?.name ?? access.projectId} - {access.permission}
                                        </Typography>
                                        <Button
                                            type="button"
                                            size="small"
                                            color="error"
                                            onClick={() => onRemoveInitialProjectAccess(access.projectId)}
                                        >
                                            Remove
                                        </Button>
                                    </Box>
                                );
                            })}
                        </Stack>
                    )}
                    <Button type="submit" variant="contained" disabled={saving}>
                        {saving ? "Creating..." : "Create user"}
                    </Button>
                </Stack>
            </Paper>

            <Paper className="admin-showcase-panel" variant="outlined">
                <Stack spacing={2}>
                    <Stack direction="row" justifyContent="space-between" alignItems="center">
                        <Box>
                            <Typography variant="h5">Deployment Users</Typography>
                            <Typography variant="body2">Users known to this MDV deployment</Typography>
                        </Box>
                        <Chip label={users.length} size="small" />
                    </Stack>
                    {users.length === 0 ? (
                        <Alert severity="info">No users found in the MDV database yet.</Alert>
                    ) : (
                        <Box className="admin-showcase-table-wrap">
                            <Table size="small">
                                <TableHead>
                                    <TableRow>
                                        <TableCell>User</TableCell>
                                        <TableCell>Email</TableCell>
                                        <TableCell>Status</TableCell>
                                        <TableCell>Admin</TableCell>
                                    </TableRow>
                                </TableHead>
                                <TableBody>
                                    {users.map((user) => (
                                        <TableRow key={user.id}>
                                            <TableCell>{formatUserName(user)}</TableCell>
                                            <TableCell>{user.email}</TableCell>
                                            <TableCell>{user.isActive ? "Active" : "Inactive"}</TableCell>
                                            <TableCell>{user.isAdmin ? "Yes" : "No"}</TableCell>
                                        </TableRow>
                                    ))}
                                </TableBody>
                            </Table>
                        </Box>
                    )}
                </Stack>
            </Paper>
        </Stack>
    );
}
