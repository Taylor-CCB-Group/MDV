import type { FormEvent } from "react";
import {
    Alert,
    Box,
    Button,
    Chip,
    CircularProgress,
    IconButton,
    MenuItem,
    Paper,
    Stack,
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableRow,
    TextField,
    Tooltip,
    Typography,
} from "@mui/material";
import { Trash2 } from "lucide-react";
import type { AdminPermission, AdminProject, AdminProjectMember, AdminUser } from "../api";
import { formatUserName, permissionFromValue, type PermissionSetter } from "../adminUtils";

type ProjectsViewProps = {
    projects: AdminProject[];
    projectMembers: AdminProjectMember[];
    selectedProject: AdminProject | undefined;
    availableProjectUsers: AdminUser[];
    memberProjectId: string;
    addMemberUserId: string;
    addMemberPermission: AdminPermission;
    memberAction: string | null;
    membersLoading: boolean;
    onMemberProjectIdChange: (value: string) => void;
    onAddMemberUserIdChange: (value: string) => void;
    onPermissionChange: (value: string, setter: PermissionSetter) => void;
    onAddMemberPermissionChange: PermissionSetter;
    onAddProjectMember: (event: FormEvent<HTMLFormElement>) => void;
    onUpdateProjectMember: (userId: number, permission: AdminPermission) => void;
    onRemoveProjectMember: (userId: number) => void;
};

export function ProjectsView({
    projects,
    projectMembers,
    selectedProject,
    availableProjectUsers,
    memberProjectId,
    addMemberUserId,
    addMemberPermission,
    memberAction,
    membersLoading,
    onMemberProjectIdChange,
    onAddMemberUserIdChange,
    onPermissionChange,
    onAddMemberPermissionChange,
    onAddProjectMember,
    onUpdateProjectMember,
    onRemoveProjectMember,
}: ProjectsViewProps) {
    return (
        <Stack spacing={2.5}>
            <Paper className="admin-showcase-panel admin-members-panel" variant="outlined">
                <Stack spacing={2}>
                    <Stack
                        direction={{ xs: "column", md: "row" }}
                        spacing={2}
                        justifyContent="space-between"
                        alignItems={{ xs: "stretch", md: "center" }}
                    >
                        <Box>
                            <Typography variant="h5">Project Access Management</Typography>
                            <Typography variant="body2">
                                {selectedProject
                                    ? `Managing permissions for ${selectedProject.name}`
                                    : "Select a project to manage access"}
                            </Typography>
                        </Box>
                        <TextField
                            select
                            label="Project"
                            value={memberProjectId}
                            onChange={(event) => onMemberProjectIdChange(event.target.value)}
                            size="small"
                            disabled={projects.length === 0}
                            sx={{ minWidth: { md: 280 } }}
                        >
                            {projects.map((project) => (
                                <MenuItem key={project.id} value={String(project.id)}>
                                    {project.name}
                                </MenuItem>
                            ))}
                        </TextField>
                    </Stack>

                    <Paper
                        component="form"
                        variant="outlined"
                        onSubmit={onAddProjectMember}
                        className="admin-showcase-inline-form"
                    >
                        <Stack
                            direction={{ xs: "column", md: "row" }}
                            spacing={1}
                            alignItems={{ xs: "stretch", md: "center" }}
                        >
                            <TextField
                                select
                                label="Existing user"
                                value={addMemberUserId}
                                onChange={(event) => onAddMemberUserIdChange(event.target.value)}
                                size="small"
                                disabled={availableProjectUsers.length === 0}
                                sx={{ minWidth: { md: 260 }, flex: 1 }}
                            >
                                {availableProjectUsers.map((user) => (
                                    <MenuItem key={user.id} value={String(user.id)}>
                                        {user.email}
                                    </MenuItem>
                                ))}
                            </TextField>
                            <TextField
                                select
                                label="Permission"
                                value={addMemberPermission}
                                onChange={(event) => onPermissionChange(event.target.value, onAddMemberPermissionChange)}
                                size="small"
                                sx={{ minWidth: { md: 140 } }}
                            >
                                <MenuItem value="view">View</MenuItem>
                                <MenuItem value="edit">Edit</MenuItem>
                                <MenuItem value="owner">Owner</MenuItem>
                            </TextField>
                            <Button
                                type="submit"
                                variant="contained"
                                disabled={memberAction === "add" || !addMemberUserId || availableProjectUsers.length === 0}
                            >
                                {memberAction === "add" ? "Adding..." : "Add access"}
                            </Button>
                        </Stack>
                    </Paper>

                    {projects.length === 0 ? (
                        <Alert severity="info">No projects found in the MDV database.</Alert>
                    ) : membersLoading ? (
                        <Stack alignItems="center" sx={{ py: 4 }}>
                            <CircularProgress size={24} />
                        </Stack>
                    ) : projectMembers.length === 0 ? (
                        <Alert severity="info">No users are assigned to this project yet.</Alert>
                    ) : (
                        <Box className="admin-showcase-table-wrap">
                            <Table size="small">
                                <TableHead>
                                    <TableRow>
                                        <TableCell>User</TableCell>
                                        <TableCell className="admin-table-col-permission">Permission</TableCell>
                                        <TableCell className="admin-table-col-bool">Read</TableCell>
                                        <TableCell className="admin-table-col-bool">Write</TableCell>
                                        <TableCell className="admin-table-col-bool">Owner</TableCell>
                                        <TableCell className="admin-table-col-actions" align="right">Actions</TableCell>
                                    </TableRow>
                                </TableHead>
                                <TableBody>
                                    {projectMembers.map((member) => (
                                        <TableRow key={`${member.projectAccess.projectId}-${member.user.id}`}>
                                            <TableCell>
                                                <Typography variant="body2">{formatUserName(member.user)}</Typography>
                                                <Typography variant="caption">{member.user.email}</Typography>
                                            </TableCell>
                                            <TableCell className="admin-table-col-permission">
                                                <TextField
                                                    select
                                                    fullWidth
                                                    className="admin-table-permission-select"
                                                    value={member.projectAccess.permission}
                                                    onChange={(event) => {
                                                        const nextPermission = permissionFromValue(event.target.value);
                                                        if (nextPermission) {
                                                            onUpdateProjectMember(member.user.id, nextPermission);
                                                        }
                                                    }}
                                                    size="small"
                                                    disabled={memberAction === `update-${member.user.id}`}
                                                >
                                                    <MenuItem value="view">View</MenuItem>
                                                    <MenuItem value="edit">Edit</MenuItem>
                                                    <MenuItem value="owner">Owner</MenuItem>
                                                </TextField>
                                            </TableCell>
                                            <TableCell className="admin-table-col-bool">{member.projectAccess.canRead ? "Yes" : "No"}</TableCell>
                                            <TableCell className="admin-table-col-bool">{member.projectAccess.canWrite ? "Yes" : "No"}</TableCell>
                                            <TableCell className="admin-table-col-bool">{member.projectAccess.isOwner ? "Yes" : "No"}</TableCell>
                                            <TableCell className="admin-table-col-actions" align="center">
                                                <Tooltip title="Remove access">
                                                    <span>
                                                        <IconButton
                                                            color="error"
                                                            size="small"
                                                            className="admin-table-action-button"
                                                            disabled={memberAction === `remove-${member.user.id}`}
                                                            aria-label="Remove access"
                                                            onClick={() => onRemoveProjectMember(member.user.id)}
                                                        >
                                                            <Trash2 size={16} strokeWidth={2.1} />
                                                        </IconButton>
                                                    </span>
                                                </Tooltip>
                                            </TableCell>
                                        </TableRow>
                                    ))}
                                </TableBody>
                            </Table>
                        </Box>
                    )}
                </Stack>
            </Paper>

            <Paper className="admin-showcase-panel" variant="outlined">
                <Stack spacing={2}>
                    <Stack direction="row" justifyContent="space-between" alignItems="center">
                        <Box>
                            <Typography variant="h5">Projects</Typography>
                            <Typography variant="body2">Existing projects available for assignment</Typography>
                        </Box>
                        <Chip label={projects.length} size="small" />
                    </Stack>
                    {projects.length === 0 ? (
                        <Alert severity="info">No projects found in the MDV database.</Alert>
                    ) : (
                        <Box className="admin-showcase-table-wrap">
                            <Table size="small">
                                <TableHead>
                                    <TableRow>
                                        <TableCell>Name</TableCell>
                                        <TableCell>ID</TableCell>
                                        <TableCell>Access</TableCell>
                                        <TableCell>Status</TableCell>
                                    </TableRow>
                                </TableHead>
                                <TableBody>
                                    {projects.map((project) => (
                                        <TableRow key={project.id}>
                                            <TableCell>{project.name}</TableCell>
                                            <TableCell>{project.id}</TableCell>
                                            <TableCell>{project.accessLevel || "unknown"}</TableCell>
                                            <TableCell>{project.isDeleted ? "Deleted" : "Active"}</TableCell>
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
