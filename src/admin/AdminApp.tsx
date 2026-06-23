import { type FormEvent, useCallback, useEffect, useState } from "react";
import {
    Alert,
    Box,
    Button,
    Chip,
    CircularProgress,
    Container,
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
import {
    adminApi,
    type AdminPermission,
    type AdminProject,
    type AdminProjectMember,
    type AdminSession,
    type AdminUser,
} from "./api";

function formatUserName(user: AdminUser) {
    const name = [user.firstName, user.lastName].filter(Boolean).join(" ");
    return name || user.email;
}

type InitialProjectAccess = {
    projectId: number;
    permission: AdminPermission;
};

export default function AdminApp() {
    const [session, setSession] = useState<AdminSession | null>(null);
    const [users, setUsers] = useState<AdminUser[]>([]);
    const [projects, setProjects] = useState<AdminProject[]>([]);
    const [projectMembers, setProjectMembers] = useState<AdminProjectMember[]>([]);
    const [loading, setLoading] = useState(true);
    const [membersLoading, setMembersLoading] = useState(false);
    const [saving, setSaving] = useState(false);
    const [memberAction, setMemberAction] = useState<string | null>(null);
    const [error, setError] = useState<string | null>(null);
    const [success, setSuccess] = useState<string | null>(null);
    const [email, setEmail] = useState("");
    const [firstName, setFirstName] = useState("");
    const [lastName, setLastName] = useState("");
    const [projectId, setProjectId] = useState("");
    const [initialProjectAccess, setInitialProjectAccess] = useState<InitialProjectAccess[]>([]);
    const [memberProjectId, setMemberProjectId] = useState("");
    const [permission, setPermission] = useState<AdminPermission>("view");
    const [addMemberUserId, setAddMemberUserId] = useState("");
    const [addMemberPermission, setAddMemberPermission] = useState<AdminPermission>("view");

    const memberUserIds = new Set(projectMembers.map((member) => member.user.id));
    const availableProjectUsers = users.filter((user) => !memberUserIds.has(user.id));

    const loadAdminData = useCallback(async () => {
        const [sessionResult, usersResult, projectsResult] = await Promise.all([
            adminApi.session(),
            adminApi.users(),
            adminApi.projects(),
        ]);
        setSession(sessionResult);
        setUsers(usersResult.users);
        setProjects(projectsResult.projects);
        if (projectsResult.projects.length > 0) {
            const firstProjectId = String(projectsResult.projects[0].id);
            setProjectId((current) => current || firstProjectId);
            setMemberProjectId((current) => current || firstProjectId);
        }
    }, []);

    const loadProjectMembers = useCallback(async (selectedProjectId: string) => {
        if (!selectedProjectId) {
            setProjectMembers([]);
            return;
        }
        const numericProjectId = Number(selectedProjectId);
        if (!Number.isInteger(numericProjectId)) {
            setProjectMembers([]);
            return;
        }
        setMembersLoading(true);
        try {
            const result = await adminApi.projectMembers(numericProjectId);
            setProjectMembers(result.members);
        } catch (err) {
            setError(err instanceof Error ? err.message : "Failed to load project members");
        } finally {
            setMembersLoading(false);
        }
    }, []);

    useEffect(() => {
        let cancelled = false;

        async function loadInitialState() {
            setLoading(true);
            setError(null);
            try {
                await loadAdminData();
                if (cancelled) return;
            } catch (err) {
                if (!cancelled) {
                    setError(err instanceof Error ? err.message : "Failed to load admin data");
                }
            } finally {
                if (!cancelled) {
                    setLoading(false);
                }
            }
        }

        void loadInitialState();
        return () => {
            cancelled = true;
        };
    }, [loadAdminData]);

    useEffect(() => {
        void loadProjectMembers(memberProjectId);
    }, [loadProjectMembers, memberProjectId]);

    function handleAddInitialProjectAccess() {
        setError(null);
        const numericProjectId = Number(projectId);
        if (!Number.isInteger(numericProjectId)) {
            setError("Select a project before adding initial access");
            return;
        }
        if (initialProjectAccess.some((access) => access.projectId === numericProjectId)) {
            setError("This project is already in the initial access list");
            return;
        }
        setInitialProjectAccess((current) => [
            ...current,
            {
                projectId: numericProjectId,
                permission,
            },
        ]);
    }

    function handleRemoveInitialProjectAccess(projectAccessId: number) {
        setInitialProjectAccess((current) =>
            current.filter((access) => access.projectId !== projectAccessId),
        );
    }

    async function handleCreateUser(event: FormEvent<HTMLFormElement>) {
        event.preventDefault();
        setSaving(true);
        setError(null);
        setSuccess(null);
        try {
            const result = await adminApi.createUser({
                email,
                firstName,
                lastName,
                projectAccess: initialProjectAccess,
            });
            setEmail("");
            setFirstName("");
            setLastName("");
            setInitialProjectAccess([]);
            const assignmentCount = result.projectAccess.length;
            setSuccess(`${result.user.email} ${result.created ? "created" : "updated"} with ${assignmentCount} project assignment${assignmentCount === 1 ? "" : "s"}`);
            await loadAdminData();
            const shouldReloadMembers = result.projectAccess.some(
                (access) => String(access.projectId) === memberProjectId,
            );
            if (shouldReloadMembers) {
                await loadProjectMembers(memberProjectId);
            } else if (!memberProjectId && result.projectAccess.length > 0) {
                setMemberProjectId(String(result.projectAccess[0].projectId));
            }
        } catch (err) {
            setError(err instanceof Error ? err.message : "Failed to create user");
        } finally {
            setSaving(false);
        }
    }

    async function handleAddProjectMember(event: FormEvent<HTMLFormElement>) {
        event.preventDefault();
        setError(null);
        setSuccess(null);
        const numericProjectId = Number(memberProjectId);
        const numericUserId = Number(addMemberUserId);
        if (!Number.isInteger(numericProjectId) || !Number.isInteger(numericUserId)) {
            setError("Select a project and user before adding access");
            return;
        }

        setMemberAction("add");
        try {
            const result = await adminApi.addProjectMember(numericProjectId, {
                userId: numericUserId,
                permission: addMemberPermission,
            });
            setSuccess(
                `${result.user.email} added with ${result.projectAccess.permission} access`,
            );
            setAddMemberUserId("");
            await loadProjectMembers(memberProjectId);
        } catch (err) {
            setError(err instanceof Error ? err.message : "Failed to add project member");
        } finally {
            setMemberAction(null);
        }
    }

    async function handleUpdateProjectMember(
        userId: number,
        nextPermission: AdminPermission,
    ) {
        setError(null);
        setSuccess(null);
        const numericProjectId = Number(memberProjectId);
        if (!Number.isInteger(numericProjectId)) {
            setError("Select a project before changing permissions");
            return;
        }

        setMemberAction(`update-${userId}`);
        try {
            const result = await adminApi.updateProjectMember(numericProjectId, userId, {
                permission: nextPermission,
            });
            setSuccess(
                `${result.user.email} changed to ${result.projectAccess.permission} access`,
            );
            await loadProjectMembers(memberProjectId);
        } catch (err) {
            setError(err instanceof Error ? err.message : "Failed to update permission");
        } finally {
            setMemberAction(null);
        }
    }

    async function handleRemoveProjectMember(userId: number) {
        setError(null);
        setSuccess(null);
        const numericProjectId = Number(memberProjectId);
        if (!Number.isInteger(numericProjectId)) {
            setError("Select a project before removing access");
            return;
        }

        setMemberAction(`remove-${userId}`);
        try {
            await adminApi.removeProjectMember(numericProjectId, userId);
            setSuccess("Project access removed");
            await loadProjectMembers(memberProjectId);
        } catch (err) {
            setError(err instanceof Error ? err.message : "Failed to remove project access");
        } finally {
            setMemberAction(null);
        }
    }

    return (
        <Box sx={{ minHeight: "100vh", bgcolor: "background.default" }}>
            <Box
                component="header"
                sx={{
                    bgcolor: "background.paper",
                    borderBottom: 1,
                    borderColor: "divider",
                    py: 2,
                }}
            >
                <Container maxWidth="xl">
                    <Stack
                        direction={{ xs: "column", md: "row" }}
                        spacing={2}
                        justifyContent="space-between"
                        alignItems={{ xs: "flex-start", md: "center" }}
                    >
                        <Box>
                            <Typography variant="h5">MDV Admin</Typography>
                            <Typography variant="body2" color="text.secondary">
                                Read-only deployment overview
                            </Typography>
                        </Box>
                        {session ? (
                            <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
                                <Chip
                                    color={session.isAdmin ? "success" : "default"}
                                    label={session.isAdmin ? "Admin session" : "No admin access"}
                                    size="small"
                                />
                                <Chip
                                    label={session.authEnabled ? "Auth enabled" : "Local dev auth"}
                                    size="small"
                                />
                                <Chip label={session.user.email} size="small" />
                            </Stack>
                        ) : null}
                    </Stack>
                </Container>
            </Box>

            <Container maxWidth="xl" sx={{ py: 3 }}>
                <Stack spacing={2}>
                    {error ? <Alert severity="error">{error}</Alert> : null}
                    {success ? <Alert severity="success">{success}</Alert> : null}
                    {loading ? (
                        <Stack alignItems="center" sx={{ py: 8 }}>
                            <CircularProgress />
                        </Stack>
                    ) : (
                        <Stack
                            direction={{ xs: "column", lg: "row" }}
                            spacing={2}
                            alignItems="flex-start"
                        >
                            <Paper
                                component="form"
                                onSubmit={handleCreateUser}
                                sx={{ p: 2, width: { xs: "100%", lg: 360 }, flexShrink: 0 }}
                            >
                                <Stack spacing={2}>
                                    <Typography variant="h6">Create Local User</Typography>
                                    <TextField
                                        label="Email"
                                        value={email}
                                        onChange={(event) => setEmail(event.target.value)}
                                        required
                                        size="small"
                                        type="email"
                                    />
                                    <Stack direction={{ xs: "column", sm: "row" }} spacing={1}>
                                        <TextField
                                            label="First name"
                                            value={firstName}
                                            onChange={(event) => setFirstName(event.target.value)}
                                            size="small"
                                            fullWidth
                                        />
                                        <TextField
                                            label="Last name"
                                            value={lastName}
                                            onChange={(event) => setLastName(event.target.value)}
                                            size="small"
                                            fullWidth
                                        />
                                    </Stack>
                                    <TextField
                                        select
                                        label="Project to add"
                                        value={projectId}
                                        onChange={(event) => setProjectId(event.target.value)}
                                        size="small"
                                        disabled={projects.length === 0}
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
                                        onChange={(event) => {
                                            const value = event.target.value;
                                            if (
                                                value === "view" ||
                                                value === "edit" ||
                                                value === "owner"
                                            ) {
                                                setPermission(value);
                                            }
                                        }}
                                        required
                                        size="small"
                                    >
                                        <MenuItem value="view">View</MenuItem>
                                        <MenuItem value="edit">Edit</MenuItem>
                                        <MenuItem value="owner">Owner</MenuItem>
                                    </TextField>
                                    {projects.length === 0 ? (
                                        <Alert severity="info">
                                            No projects are available. The user can still be
                                            created without project access.
                                        </Alert>
                                    ) : null}
                                    <Button
                                        type="button"
                                        variant="outlined"
                                        disabled={projects.length === 0}
                                        onClick={handleAddInitialProjectAccess}
                                    >
                                        Add initial access
                                    </Button>
                                    {initialProjectAccess.length === 0 ? (
                                        <Alert severity="info">
                                            This user will be created without project access.
                                        </Alert>
                                    ) : (
                                        <Stack spacing={1}>
                                            {initialProjectAccess.map((access) => {
                                                const project = projects.find(
                                                    (item) => item.id === access.projectId,
                                                );
                                                return (
                                                    <Stack
                                                        key={access.projectId}
                                                        direction="row"
                                                        spacing={1}
                                                        alignItems="center"
                                                        justifyContent="space-between"
                                                    >
                                                        <Typography variant="body2">
                                                            {project?.name ?? access.projectId} -{" "}
                                                            {access.permission}
                                                        </Typography>
                                                        <Button
                                                            type="button"
                                                            size="small"
                                                            color="error"
                                                            onClick={() =>
                                                                handleRemoveInitialProjectAccess(
                                                                    access.projectId,
                                                                )
                                                            }
                                                        >
                                                            Remove
                                                        </Button>
                                                    </Stack>
                                                );
                                            })}
                                        </Stack>
                                    )}
                                    <Button
                                        type="submit"
                                        variant="contained"
                                        disabled={saving}
                                    >
                                        {saving ? "Creating..." : "Create user"}
                                    </Button>
                                </Stack>
                            </Paper>

                            <Paper sx={{ p: 2, flex: 1, width: "100%" }}>
                                <Stack spacing={2}>
                                    <Stack
                                        direction="row"
                                        justifyContent="space-between"
                                        alignItems="center"
                                    >
                                        <Typography variant="h6">Projects</Typography>
                                        <Chip label={projects.length} size="small" />
                                    </Stack>
                                    {projects.length === 0 ? (
                                        <Alert severity="info">
                                            No projects found in the MDV database.
                                        </Alert>
                                    ) : (
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
                                                        <TableCell>
                                                            {project.accessLevel || "unknown"}
                                                        </TableCell>
                                                        <TableCell>
                                                            {project.isDeleted
                                                                ? "Deleted"
                                                                : "Active"}
                                                        </TableCell>
                                                    </TableRow>
                                                ))}
                                            </TableBody>
                                        </Table>
                                    )}
                                </Stack>
                            </Paper>

                            <Paper sx={{ p: 2, flex: 1, width: "100%" }}>
                                <Stack spacing={2}>
                                    <Stack
                                        direction="row"
                                        justifyContent="space-between"
                                        alignItems="center"
                                    >
                                        <Typography variant="h6">Users</Typography>
                                        <Chip label={users.length} size="small" />
                                    </Stack>
                                    {users.length === 0 ? (
                                        <Alert severity="info">
                                            No users found in the MDV database yet.
                                        </Alert>
                                    ) : (
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
                                                        <TableCell>
                                                            {user.isActive ? "Active" : "Inactive"}
                                                        </TableCell>
                                                        <TableCell>
                                                            {user.isAdmin ? "Yes" : "No"}
                                                        </TableCell>
                                                    </TableRow>
                                                ))}
                                            </TableBody>
                                        </Table>
                                    )}
                                </Stack>
                            </Paper>
                        </Stack>
                    )}
                    {!loading ? (
                        <Paper sx={{ p: 2 }}>
                            <Stack spacing={2}>
                                <Stack
                                    direction={{ xs: "column", md: "row" }}
                                    spacing={2}
                                    justifyContent="space-between"
                                    alignItems={{ xs: "stretch", md: "center" }}
                                >
                                    <Stack direction="row" spacing={1} alignItems="center">
                                        <Typography variant="h6">Project Members</Typography>
                                        <Chip label={projectMembers.length} size="small" />
                                    </Stack>
                                    <TextField
                                        select
                                        label="Project"
                                        value={memberProjectId}
                                        onChange={(event) =>
                                            setMemberProjectId(event.target.value)
                                        }
                                        size="small"
                                        disabled={projects.length === 0}
                                        sx={{ minWidth: { xs: "100%", md: 280 } }}
                                    >
                                        {projects.map((project) => (
                                            <MenuItem key={project.id} value={String(project.id)}>
                                                {project.name}
                                            </MenuItem>
                                        ))}
                                    </TextField>
                                </Stack>
                                {projects.length === 0 ? (
                                    <Alert severity="info">
                                        No projects found in the MDV database.
                                    </Alert>
                                ) : (
                                    <Paper
                                        component="form"
                                        variant="outlined"
                                        onSubmit={handleAddProjectMember}
                                        sx={{ p: 2 }}
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
                                                onChange={(event) =>
                                                    setAddMemberUserId(event.target.value)
                                                }
                                                size="small"
                                                disabled={availableProjectUsers.length === 0}
                                                sx={{ minWidth: { md: 260 } }}
                                            >
                                                {availableProjectUsers.map((user) => (
                                                    <MenuItem
                                                        key={user.id}
                                                        value={String(user.id)}
                                                    >
                                                        {user.email}
                                                    </MenuItem>
                                                ))}
                                            </TextField>
                                            <TextField
                                                select
                                                label="Permission"
                                                value={addMemberPermission}
                                                onChange={(event) => {
                                                    const value = event.target.value;
                                                    if (
                                                        value === "view" ||
                                                        value === "edit" ||
                                                        value === "owner"
                                                    ) {
                                                        setAddMemberPermission(value);
                                                    }
                                                }}
                                                size="small"
                                                sx={{ minWidth: { md: 140 } }}
                                            >
                                                <MenuItem value="view">View</MenuItem>
                                                <MenuItem value="edit">Edit</MenuItem>
                                                <MenuItem value="owner">Owner</MenuItem>
                                            </TextField>
                                            <Button
                                                type="submit"
                                                variant="outlined"
                                                disabled={
                                                    memberAction === "add" ||
                                                    !addMemberUserId ||
                                                    availableProjectUsers.length === 0
                                                }
                                            >
                                                {memberAction === "add" ? "Adding..." : "Add"}
                                            </Button>
                                        </Stack>
                                        {availableProjectUsers.length === 0 ? (
                                            <Alert severity="info" sx={{ mt: 1 }}>
                                                Every known user already has access to this
                                                project.
                                            </Alert>
                                        ) : null}
                                    </Paper>
                                )}
                                {projects.length === 0 ? null : membersLoading ? (
                                    <Stack alignItems="center" sx={{ py: 4 }}>
                                        <CircularProgress size={24} />
                                    </Stack>
                                ) : projectMembers.length === 0 ? (
                                    <Alert severity="info">
                                        No users are assigned to this project yet.
                                    </Alert>
                                ) : (
                                    <Table size="small">
                                        <TableHead>
                                            <TableRow>
                                                <TableCell>User</TableCell>
                                                <TableCell>Email</TableCell>
                                                <TableCell>Permission</TableCell>
                                                <TableCell>Read</TableCell>
                                                <TableCell>Write</TableCell>
                                                <TableCell>Owner</TableCell>
                                                <TableCell align="right">Actions</TableCell>
                                            </TableRow>
                                        </TableHead>
                                        <TableBody>
                                            {projectMembers.map((member) => (
                                                <TableRow
                                                    key={`${member.projectAccess.projectId}-${member.user.id}`}
                                                >
                                                    <TableCell>
                                                        {formatUserName(member.user)}
                                                    </TableCell>
                                                    <TableCell>{member.user.email}</TableCell>
                                                    <TableCell>
                                                        <TextField
                                                            select
                                                            value={
                                                                member.projectAccess.permission
                                                            }
                                                            onChange={(event) => {
                                                                const value =
                                                                    event.target.value;
                                                                if (
                                                                    value === "view" ||
                                                                    value === "edit" ||
                                                                    value === "owner"
                                                                ) {
                                                                    void handleUpdateProjectMember(
                                                                        member.user.id,
                                                                        value,
                                                                    );
                                                                }
                                                            }}
                                                            size="small"
                                                            disabled={
                                                                memberAction ===
                                                                `update-${member.user.id}`
                                                            }
                                                        >
                                                            <MenuItem value="view">View</MenuItem>
                                                            <MenuItem value="edit">Edit</MenuItem>
                                                            <MenuItem value="owner">
                                                                Owner
                                                            </MenuItem>
                                                        </TextField>
                                                    </TableCell>
                                                    <TableCell>
                                                        {member.projectAccess.canRead ? "Yes" : "No"}
                                                    </TableCell>
                                                    <TableCell>
                                                        {member.projectAccess.canWrite
                                                            ? "Yes"
                                                            : "No"}
                                                    </TableCell>
                                                    <TableCell>
                                                        {member.projectAccess.isOwner
                                                            ? "Yes"
                                                            : "No"}
                                                    </TableCell>
                                                    <TableCell align="right">
                                                        <Button
                                                            color="error"
                                                            size="small"
                                                            disabled={
                                                                memberAction ===
                                                                `remove-${member.user.id}`
                                                            }
                                                            onClick={() =>
                                                                void handleRemoveProjectMember(
                                                                    member.user.id,
                                                                )
                                                            }
                                                        >
                                                            {memberAction ===
                                                            `remove-${member.user.id}`
                                                                ? "Removing..."
                                                                : "Remove"}
                                                        </Button>
                                                    </TableCell>
                                                </TableRow>
                                            ))}
                                        </TableBody>
                                    </Table>
                                )}
                            </Stack>
                        </Paper>
                    ) : null}
                </Stack>
            </Container>
        </Box>
    );
}
