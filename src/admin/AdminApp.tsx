import { type FormEvent, useCallback, useEffect, useState } from "react";
import {
    Alert,
    Box,
    Button,
    Chip,
    CircularProgress,
    Container,
    Drawer,
    List,
    ListItemButton,
    ListItemIcon,
    ListItemText,
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
import { FolderKanban, UserPlus, UserRound, type LucideIcon } from "lucide-react";
import {
    AdminApiError,
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

type AdminAccessState = "login_required" | "admin_required" | null;
type AdminRoute = "users" | "projects" | "profile";

const adminRoutes: Array<{ id: AdminRoute; label: string; description: string; icon: LucideIcon }> = [
    {
        id: "users",
        label: "Manage Users",
        description: "Create deployment users and review account status",
        icon: UserPlus,
    },
    {
        id: "projects",
        label: "Manage Projects",
        description: "Manage project membership and permission levels",
        icon: FolderKanban,
    },
    {
        id: "profile",
        label: "Profile",
        description: "Review the current admin session",
        icon: UserRound,
    },
];

function getRouteFromHash(): AdminRoute {
    const hash = window.location.hash.replace(/^#\/?/, "");
    if (hash === "projects") return "projects";
    if (hash === "profile") return "profile";
    return "users";
}

function getRouteDescription(route: AdminRoute) {
    return adminRoutes.find((item) => item.id === route)?.description ?? adminRoutes[0].description;
}

function getRouteLabel(route: AdminRoute) {
    return adminRoutes.find((item) => item.id === route)?.label ?? adminRoutes[0].label;
}

function getAccessState(err: unknown): AdminAccessState {
    if (err instanceof AdminApiError) {
        if (err.status === 401) return "login_required";
        if (err.status === 403) return "admin_required";
    }
    return null;
}

function permissionFromValue(value: string): AdminPermission | null {
    if (value === "view" || value === "edit" || value === "owner") {
        return value;
    }
    return null;
}

export default function AdminApp() {
    const [activeRoute, setActiveRoute] = useState<AdminRoute>(() => getRouteFromHash());
    const [session, setSession] = useState<AdminSession | null>(null);
    const [users, setUsers] = useState<AdminUser[]>([]);
    const [projects, setProjects] = useState<AdminProject[]>([]);
    const [projectMembers, setProjectMembers] = useState<AdminProjectMember[]>([]);
    const [loading, setLoading] = useState(true);
    const [membersLoading, setMembersLoading] = useState(false);
    const [saving, setSaving] = useState(false);
    const [memberAction, setMemberAction] = useState<string | null>(null);
    const [accessState, setAccessState] = useState<AdminAccessState>(null);
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
    const selectedMemberProject = projects.find((project) => String(project.id) === memberProjectId);

    useEffect(() => {
        function handleHashChange() {
            setActiveRoute(getRouteFromHash());
        }

        if (!window.location.hash) {
            window.location.hash = "/users";
        }
        window.addEventListener("hashchange", handleHashChange);
        handleHashChange();
        return () => window.removeEventListener("hashchange", handleHashChange);
    }, []);

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
            setAccessState(null);
            try {
                await loadAdminData();
                if (cancelled) return;
            } catch (err) {
                if (!cancelled) {
                    const nextAccessState = getAccessState(err);
                    if (nextAccessState) {
                        setAccessState(nextAccessState);
                    } else {
                        setError(err instanceof Error ? err.message : "Failed to load admin data");
                    }
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

    function handlePermissionChange(value: string, setter: (next: AdminPermission) => void) {
        const parsed = permissionFromValue(value);
        if (parsed) setter(parsed);
    }

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
            setSuccess(`${result.user.email} added with ${result.projectAccess.permission} access`);
            setAddMemberUserId("");
            await loadProjectMembers(memberProjectId);
        } catch (err) {
            setError(err instanceof Error ? err.message : "Failed to add project member");
        } finally {
            setMemberAction(null);
        }
    }

    async function handleUpdateProjectMember(userId: number, nextPermission: AdminPermission) {
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
            setSuccess(`${result.user.email} changed to ${result.projectAccess.permission} access`);
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

    function navigateAdmin(route: AdminRoute) {
        window.location.hash = `/${route}`;
        setActiveRoute(route);
    }

    return (
        <Box className="admin-showcase-page">
            <Drawer
                variant="permanent"
                anchor="left"
                className="admin-drawer"
                PaperProps={{ className: "admin-drawer-paper" }}
            >
                <Box className="admin-drawer-content">
                    <Box className="admin-drawer-brand">
                        <Typography variant="h5">MDV Admin Portal</Typography>
                    </Box>
                    <List className="admin-drawer-list" disablePadding>
                        {adminRoutes.map((route) => {
                            const Icon = route.icon;
                            return (
                                <ListItemButton
                                    key={route.id}
                                    selected={activeRoute === route.id}
                                    onClick={() => navigateAdmin(route.id)}
                                >
                                    <ListItemIcon>
                                        <Icon size={18} strokeWidth={2.1} />
                                    </ListItemIcon>
                                    <ListItemText primary={route.label} />
                                </ListItemButton>
                            );
                        })}
                    </List>
                </Box>
            </Drawer>

            <Box className="admin-main">
                <Box component="header" className="admin-showcase-header">
                    <Container maxWidth="xl">
                        <Stack
                            direction={{ xs: "column", md: "row" }}
                            spacing={2}
                            justifyContent="space-between"
                            alignItems={{ xs: "flex-start", md: "center" }}
                        >
                            <Box>
                                <Typography variant="h4">{getRouteLabel(activeRoute)}</Typography>
                                <Typography variant="body2">
                                    {getRouteDescription(activeRoute)}
                                </Typography>
                            </Box>
                            {session ? (
                                <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
                                    <Chip
                                        label={session.authEnabled ? "Auth0" : "Local"}
                                        size="small"
                                    />
                                    <Chip
                                        color={session.isAdmin ? "success" : "default"}
                                        label={session.isAdmin ? "Admin" : "No admin access"}
                                        size="small"
                                    />
                                    <Chip className="admin-profile-chip" label={session.user.email} size="small" />
                                </Stack>
                            ) : null}
                        </Stack>
                    </Container>
                </Box>

                <Container maxWidth="xl" className="admin-content" sx={{ py: 3 }}>
                    <Stack spacing={2.5}>
                        {error ? <Alert severity="error">{error}</Alert> : null}
                        {success ? <Alert severity="success">{success}</Alert> : null}
                        {loading ? (
                            <Stack alignItems="center" sx={{ py: 8 }}>
                                <CircularProgress />
                            </Stack>
                        ) : accessState ? (
                            <Paper className="admin-showcase-panel" variant="outlined">
                                <Stack spacing={1.5}>
                                    <Typography variant="h6">
                                        {accessState === "login_required"
                                            ? "Sign in required"
                                            : "Admin access required"}
                                    </Typography>
                                    <Typography color="text.secondary">
                                        {accessState === "login_required"
                                            ? "You need to sign in before using the MDV Admin Portal."
                                            : "Your account is signed in, but it is not marked as an MDV admin for this deployment."}
                                    </Typography>
                                    <Alert severity={accessState === "login_required" ? "info" : "warning"}>
                                        {accessState === "login_required"
                                            ? "After signing in, reload this page to continue."
                                            : "Ask an existing admin or deployment owner to grant admin access before continuing."}
                                    </Alert>
                                </Stack>
                            </Paper>
                        ) : activeRoute === "users" ? (
                            <Stack spacing={2.5}>
                                <Paper
                                    component="form"
                                    onSubmit={handleCreateUser}
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
                                        <Stack direction={{ xs: "column", sm: "row" }} spacing={1}>
                                            <TextField
                                                select
                                                label="Project"
                                                value={projectId}
                                                onChange={(event) => setProjectId(event.target.value)}
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
                                                onChange={(event) => handlePermissionChange(event.target.value, setPermission)}
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
                                            onClick={handleAddInitialProjectAccess}
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
                                                                onClick={() => handleRemoveInitialProjectAccess(access.projectId)}
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
                        ) : activeRoute === "projects" ? (
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
                                                    {selectedMemberProject
                                                        ? `Managing permissions for ${selectedMemberProject.name}`
                                                        : "Select a project to manage access"}
                                                </Typography>
                                            </Box>
                                            <TextField
                                                select
                                                label="Project"
                                                value={memberProjectId}
                                                onChange={(event) => setMemberProjectId(event.target.value)}
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
                                            onSubmit={handleAddProjectMember}
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
                                                    onChange={(event) => setAddMemberUserId(event.target.value)}
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
                                                    onChange={(event) => handlePermissionChange(event.target.value, setAddMemberPermission)}
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
                                                            <TableCell>Permission</TableCell>
                                                            <TableCell>Read</TableCell>
                                                            <TableCell>Write</TableCell>
                                                            <TableCell>Owner</TableCell>
                                                            <TableCell align="right">Actions</TableCell>
                                                        </TableRow>
                                                    </TableHead>
                                                    <TableBody>
                                                        {projectMembers.map((member) => (
                                                            <TableRow key={`${member.projectAccess.projectId}-${member.user.id}`}>
                                                                <TableCell>
                                                                    <Typography variant="body2">{formatUserName(member.user)}</Typography>
                                                                    <Typography variant="caption">{member.user.email}</Typography>
                                                                </TableCell>
                                                                <TableCell>
                                                                    <TextField
                                                                        select
                                                                        value={member.projectAccess.permission}
                                                                        onChange={(event) => {
                                                                            const nextPermission = permissionFromValue(event.target.value);
                                                                            if (nextPermission) {
                                                                                void handleUpdateProjectMember(member.user.id, nextPermission);
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
                                                                <TableCell>{member.projectAccess.canRead ? "Yes" : "No"}</TableCell>
                                                                <TableCell>{member.projectAccess.canWrite ? "Yes" : "No"}</TableCell>
                                                                <TableCell>{member.projectAccess.isOwner ? "Yes" : "No"}</TableCell>
                                                                <TableCell align="right">
                                                                    <Button
                                                                        color="error"
                                                                        size="small"
                                                                        disabled={memberAction === `remove-${member.user.id}`}
                                                                        onClick={() => void handleRemoveProjectMember(member.user.id)}
                                                                    >
                                                                        {memberAction === `remove-${member.user.id}` ? "Removing..." : "Remove"}
                                                                    </Button>
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
                        ) : (
                            <Paper className="admin-showcase-panel admin-profile-panel" variant="outlined">
                                <Stack spacing={2.5}>
                                    <Box>
                                        <Typography variant="h5">Admin Profile</Typography>
                                        <Typography variant="body2">
                                            Mock profile details for the current admin session.
                                        </Typography>
                                    </Box>
                                    {session ? (
                                        <Stack spacing={1.5}>
                                            <Box className="admin-profile-row">
                                                <Typography variant="body2">Email</Typography>
                                                <Typography>{session.user.email}</Typography>
                                            </Box>
                                            <Box className="admin-profile-row">
                                                <Typography variant="body2">Role</Typography>
                                                <Typography>{session.isAdmin ? "Admin" : "No admin access"}</Typography>
                                            </Box>
                                            <Box className="admin-profile-row">
                                                <Typography variant="body2">Authentication</Typography>
                                                <Typography>{session.authEnabled ? "Auth0 enabled" : "Local dev mode"}</Typography>
                                            </Box>
                                            <Box className="admin-profile-row">
                                                <Typography variant="body2">Assigned workspace</Typography>
                                                <Typography>MDV deployment administrator</Typography>
                                            </Box>
                                            <Alert severity="info">
                                                Profile data is mocked for now and can be connected to real account metadata later.
                                            </Alert>
                                        </Stack>
                                    ) : null}
                                </Stack>
                            </Paper>
                        )}
                    </Stack>
                </Container>
            </Box>
        </Box>
    );
}
