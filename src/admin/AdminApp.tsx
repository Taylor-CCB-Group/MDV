import { type FormEvent, useCallback, useEffect, useState } from "react";
import {
    Alert,
    Box,
    CircularProgress,
    Container,
    Stack,
    Typography,
} from "@mui/material";
import {
    AdminApiError,
    adminApi,
    type AdminPermission,
    type AdminProject,
    type AdminProjectMember,
    type AdminSession,
    type AdminUser,
} from "./api";
import {
    getRouteDescription,
    permissionFromValue,
    type AdminAccessState,
    type AdminRoute,
    type InitialProjectAccess,
    type PermissionSetter,
} from "./adminUtils";
import { AccessStatePanel } from "./components/AccessStatePanel";
import { AdminHeader } from "./components/AdminHeader";
import { ProfileView } from "./components/ProfileView";
import { ProjectsView } from "./components/ProjectsView";
import { UsersView } from "./components/UsersView";
import { buildApiUrl, buildDashboardUrl } from "@/utils/mdvRouting";

function getRouteFromHash(): AdminRoute {
    const hash = window.location.hash.replace(/^#\/?/, "");
    if (hash === "projects") return "projects";
    if (hash === "profile") return "profile";
    return "users";
}

function getAccessState(err: unknown): AdminAccessState {
    if (err instanceof AdminApiError) {
        if (err.status === 401) return "login_required";
        if (err.status === 403) return "admin_required";
    }
    return null;
}

function handleSignIn() {
    window.location.href = buildApiUrl("login");
}

function handleSignOut() {
    window.location.href = buildApiUrl("logout");
}

function handleGoToProjects() {
    window.location.href = buildDashboardUrl();
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
    const [syncingUsers, setSyncingUsers] = useState(false);
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

    const authEnabled = session?.authEnabled ?? false;
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

    function handlePermissionChange(value: string, setter: PermissionSetter) {
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

    async function handleSyncUsers() {
        setError(null);
        setSuccess(null);
        setSyncingUsers(true);
        try {
            const result = await adminApi.syncUsers();
            await loadAdminData();
            if (memberProjectId) {
                await loadProjectMembers(memberProjectId);
            }
            const userDelta = result.usersAfter - result.usersBefore;
            const adminDelta = result.adminsAfter - result.adminsBefore;
            setSuccess(`${result.message} Users: ${result.usersBefore} -> ${result.usersAfter} (${userDelta >= 0 ? "+" : ""}${userDelta}). Admins: ${result.adminsBefore} -> ${result.adminsAfter} (${adminDelta >= 0 ? "+" : ""}${adminDelta}).`);
        } catch (err) {
            setError(err instanceof Error ? err.message : "Failed to sync users from Auth0");
        } finally {
            setSyncingUsers(false);
        }
    }

    function navigateAdmin(route: AdminRoute) {
        window.location.hash = `/${route}`;
        setActiveRoute(route);
    }

    function renderContent() {
        if (loading) {
            return (
                <Stack alignItems="center" sx={{ py: 8 }}>
                    <CircularProgress />
                </Stack>
            );
        }

        if (accessState) {
            return (
                <AccessStatePanel
                    accessState={accessState}
                    onSignIn={handleSignIn}
                    onSignOut={handleSignOut}
                />
            );
        }

        if (activeRoute === "users") {
            return (
                <UsersView
                    users={users}
                    projects={projects}
                    email={email}
                    firstName={firstName}
                    lastName={lastName}
                    projectId={projectId}
                    permission={permission}
                    initialProjectAccess={initialProjectAccess}
                    saving={saving}
                    onEmailChange={setEmail}
                    onFirstNameChange={setFirstName}
                    onLastNameChange={setLastName}
                    onProjectIdChange={setProjectId}
                    onPermissionChange={handlePermissionChange}
                    onPermissionStateChange={setPermission}
                    onAddInitialProjectAccess={handleAddInitialProjectAccess}
                    onRemoveInitialProjectAccess={handleRemoveInitialProjectAccess}
                    onCreateUser={(event) => void handleCreateUser(event)}
                />
            );
        }

        if (activeRoute === "projects") {
            return (
                <ProjectsView
                    projects={projects}
                    projectMembers={projectMembers}
                    selectedProject={selectedMemberProject}
                    availableProjectUsers={availableProjectUsers}
                    memberProjectId={memberProjectId}
                    addMemberUserId={addMemberUserId}
                    addMemberPermission={addMemberPermission}
                    memberAction={memberAction}
                    membersLoading={membersLoading}
                    onMemberProjectIdChange={setMemberProjectId}
                    onAddMemberUserIdChange={setAddMemberUserId}
                    onPermissionChange={handlePermissionChange}
                    onAddMemberPermissionChange={setAddMemberPermission}
                    onAddProjectMember={(event) => void handleAddProjectMember(event)}
                    onUpdateProjectMember={(userId, nextPermission) => void handleUpdateProjectMember(userId, nextPermission)}
                    onRemoveProjectMember={(userId) => void handleRemoveProjectMember(userId)}
                />
            );
        }

        return <ProfileView session={session} />;
    }

    return (
        <Box className="admin-showcase-page">
            <Box className="admin-main">
                <AdminHeader
                    activeRoute={activeRoute}
                    accessState={accessState}
                    authEnabled={authEnabled}
                    session={session}
                    syncingUsers={syncingUsers}
                    onNavigate={navigateAdmin}
                    onGoToProjects={handleGoToProjects}
                    onSignIn={handleSignIn}
                    onSignOut={handleSignOut}
                    onSyncUsers={() => void handleSyncUsers()}
                />

                <Container maxWidth="xl" className="admin-content" sx={{ py: 3 }}>
                    <Box className="admin-page-heading" sx={{ mb: 2.5 }}>
                        <Typography variant="body1" color="text.secondary">
                            {getRouteDescription(activeRoute)}
                        </Typography>
                    </Box>
                    <Stack spacing={2.5}>
                        {error ? (
                            <Alert severity="error" onClose={() => setError(null)}>
                                {error}
                            </Alert>
                        ) : null}
                        {success ? (
                            <Alert severity="success" onClose={() => setSuccess(null)}>
                                {success}
                            </Alert>
                        ) : null}
                        {renderContent()}
                    </Stack>
                </Container>
            </Box>
        </Box>
    );
}
